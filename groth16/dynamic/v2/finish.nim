
#
# Finishing a Dynark-style proof
#
# We are implementing the cross-term differently from the paper:
# 
# - in the paper there are only two convolutions (4 field FFT-s), but 2 MSM-s
# - we have three convolutions (6 field FFT-s), but only 1 MSM. 
#
# In practice our version is significantly faster, as MSM-s are much more expensive
# than field FFTs (at least for practical sizes)
#

import std/options

import taskpools

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/misc

import groth16/zkey_types
import groth16/files/witness

import groth16/math/domain
import groth16/math/ntt
import groth16/math/poly
import groth16/math/convolution

import groth16/prover/types
import groth16/prover/shared
import groth16/partial/finish

import groth16/dynamic/shared
import groth16/dynamic/v2/types
import groth16/dynamic/v2/setup

#-------------------------------------------------------------------------------

proc finishDynaProofWithMaskV2*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV2, mask: Mask, pool: Taskpool, printTimings: bool): Proof =

  let N = zkey.header.domainSize
  let M = zkey.header.nvars
  let D = createDomain( N )

  let partialMask = dynaPreProof.partialProof.partial_mask

  let wvec    = dynaPreProof.dynaSetup.weightVec 
  let wvecRev = fftReverseVec( wvec )

  var deltaAB: OnlyAB
  var witnessDelta: seq[Option[F]] = newSeq[Option[F]]( M )
  withMeasureTime(printTimings,"building deltaAz, deltaBz"):
    for j in 0..<M:
      if partialMask[j]:
        witnessDelta[j] = none(F)
      else:
        witnessDelta[j] = some(wtns.values[j])
    deltaAB = buildOnlyAB( zkey , witnessDelta )

  var proof: Proof
  withMeasureTime(printTimings,"finishing the linear terms"):
    proof = finishPartialProofWithMaskGeneric( false, zkey, wtns, dynaPreProof.partialProof, mask, pool, false )

  var nonlin: G1 = infG1

#[ 

  --- this only makes sense if we can restrict to a subgroup ---
  
  var cs: seq[F]
  withMeasureTime(printTimings,"computing the nonlinear \"cross\" term"):
    withMeasureTime(printTimings," - computing the cross coeffs took"):
      cs = crossTermCoeffs(D , deltaAB.valuesAz , deltaAB.valuesBz )
    withMeasureTime(printTimings," - computing the cross MSM took"):
      nonlin = msmMultiThreadedG1( cs , dynaPreProof.dynaSetup.pointsDeltaLZ , pool )
  
  echo " - nonzero coefficients in deltaA = " & $countNonZerosFr(deltaAB.valuesAz)
  echo " - nonzero coefficients in deltaB = " & $countNonZerosFr(deltaAB.valuesBz)
  echo " - nonzero coefficients in the cross-term = " & $countNonZerosFr(cs)

]#

  # following the Dynark paper
  withMeasureTime(printTimings,"computing the nonlinear \"cross\" term (Dynark-style)"):

    let idxs = trueIndices( dynaPreProof.deltaImages.imageAB ) 
    let K = idxs.len

    withMeasureTime(printTimings," - computing the diagonal part of the cross-term"):
      var cs: seq[F]  = newSeq[F ]( K )
      var ps: seq[G1] = newSeq[G1]( K )
      for (k,i) in idxs.pairs:
        cs[k] = deltaAB.valuesAz[i] * deltaAB.valuesBz[i] 
        ps[k] = dynaPreProof.dynaSetup.diagPhiPoints[i]
      nonlin += msmMultiThreadedG1( cs , ps , pool ) 

    withMeasureTime(printTimings," - computing the off-diagonal part of the cross-term"):
      let As = deltaAB.valuesAz 
      let Bs = deltaAB.valuesBz 
      let AstarW = fieldConvolution( As , wvecRev )
      let BstarW = fieldConvolution( Bs , wvecRev )     # TODO: make this sparse somehow??!

      var ds: seq[F] = newSeq[F]( K )
      for (k,i) in idxs.pairs:
        ds[k] = As[i] * BstarW[i] + Bs[i] * AstarW[i]    

      let ls = selectTrues( dynaPreProof.deltaImages.imageAB , dynaPreProof.dynaSetup.pointsDeltaLZ )
      nonlin += msmMultiThreadedG1( ds , ls , pool ) 

#[
      var ds: seq[F]  = newSeq[F ]( K )
      var qs: seq[G1] = newSeq[G1]( K )
      for (k,i) in idxs.pairs:
        var x: F = zeroFr
        for j in idxs: 
          if (i != j):
            let ab = deltaAB.valuesAz[i] * deltaAB.valuesBz[j]  +
                     deltaAB.valuesAz[j] * deltaAB.valuesBz[i]
            x += ab * wvec[ safeMod( j - i , N ) ]
        ds[k] = x
        qs[k] = dynaPreProof.dynaSetup.pointsDeltaLZ[i]
      nonlin += msmMultiThreadedG1( ds , qs , pool ) 
 ]#

  withMeasureTime(printTimings,"computing the nonlinear \"projection\" terms"):
    let zs = selectFalses( partialMask , wtns.values )
    nonlin += msmMultiThreadedG1( zs , dynaPreProof.dynaPreprocess.unified , pool)

  proof.pi_c += nonlin

  return proof

proc finishDynaProofV2*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV2 , pool: Taskpool, printTimings = false ): Proof =
  let mask = randomMask()
  return finishDynaProofWithMaskV2( zkey, wtns, dynaPreProof, mask, pool, printTimings )

#-------------------------------------------------------------------------------
