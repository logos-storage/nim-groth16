
#
# Finishing a Dynark-style proof
#
# This version more-or-less follows the Dynark paper
#

import std/options

import taskpools

#import constantine/math/arithmetic
#import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/misc

import groth16/zkey_types
import groth16/files/witness

import groth16/math/domain
import groth16/math/convolution
#import groth16/math/poly
#import groth16/math/ntt

import groth16/prover/types
import groth16/prover/shared
import groth16/partial/finish

import groth16/dynamic/shared
import groth16/dynamic/v1/types

#-------------------------------------------------------------------------------

proc finishDynaProofWithMaskV1*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV1, mask: Mask, pool: Taskpool, printTimings: bool): Proof =

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

  withMeasureTime(printTimings,"computing the nonlinear \"projection\" terms"):
    let comprAz = selectTrues( dynaPreProof.deltaImages.imageA , deltaAB.valuesAz )
    let comprBz = selectTrues( dynaPreProof.deltaImages.imageB , deltaAB.valuesBz )

    # echo "comprAz.len = " & $comprAz.len
    # echo "comprBz.len = " & $comprBz.len
    # echo "projBA.len = " & $dynaPreProof.dynaPreprocess.projA0.len
    # echo "projB0.len = " & $dynaPreProof.dynaPreprocess.projB0.len

    nonlin += msmMultiThreadedG1( comprAz , dynaPreProof.dynaPreprocess.projB0 , pool ) 
    nonlin += msmMultiThreadedG1( comprBz , dynaPreProof.dynaPreprocess.projA0 , pool ) 

  proof.pi_c += nonlin

  return proof

proc finishDynaProofV1*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV1 , pool: Taskpool, printTimings = false ): Proof =
  let mask = randomMask()
  return finishDynaProofWithMaskV1( zkey, wtns, dynaPreProof, mask, pool, printTimings )

#-------------------------------------------------------------------------------
