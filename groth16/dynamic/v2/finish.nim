
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
import std/sequtils

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

type LinearPart = object
  unfinished_proof : Proof     # pi_A and pi_B are done, but pi_C is still mostly empty
  rho_coeffs       : seq[F]    # r * witness
  rho_points       : seq[G1]   # just pointsB1

# you can turn on/off the computation of the nonlinear term (useful for the "dynark" version)
proc finishLinearTerms( zkey: ZKey, wtns: Witness, partialProof: PartialProof, mask: Mask, pool: Taskpool, printTimings: bool): LinearPart =

  assert( zkey.header.curve == wtns.curve )

  let witness = wtns.values

  let hdr  : GrothHeader  = zkey.header
  let spec : SpecPoints   = zkey.specPoints
  let pts  : ProverPoints = zkey.pPoints     

  let nvars = hdr.nvars
  let npubs = hdr.npubs

  assert( nvars == witness.len , "wrong witness length" )

  let partial_mask = partialProof.partial_mask

  # remark: with the special variable "1" we actuall have (npub+1) public IO variables
  var pubIO = newSeq[Fr[BN254_Snarks]]( npubs + 1)
  for i in 0..npubs: pubIO[i] = witness[i]             

  # note: we have to ignore public inputs (plus 1 for the special first entry)
  let startIdx : int = npubs + 1
  var count    : int = countFalses(partial_mask)       # number of NEW witness elements 

  # compactify the stuff so we can call normal vector MSM
  var compact_witness:  seq[F]  = newSeq[F](count)
  var compact_pointsA1: seq[G1] = newSeq[G1](count)
  var compact_pointsB1: seq[G1] = newSeq[G1](count)
  var compact_pointsB2: seq[G2] = newSeq[G2](count)
  block:
    var j: int = 0;
    for i in 0..<nvars:
      if (not partial_mask[i]):
        compact_witness[j ] = witness[i]
        compact_pointsA1[j] = pts.pointsA1[i]
        compact_pointsB1[j] = pts.pointsB1[i]
        compact_pointsB2[j] = pts.pointsB2[i]
        j += 1

  # masking coeffs
  let r = mask.r
  let s = mask.s

  assert( witness.len == pts.pointsA1.len )
  assert( witness.len == pts.pointsB1.len )
  assert( witness.len == pts.pointsB2.len )
  assert( hdr.domainSize == pts.pointsH1.len )

  # echo "count = " & $count
  # for x in compact_witness:
  #   echo toDecimalFr(x)

  var pi_a : G1 = partialProof.partial_pi_a
  withMeasureTime(printTimings,"computing pi_A (G1 MSM)"):
    pi_a += r ** spec.delta1
    pi_a += msmMultiThreadedG1( compact_witness , compact_pointsA1, pool )

  # var rho : G1 = partialProof.partial_rho
  # withMeasureTime(printTimings,"computing rho (G1 MSM)"):
  #   rho += msmMultiThreadedG1( compact_witness , compact_pointsB1, pool )

  var pi_b : G2 = partialProof.partial_pi_b
  withMeasureTime(printTimings,"computing pi_B (G2 MSM)"):
    pi_b += s ** spec.delta2
    pi_b += msmMultiThreadedG2( compact_witness , compact_pointsB2, pool )

  var pi_c : G1 = partialProof.partial_pi_c
  pi_c += s ** pi_a + r ** partialProof.partial_rho   

  var rho_coeffs: seq[F] = newSeq[F]( count )
  withMeasureTime(printTimings,"rho coefficients"):
    for j in 0..<count:
      rho_coeffs[j] = r * compact_witness[j]

  let prf = Proof( curve:"bn128", publicIO:pubIO, pi_a:pi_a, pi_b:pi_b, pi_c:pi_c )

  return LinearPart( unfinished_proof : prf              ,
                     rho_coeffs       : rho_coeffs       ,
                     rho_points       : compact_pointsB1 )

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

  var linpart: LinearPart
  withMeasureTime(printTimings,"finishing the linear terms"):
    linpart = finishLinearTerms( zkey, wtns, dynaPreProof.partialProof, mask, pool, true ) # false )

  var proof  : Proof = linpart.unfinished_proof
  var nonlin : G1    = infG1

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

  var cross_cs: seq[F]
  var cross_ps: seq[G1]

  # following the Dynark paper
  withMeasureTime(printTimings,"computing the nonlinear \"cross\" term (Dynark-style)"):

    let idxs = trueIndices( dynaPreProof.deltaImages.imageAB ) 
    let K = idxs.len

    withMeasureTime(false," - computing the diagonal part of the cross-term"):
      var cs: seq[F]  = newSeq[F ]( K )
      var ps: seq[G1] = newSeq[G1]( K )
      for (k,i) in idxs.pairs:
        cs[k] = deltaAB.valuesAz[i] * deltaAB.valuesBz[i] 
        ps[k] = dynaPreProof.dynaSetup.diagPhiPoints[i]

      # nonlin += msmMultiThreadedG1( cs , ps , pool ) 

      cross_cs = cs
      cross_ps = ps

    withMeasureTime(false," - computing the off-diagonal part of the cross-term"):
      let As = deltaAB.valuesAz 
      let Bs = deltaAB.valuesBz 
      let AstarW = fieldConvolution( As , wvecRev )
      let BstarW = fieldConvolution( Bs , wvecRev )     # TODO: make this sparse somehow??!

      var ds: seq[F] = newSeq[F]( K )
      for (k,i) in idxs.pairs:
        ds[k] = As[i] * BstarW[i] + Bs[i] * AstarW[i]    

      let ls = selectTrues( dynaPreProof.deltaImages.imageAB , dynaPreProof.dynaSetup.pointsDeltaLZ )

      # nonlin += msmMultiThreadedG1( ds , ls , pool ) 

      cross_cs = concat( cross_cs , ds )
      cross_ps = concat( cross_ps , ls )

  withMeasureTime(printTimings,"computing all the remaining MSMs in one go"):
    var zs = selectFalses( partialMask , wtns.values )
    var ps = dynaPreProof.dynaPreprocess.unified

    # include the rho term too - one longer MSM is faster than two shorter
    zs = concat( zs , linpart.rho_coeffs , cross_cs )
    ps = concat( ps , linpart.rho_points , cross_ps )

    nonlin += msmMultiThreadedG1( zs , ps , pool)

  proof.pi_c += nonlin

  return proof

proc finishDynaProofV2*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV2 , pool: Taskpool, printTimings = false ): Proof =
  let mask = randomMask()
  return finishDynaProofWithMaskV2( zkey, wtns, dynaPreProof, mask, pool, printTimings )

#-------------------------------------------------------------------------------
