
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

import groth16/prover/types
import groth16/prover/shared
import groth16/partial/finish

import groth16/dynamic/types
import groth16/dynamic/setup
import groth16/dynamic/shared

#-------------------------------------------------------------------------------

#
# Computes the expansion `f(x)` where `f(x) = A(x)B(x) mod (x^N - 1)` in
# terms of `(x^N-1)*L_i(x)`
#
# The inputs `As` and `Bs` are the Lagrange-basis coefficients of A(x) and B(x).
# 
# Note: the remainder modulo `Z(x) := x^N-1` is `sum_k a[k]*b[k]*L_k(x)`
# 
func crossTermCoeffs*(D: Domain, As: seq[F], Bs: seq[F]) : seq[F] =

  let N = D.domainSize
  assert( N == As.len )
  assert( N == Bs.len )

  let ABs = pointwiseProdFr( As, Bs )

  let Aconv  = fieldConvolveWithWVecBar( D , As  )
  let Bconv  = fieldConvolveWithWVecBar( D , Bs  )
  let ABconv = fieldConvolveWithWVecBar( D , ABs )

  let sumW = sumOfWVec( N )

  var output: seq[F] = newSeq[F]( N )
 
  for k in 0..<N:
    output[k] = As[k]*Bconv[k] + Bs[k]*Aconv[k] - ABconv[k] - ABs[k]*sumW

  return output

#-------------------------------------------------------------------------------

proc testCrossTermCoeffs*(N : int): bool = 

  let D = createDomain( N )

  let As = randFrSeq(N)
  let Bs = randFrSeq(N)

  let tau  = randFr()
  let ztau = smallPowFr(tau,N) - oneFr                      # Z(tau) = tau^N - 1

  # reference `A(tau)*B(tau)`
  var Atau: F = zeroFr
  var Btau: F = zeroFr
  for i in 0..<N:
    Atau += As[i] * evalLagrangePolyAt( D , i , tau )       # A(x) = sum_i A_i * L_i(x)
    Btau += Bs[i] * evalLagrangePolyAt( D , i , tau )
  var reference = Atau * Btau

  # correction term (because of the modulo Z(x) behaviour of what we test)
  for k in 0..<N:
    reference -= As[k] * Bs[k] * evalLagrangePolyAt( D , k , tau )

  # the thing we want to test
  let coeffs = crossTermCoeffs( D , As , Bs )
  var smart: F = zeroFr
  for i in 0..<N:
    let lztau = zTau * evalLagrangePolyAt( D , i , tau ) 
    smart += coeffs[i] * lztau

  return (smart === reference)

#-------------------------------------------------------------------------------

proc finishDynaProofWithMaskV1*( zkey: ZKey, wtns: Witness, dynaPreProof: DynaPreProofV1, mask: Mask, pool: Taskpool, printTimings: bool): Proof =

  let N = zkey.header.domainSize
  let M = zkey.header.nvars
  let D = createDomain( N )

  let partialMask = dynaPreProof.partialProof.partial_mask

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

  var nonlin: G1
  var cs: seq[F]
  withMeasureTime(printTimings,"computing the nonlinear \"cross\" term"):
    withMeasureTime(printTimings," - computing the cross coeffs took"):
      cs = crossTermCoeffs(D , deltaAB.valuesAz , deltaAB.valuesBz )
    withMeasureTime(printTimings," - computing the cross MSM took"):
      nonlin = msmMultiThreadedG1( cs , dynaPreProof.dynaSetup.pointsDeltaLZ , pool )

  echo " - nonzero coefficients in deltaA = " & $countNonZerosFr(deltaAB.valuesAz)
  echo " - nonzero coefficients in deltaB = " & $countNonZerosFr(deltaAB.valuesBz)
  echo " - nonzero coefficients in the cross-term = " & $countNonZerosFr(cs)

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
