
#
# the goal of Dynark-style preprocessing is to compute the `2d` group elements
#
#  U_k := L_k(\tau) * A0(tau) * g1 
#  V_k := L_k(\tau) * B0(tau) * g1 
#
# where A0(x), B0(x) are the polynomials corresponding to the partial witness
#
# this is V1, the version closest to the Dynark paper
#

import std/options

import taskpools

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/misc

import groth16/zkey_types

import groth16/math/domain
import groth16/math/convolution
import groth16/math/poly

import groth16/partial/precalc

import groth16/dynamic/types
import groth16/dynamic/setup
import groth16/dynamic/shared
import groth16/dynamic/finish


#-------------------------------------------------------------------------------

func projectionElementsV1*(setup: DynaSetupV1, D: Domain, As: seq[F], complImage: seq[bool]): seq[G1] =
  let N = D.domainSize
  var Us: seq[G1] = newSeq[G1]( N )

  # # reverse indexed wvec: wvecBar[i] = wvec[-i]
  # let wvecBar: seq[F] = fftReverseVec( setup.weightVec )

  let fldN : F = intToFr( N ) 
  let sumW : F = sumOfWVec( N ) 

  let WBarStarA : seq[F]  = fieldConvolveWithWVecBar( D , As ) 
  let ALstarW   : seq[G1] = groupConvolveWithWVec( D , pointwiseScaleG1( As , setup.pointsDeltaLZ ) )

  # 2*d scalar multiplications + the group convolution above
  for k in 0..<N: 
 
    # we only compute for the _image of_ the complementer of the partial witness
    if complImage[k]:         
      let cf : F = WBarStarA[k] - As[k] * sumW
      Us[k] = ALstarW[k] - (As[k] ** setup.wConvDeltaLZ[k]) + (cf ** setup.pointsDeltaLZ[k])

  return selectTrues( complImage , Us )

#---------------------------------------

# for testing purposes
func simulateProjectionElementsV1*( D: Domain, tau: F, delta: F, As: seq[F], complImage: seq[bool] ): seq[G1] =

  let N = D.domainSize
 
  let setup = simulateDynaSetupV1( D, tau, delta )

  let deltaInv : F = invFr(delta)

  # # reverse indexed wvec: wvecBar[i] = wvec[-i]
  # let wvecBar: seq[F] = fftReverseVec( calculateWVec( D ) )

  # sum_i A[i]*L_i(tau)
  var asLTau: F = zeroFr
  for i in 0..<N:
    asLTau += As[i] * evalLagrangePolyAt( D, i, tau )

  # delta^-1 * sum A[i] * L_i(tau)
  var deltaAsLTau: F = deltaInv * asLTau 
 
  #
  # a note about the correction term
  #
  # so in the actual protocol, we simply calculate things modulo `(x^N - 1)`. 
  # This very useful because otherwise we would have to change the trusted setup ceremony...
  # however, to make this simulation compatible, we have to do the same here!
  #
  # fortunately, here we know explicitly that the remainder modulo `(x^N - 1)` 
  # in `U[k]` is `delta^-1 * As[k] * L_k(x)`
  #

  var Us: seq[G1] = newSeq[G1]( N )
  for k in 0..<N:
    if complImage[k]:
      let Lk_tau : F = evalLagrangePolyAt(D, k, tau)
      let y      : F = deltaAsLTau * Lk_tau
      let corr   : F = deltaInv * As[k] * Lk_tau
      Us[k] = (y - corr) ** gen1
 
  return selectTrues( complImage , Us )

#---------------------------------------

proc testProjectionElementsV1*(N: int, tau: F, delta: F ): bool =

  let D  = createDomain( N )
  let As = randFrSeq( N )

  # image of the changes mask, let's just compute everything, easier for testing
  var trues: seq[bool] = newSeq[bool](N) ; for i in 0..<N: trues[i] = true

  let setup = simulateDynaSetupV1( D, tau, delta )

  let simulated = simulateProjectionElementsV1(D, tau, delta, As, trues)
  let computed  = projectionElementsV1( setup, D, As, trues )

  return isEqualG1Seq( simulated , computed )

#---------------------------------------

func dynaPreprocessV1*(zkey: ZKey, setup: DynaSetupV1, partialAB: PartialAB): DynaPreprocessV1 =
  let N = zkey.header.domainSize
  let D = createDomain(N)

  # note: as we want to compute A0(x)*B'(x), the image masks should be the other ones!!
  let projA = projectionElementsV1( setup , D , partialAB.valuesAz , partialAB.complImageB )
  let projB = projectionElementsV1( setup , D , partialAB.valuesBz , partialAB.complImageA )
  
  return DynaPreprocessV1( projA0: projA, projB0: projB )

#-------------------------------------------------------------------------------

proc dynaPreProofV1*(zkey: ZKey, setup: DynaSetupV1, partialWitness: PartialWitness, pool: Taskpool, printTimings: bool): DynaPreProofV1 =

  let N = zkey.header.domainSize
  let D = createDomain(N)

  var partialAB : PartialAB
  withMeasureTime(printTimings,"build partial AB"):
    partialAB = buildPartialAB( zkey, partialWitness.values )
 
  let imageAB = orBoolSeqs( partialAB.complImageA , partialAB.complImageB )
  let deltaImages = DeltaImages( imageA  : partialAB.complImageA ,
                                 imageB  : partialAB.complImageB ,
                                 imageAB : imageAB               )

  var partialProof : PartialProof
  withMeasureTime(printTimings,"precomputing the linear terms"):
    partialProof = generatePartialProof( zkey, partialWitness, pool, false )

  var nonlin: G1
  withMeasureTime(printTimings,"precomputing the nonlinear term `A0(tau)*B0(tau)` "):
    let cs = crossTermCoeffs(D , partialAB.valuesAz , partialAB.valuesBz )
    nonlin = msmMultiThreadedG1( cs , setup.pointsDeltaLZ , pool )

  partialProof.partial_pi_c += nonlin

  var preprocess : DynaPreprocessV1
  withMeasureTime(printTimings,"precomputing the \"projection\" points"):
    preprocess = dynaPreprocessV1( zkey, setup, partialAB )

  return DynaPreProofV1( partialProof   : partialProof ,
                         deltaImages    : deltaImages  ,
                         dynaSetup      : setup        ,
                         dynaPreprocess : preprocess   )

#-------------------------------------------------------------------------------



