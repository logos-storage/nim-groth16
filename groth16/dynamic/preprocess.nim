
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

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

import groth16/zkey_types

import groth16/math/domain
import groth16/math/convolution
import groth16/math/poly

import groth16/dynamic/types
import groth16/dynamic/setup
import groth16/dynamic/shared

#-------------------------------------------------------------------------------

# computes the vectors A*z, B*z (but skips C*z)
func buildPartialAB*( zkey: ZKey, pwitness: seq[Option[Fr[BN254_Snarks]]] ): PartialAB =
  let hdr: GrothHeader = zkey.header
  let domSize = hdr.domainSize

  var valuesAz = newSeq[Fr[BN254_Snarks]](domSize)
  var valuesBz = newSeq[Fr[BN254_Snarks]](domSize)
 
  # we also compute the image of the complement of the partial witness under A and B
  var complImageA = newSeq[bool](domSize) 
  var complImageB = newSeq[bool](domSize) 
  for i in 0..<domSize: 
    complImageA[i] = false
    complImageB[i] = false

  for entry in zkey.coeffs:
    case entry.matrix 

      of MatrixA: 
        if isSome(pwitness[entry.col]):
          valuesAz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()
        else:
          complImageA[entry.row] = true

      of MatrixB: 
        if isSome(pwitness[entry.col]):
          valuesBz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()
        else:
          complImageB[entry.row] = true

      else: raise newException(AssertionDefect, "fatal error")

  return PartialAB( valuesAz:valuesAz,
                    valuesBz:valuesBz, 
                    complImageA:complImageA,
                    complImageB:complImageB )

#-------------------------------------------------------------------------------

func projectionElementsV1*(setup: DynaSetupV1, D: Domain, As: seq[F], complImage: seq[bool]): seq[G1] =
  let N = D.domainSize
  var Us: seq[G1] = newSeq[G1]( N )

  # reverse indexed wvec: wvecBar[i] = wvec[-i]
  let wvecBar: seq[F] = fftReverseVec( setup.weightVec )

  let fldN    : F = intToFr( N ) 
  let negSumW : F = (fldN - oneFr) / (fldN + fldN) 

  let WBarStarA : seq[F]  = fieldConvolution( As , wvecBar ) 
  let ALstarW   : seq[G1] = groupConvolution( setup.weightVec , pointwiseScaleG1( As , setup.pointsDeltaLZ ) )

  # 2*d scalar multiplications + the group convolution above
  for k in 0..<N: 
 
    # we only compute for the _image of_ the complementer of the partial witness
    if complImage[k]:         
      let cf : F = WBarStarA[k] - As[k] * negSumW
      Us[k] = ALstarW[k] - (As[k] ** setup.wConvDeltaLZ[k]) + (cf ** setup.pointsDeltaLZ[k])

  return Us

#---------------------------------------

# for testing purposes
func simulateProjectionElementsV1*( D: Domain, tau: F, delta: F, As: seq[F], complImage: seq[bool] ): seq[G1] =

  let N = D.domainSize
 
  let setup = simulateDynaSetupV1( D, tau, delta )

  let deltaInv : F = invFr(delta)

  # reverse indexed wvec: wvecBar[i] = wvec[-i]
  let wvecBar: seq[F] = fftReverseVec( setup.weightVec )

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
 
  return Us

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

func dynaPreprocessV1*(zkey: ZKey, setup: DynaSetupV1, partialWitness: PartialWitness): DynaPreprocess =
  let N = zkey.header.domainSize
  let D = createDomain(N)

  # let wtnsMask = partialWitnessMask(partialWitness)

  let partialAB = buildPartialAB( zkey, partialWitness.values )

  let projA = projectionElementsV1( setup , D , partialAB.valuesAz , partialAB.complImageA )
  let projB = projectionElementsV1( setup , D , partialAB.valuesBz , partialAB.complImageB )
  
  return DynaPreprocess( projA0: projA, projB0: projB )

#-------------------------------------------------------------------------------





