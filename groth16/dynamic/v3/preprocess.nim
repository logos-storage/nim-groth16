
#
# this is V3, we compute the projection terms differently
#

import std/options
import std/tables             

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
import groth16/math/matrix  

import groth16/partial/precalc

import groth16/dynamic/shared
import groth16/dynamic/permute
# import groth16/dynamic/v2/types as v2types
import groth16/dynamic/v2/preprocess
import groth16/dynamic/v3/types
import groth16/dynamic/v3/setup

#-------------------------------------------------------------------------------

func projectionElementsV3*(setup: DynaSetupV3, As: seq[F] ): seq[G1] =
  let sg = setup.imageSubgroup
  let D  = sg.bigDomain
  let N  = D.domainSize
  let K  = sg.smallDomain.domainSize
  let ell = N div K

  var Us: seq[G1] = newSeq[G1]( K )

  let fldN : F = intToFr( N ) 
  let sumW : F = sumOfWVec( N ) 

  let WBarStarA : seq[F]  = fieldConvolveWithWVecBarOnSubgroup( sg , As ) 
  let ALstarW   : seq[G1] = groupConvolveWithWVecOnSubgroup( sg , pointwiseScaleG1( As , setup.pointsDeltaLZ ) )

  # 2*d scalar multiplications + the group convolution above
  for k in 0..<K: 
 
    let ellk = ell*k
    let cf : F = WBarStarA[k] - As[ellk] * sumW
    Us[k] = ALstarW[k] - (As[ellk] ** setup.wConvDeltaLZ[ellk]) + (cf ** setup.pointsDeltaLZ[ellk])

  return Us

#-------------------------------------------------------------------------------

proc dynaPreprocessV3*(zkey: ZKey, setup: DynaSetupV3, partialAB: PartialAB, partial_mask: seq[bool], zdeltaMask: seq[bool] ): DynaPreprocessV3 =
  let sg  = setup.imageSubgroup
  let D   = sg.bigDomain
  let K   = sg.smallDomain.domainSize
  let N   = D.domainSize
  let ell = N div K

  var projA_mini: seq[G1]
  var projB_mini: seq[G1]
  withMeasureTime(true," + the \"projection\" points U_k, V_k"):
    projA_mini = projectionElementsV3( setup , partialAB.valuesAz )
    projB_mini = projectionElementsV3( setup , partialAB.valuesBz )

  let m = countTrues(zdeltaMask)
  var Xs: seq[G1] = newSeq[G1]( m )

  let points_C1 = compactifyPointsC1( zkey, partial_mask )

  withMeasureTime(true," + the \"unified\" points X_j"):

    let matABC = zkeyToSparseMatrices(zkey)
    let colsA = selectTrues( zdeltaMask , matABC.A.columns )
    let colsB = selectTrues( zdeltaMask , matABC.B.columns )

    for j in 0..<m:
      var h : G1 = points_C1[j]
      for (k,a) in colsA[j].pairs: h += (a ** projB_mini[k div ell])
      for (k,b) in colsB[j].pairs: h += (b ** projA_mini[k div ell])
      Xs[j] = h

  return DynaPreprocessV3( unified: Xs )

#-------------------------------------------------------------------------------

proc dynaPreProofV3*(zkey: ZKey, setup: DynaSetupV3, partialWitness: PartialWitness, pool: Taskpool, printTimings: bool): DynaPreProofV3 =

  let N = zkey.header.domainSize
  let D = createDomain(N)

  let partialMask = partialWitnessMask( partialWitness )
  let zdeltaMask  = notBoolSeq( partialMask )

  var partialAB : PartialAB
  withMeasureTime(printTimings,"build partial AB"):
    partialAB = buildPartialAB( zkey, partialWitness.values )
 
  let imageAB = orBoolSeqs( partialAB.complImageA , partialAB.complImageB )
  let deltaImages = DeltaImages( imageA  : partialAB.complImageA ,
                                 imageB  : partialAB.complImageB ,
                                 imageAB : imageAB               )

  # guess the smaller domain
  let k0    = countTrues(imageAB)
  let log2k = ceilingLog2(k0)
  let K     = (1 shl log2k)
  let imgSubgroup = createSubgroup( D , K )

  echo "size of the delta image = " & $countTrues(imageAB)
  # echo $trueIndices( imageAB )

  assert( liesInSubgroup( imgSubgroup , imageAB ) , "in V3, we expect the delta image to lie in the proper subgroup !!" )

  var partialProof : PartialProof
  withMeasureTime(printTimings,"precomputing the linear terms"):
    partialProof = generatePartialProof( zkey, partialWitness, pool, false )

  var nonlin: G1
  withMeasureTime(printTimings,"precomputing the nonlinear term `A0(tau)*B0(tau)` "):
    let cs = crossTermCoeffs(D , partialAB.valuesAz , partialAB.valuesBz )
    nonlin = msmMultiThreadedG1( cs , setup.pointsDeltaLZ , pool )

  partialProof.partial_pi_c += nonlin

  var preprocess : DynaPreprocessV3
  withMeasureTime(printTimings,"precomputing the \"projection\" and then the unified points"):
    preprocess = dynaPreprocessV3( zkey, setup, partialAB, partialMask, zdeltaMask )

  return DynaPreProofV3( partialProof   : partialProof ,
                         deltaImages    : deltaImages  ,
                         dynaSetup      : setup        ,
                         dynaPreprocess : preprocess   )

#-------------------------------------------------------------------------------



