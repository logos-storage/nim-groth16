
#
# this is V2, we compute the projection terms differently
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
import groth16/dynamic/v2/types
import groth16/dynamic/v2/setup

#-------------------------------------------------------------------------------

func projectionElementsV2*(setup: DynaSetupV2, D: Domain, As: seq[F], complImage: seq[bool]): seq[G1] =
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

  return Us

#---------------------------------------

proc dynaPreprocessV2*(zkey: ZKey, setup: DynaSetupV2, partialAB: PartialAB, zdeltaMask: seq[bool] ): DynaPreprocessV2 =
  let N = zkey.header.domainSize
  let D = createDomain(N)

  var projA_full: seq[G1]
  var projB_full: seq[G1]
  withMeasureTime(true," + the \"projection\" points U_k, V_k"):
    # note: as we want to compute A0(x)*B'(x), the image masks should be the other ones!!
    projA_full = projectionElementsV2( setup , D , partialAB.valuesAz , partialAB.complImageB )
    projB_full = projectionElementsV2( setup , D , partialAB.valuesBz , partialAB.complImageA )

  let m = countTrues(zdeltaMask)
  var Xs: seq[G1] = newSeq[G1]( m )

  withMeasureTime(true," + the \"unified\" points X_j"):

    let matABC = zkeyToSparseMatrices(zkey)
    let colsA = selectTrues( zdeltaMask , matABC.A.columns )
    let colsB = selectTrues( zdeltaMask , matABC.B.columns )

    for j in 0..<m:
      var h : G1 = infG1
      for (k,a) in colsA[j].pairs: h += (a ** projB_full[k])
      for (k,b) in colsB[j].pairs: h += (b ** projA_full[k])
      Xs[j] = h

  return DynaPreprocessV2( unified: Xs )

#[

  # just tmp testing
  let M = zkey.header.nvars
  let matABC = zkeyToSparseMatrices(zkey)
  let dmask = notBoolSeq( partialWitnessMask(partialWitness) )
  let cntA = selectTrues( dmask , sparseMatrixColumnCounts(matABC.A) )
  let cntB = selectTrues( dmask , sparseMatrixColumnCounts(matABC.B) )
  # echo $cntA
  # echo $cntB
  echo $sumIntSeq(cntA)
  echo $sumIntSeq(cntB)
  echo $countTrues(dmask)
  for j in 0..<M: 
    if dmask[j]:
      for (i,x) in matABC.A.columns[j].pairs:
        echo toDecimalFr(x)
  for j in 0..<M: 
    if dmask[j]:
      for (i,x) in matABC.B.columns[j].pairs:
        echo toDecimalFr(x)
]#

#-------------------------------------------------------------------------------

proc dynaPreProofV2*(zkey: ZKey, setup: DynaSetupV2, partialWitness: PartialWitness, pool: Taskpool, printTimings: bool): DynaPreProofV2 =

  let N = zkey.header.domainSize
  let D = createDomain(N)

  let zdeltaMask = notBoolSeq( partialWitnessMask( partialWitness ) )

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

  var preprocess : DynaPreprocessV2
  withMeasureTime(printTimings,"precomputing the \"projection\" points"):
    preprocess = dynaPreprocessV2( zkey, setup, partialAB, zdeltaMask )

  return DynaPreProofV2( partialProof   : partialProof ,
                         deltaImages    : deltaImages  ,
                         dynaSetup      : setup        ,
                         dynaPreprocess : preprocess   )

#-------------------------------------------------------------------------------



