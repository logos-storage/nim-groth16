
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
import groth16/dynamic/v2/types as v2types
import groth16/dynamic/v2/preprocess
import groth16/dynamic/v3/types
import groth16/dynamic/v3/setup

#-------------------------------------------------------------------------------

func convert_preprocess_V2_to_V3( v2: DynaPreprocessV2 ): DynaPreprocessV3 = 
  return DynaPreprocessV3( unified: v2.unified )

func convert_setup_V3_to_V2( v3: DynaSetupV3 ): DynaSetupV2 =
  let fakePhi: seq[G1] = newSeq[G1]()
  return DynaSetupV2( weightVec     : v3.weightVec     ,
                      pointsDeltaLZ : v3.pointsDeltaLZ ,
                      wConvDeltaLZ  : v3.wConvDeltaLZ  ,
                      diagPhiPoints : fakePhi          )

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

  var preprocessV3 : DynaPreprocessV3
  withMeasureTime(printTimings,"precomputing the \"projection\" and then the unified points"):
    let setupV2   = convert_setup_V3_to_V2( setup )
    preprocessV3  = convert_preprocess_V2_to_V3( dynaPreprocessV2( zkey, setupV2, partialAB, partialMask, zdeltaMask ))

  return DynaPreProofV3( partialProof   : partialProof ,
                         deltaImages    : deltaImages  ,
                         dynaSetup      : setup        ,
                         dynaPreprocess : preprocessV3 )

#-------------------------------------------------------------------------------



