
{.push raises:[].}

import constantine/math/arithmetic
# import constantine/named/properties_fields

import taskpools

import groth16/bn128
import groth16/bn128/arrays
import groth16/misc

import groth16/math/domain
import groth16/math/group_fft
import groth16/math/poly
import groth16/math/convolution
# import groth16/math/convert
# import groth16/math/ntt

import groth16/zkey_types

import groth16/dynamic/shared
import groth16/dynamic/v3/types

#-------------------------------------------------------------------------------

# does the setup from the ZKey (prover key)
proc dynaSetupV3FromZKey*(zkey: Zkey, subgroupSize: int, pool: Taskpool ): DynaSetupV3 = 

  assert( isPowerOfTwo(subgroupSize) , "subgroup size should be a power of two!")

  let N  = zkey.header.domainSize
  let K  = subgroupSize
  let D  = createDomain(N)
  let sg = createSubgroup( D , K)

  let deltaZTau = getDeltaZTauPows(zkey)

  let sumW    = sumOfWVec( N )
  let wvec    = calculateWVec( D )
  let deltaLZ = inverseGroupFFT( deltaZTau , D )
  let wconv   = groupConvolution( wvec , deltaLZ )

  let miniPts = selectOnSubgroup( sg , calculateDiagPhiFFT1( wvec , deltaLZ , wconv ) ) 

  return DynaSetupV3( imageSubgroup  : sg         ,
                      weightVec      : wvec       ,
                      pointsDeltaLZ  : deltaLZ    , 
                      wConvDeltaLZ   : wconv      ,
                      miniDiagPoints : miniPts    )

 #-------------------------------------------------------------------------------
