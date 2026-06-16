
{.push raises:[].}

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
import groth16/dynamic/v2/types

#-------------------------------------------------------------------------------

# does the setup from the ZKey (prover key)
proc dynaSetupV2FromZKey*(zkey: Zkey, pool: Taskpool ): DynaSetupV2 = 

  let N = zkey.header.domainSize
  let D = createDomain(N)

  let deltaZTau = getDeltaZTauPows(zkey)
  
  let sumW    = sumOfWVec( N )
  let wvec    = calculateWVec( D )
  let deltaLZ = inverseGroupFFT( deltaZTau , D )
  let conv    = groupConvolution( wvec , deltaLZ )

  var diagPhi  : seq[G1]
  withMeasureTime(true,"phi diagonal took"):
    diagPhi = calculateDiagPhiFFT( wvec , deltaLZ )

  return DynaSetupV2( weightVec     : wvec    ,
                      pointsDeltaLZ : deltaLZ , 
                      wConvDeltaLZ  : conv    ,
                      diagPhiPoints : diagPhi )

#-------------------------------------------------------------------------------
