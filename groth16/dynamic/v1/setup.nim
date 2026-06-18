
{.push raises:[].}

# import constantine/named/properties_fields

import taskpools

import groth16/bn128
import groth16/misc
#import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/group_fft
import groth16/math/poly
import groth16/math/convolution
# import groth16/math/convert
# import groth16/math/ntt

import groth16/zkey_types

import groth16/dynamic/shared
import groth16/dynamic/v1/types

#-------------------------------------------------------------------------------

# does the setup from the ZKey (prover key)
proc dynaSetupV1FromZKey*(zkey: Zkey, pool: Taskpool ): DynaSetupV1 = 

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

  return DynaSetupV1( weightVec     : wvec    ,
                      pointsDeltaLZ : deltaLZ , 
                      wConvDeltaLZ  : conv    ,
                      diagPhiPoints : diagPhi )

#---------------------------------------

# simulates a setup (from the "toxic waste" values `tau` and `delta`), so that 
# we can test components of the system
func simulateDynaSetupV1*( D: Domain, tau: F, delta: F ): DynaSetupV1 = 

  let N = D.domainSize
  let ztau:      F = smallPowFr(tau,N) - oneFr
  let deltaZTau: F = ztau / delta

  # compute `delta^-1 * L_i(tau) * (tau^N - 1) ** g1 
  var deltaLZ: seq[G1] = newSeq[G1]( N ) 
  for i in 0..<N: 
    let y: F = deltaZTau * evalLagrangePolyAt( D, i, tau )    
    deltaLZ[i] = y ** gen1

  let wvec    = calculateWVec( D )
  let conv    = groupConvolution(    wvec , deltaLZ )
  let diagPhi = calculateDiagPhiFFT( wvec , deltaLZ )

  return DynaSetupV1( weightVec     : wvec    ,
                      pointsDeltaLZ : deltaLZ , 
                      wConvDeltaLZ  : conv    ,
                      diagPhiPoints : diagPhi )

#-------------------------------------------------------------------------------
