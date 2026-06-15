
{.push raises:[].}

# import constantine/named/properties_fields

import groth16/bn128
import groth16/misc

import groth16/math/domain
import groth16/math/group_fft
import groth16/math/poly
import groth16/math/convolution
import groth16/math/convert
# import groth16/math/ntt

import groth16/zkey_types

import groth16/dynamic/types
import groth16/dynamic/shared

#-------------------------------------------------------------------------------

#  let deltaImgA  = sparseMatrixImage( A , delta_mask )
#  let deltaImgB  = sparseMatrixImage( B , delta_mask )
#  let deltaImgAB = orBoolSeqs( deltaImgA , deltaImgB )

#-------------------------------------------------------------------------------

# does the setup from the ZKey (prover key)
proc dynaSetupV1FromZKey*(zkey: Zkey): DynaSetupV1 = 

  let N = zkey.header.domainSize
  let D = createDomain(N)

  # assert( zkey.header.flavour == JensGroth , "DynaSetupV1 requires classic (quotient) flavour, not Jordi's one!" )

  var deltaZTau : seq[G1]                     # the points `delta^-1 * (tau^N-1) * tau^i * g1`
  case zkey.header.flavour     

    of JensGroth:
      deltaZTau = zkey.pPoints.pointsH1

    of Snarkjs:
      echo "Jordi-style .zkey detected; converting points! (slow...)"
      withMeasureTime(true,"Jordi-to-Jens conversion"):
        deltaZTau = convertPointsFromJordi(D , zkey.pPoints.pointsH1)
  
  let wvec    = calculateWVec( D )
  let deltaLZ = inverseGroupFFT( deltaZTau , D )
  let conv    = groupConvolution( wvec , deltaLZ )

  return DynaSetupV1( pointsDeltaLZ : deltaLZ , 
                      wConvDeltaLZ  : conv    )

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

  let wvec   = calculateWVec( D )
  let conv   = groupConvolution( wvec , deltaLZ )

  return DynaSetupV1( pointsDeltaLZ : deltaLZ , 
                      wConvDeltaLZ  : conv    )

#-------------------------------------------------------------------------------
