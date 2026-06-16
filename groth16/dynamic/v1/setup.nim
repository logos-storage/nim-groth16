
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
import groth16/math/convert
# import groth16/math/ntt

import groth16/zkey_types

import groth16/dynamic/shared
import groth16/dynamic/v1/types

#-------------------------------------------------------------------------------

#  let deltaImgA  = sparseMatrixImage( A , delta_mask )
#  let deltaImgB  = sparseMatrixImage( B , delta_mask )
#  let deltaImgAB = orBoolSeqs( deltaImgA , deltaImgB )

#-------------------------------------------------------------------------------

# the points `delta^{-1} * phi_ii(tau) * (tau^N - 1) * g1` as in the Dynark paper
#
# where
#
# > phi_ii = - sum_j phi_ij = - sum_j ( W[j-i]*L[i] + W[i-j]*L[j] )
# > phi_ij = W[j-i]*L[i] + W[i-j]*L[j] 
#
func calculateDiagPhiFFT*( wvec: seq[F], deltaLZ: seq[G1] ): seq[G1] =
  let N = wvec.len
  assert( N == deltaLZ.len )

  let sumW = sumSeqFr( wvec )

  var hs: seq[G1] = groupConvolution( wvec , deltaLZ )
  for i in 0..<N:
    hs[i] += (sumW ** deltaLZ[i])
    hs[i] =  negG1(hs[i])

  return hs

# NOTE: this is EXTREMELY SLOW
proc calculateDiagPhiNaive*( wvec: seq[F], deltaLZ: seq[G1], pool: Taskpool ): seq[G1] =
  let N = wvec.len
  assert( N == deltaLZ.len )
  
  let sumW = sumSeqFr( wvec )

  var diagPhi: seq[G1] = newSeq[G1]( N )  
  for i in 0..<N:
    var ws: seq[F] = newSeq[F]( N )
    for j in 0..<N:
      if i != j:
        ws[j] = wvec[ safeMod(i-j , N) ] 
      else:
        ws[j] = sumW
    diagPhi[i] = negG1( msmMultiThreadedG1( ws, deltaLZ, pool ) )

  return diagPhi

#-------------------------------------------------------------------------------

# does the setup from the ZKey (prover key)
proc dynaSetupV1FromZKey*(zkey: Zkey, pool: Taskpool ): DynaSetupV1 = 

  let N = zkey.header.domainSize
  let D = createDomain(N)

  var deltaZTau : seq[G1]      # the points `delta^-1 * (tau^N-1) * tau^i * g1`
  case zkey.header.flavour     

    of JensGroth:
      deltaZTau = zkey.pPoints.pointsH1

    of Snarkjs:
      echo "Jordi-style .zkey detected; converting points! (slow...)"
      withMeasureTime(true,"Jordi-to-Jens conversion"):
        deltaZTau = convertPointsFromJordi(D , zkey.pPoints.pointsH1)
  
  let sumW    = sumOfWVec( N )
  let wvec    = calculateWVec( D )
  let deltaLZ = inverseGroupFFT( deltaZTau , D )
  let conv    = groupConvolution( wvec , deltaLZ )

  var diagPhi  : seq[G1]
  var diagPhi2 : seq[G1]
  withMeasureTime(true,"phi diagonal took"):
    diagPhi = calculateDiagPhiFFT( wvec , deltaLZ )

  echo "agree = " & $isEqualG1Seq( diagPhi , diagPhi2 ) 

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
