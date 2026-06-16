
import std/options

import taskpools

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/misc

import groth16/math/domain
import groth16/math/ntt
import groth16/math/group_fft
import groth16/math/poly
import groth16/math/convolution
import groth16/math/convert

import groth16/zkey_types
import groth16/dynamic/types

#-------------------------------------------------------------------------------

#
# the points `delta^-1 * (tau^N-1) * tau^i * g1`
#
# - in Jens-syle setup, these are simply part of the prover key
# - in Jordi-style setup, you have to do an group IFFT
#
proc getDeltaZTauPows*(zkey: ZKey): seq[G1] = 
  var deltaZTau : seq[G1]      
  case zkey.header.flavour     

    of JensGroth:
      deltaZTau = zkey.pPoints.pointsH1

    of Snarkjs:
      echo "Jordi-style .zkey detected; converting points! (slow...)"
      withMeasureTime(true,"Jordi-to-Jens conversion"):
        let N = zkey.header.domainSize
        let D = createDomain( N )
        deltaZTau = convertPointsFromJordi(D , zkey.pPoints.pointsH1)
  
  return deltaZTau

#-------------------------------------------------------------------------------

# reverse indexing vecBar[i] = vec[-i] 
func fftReverseVec*[T]( vec: seq[T] ): seq[T] =
  let N = vec.len
  var vecBar: seq[T] = newSeq[T]( N )
  vecBar[0] = vec[0]
  for i in 1..<N:
    vecBar[N-i] = vec[i]
  return vecBar

#-------------------------------------------------------------------------------

# the k-th element `W_k` of the "weight vector"
func Wcoeff*( D: Domain, k0: int ): F = 
  let N = D.domainSize
  let k = safeMod( k0 , N )
  let invOmega : F = invFr( D.domainGen )
  return invFr( intToFr(N) * ( smallPowFr(invOmega,k) - oneFr ) )

# the "weight vector" from the Dynark paper
#
# these weights appear in the expansion of the
# product of Lagrange polynomials `L_i(x)L_k(x)`)
#
func calculateWVec*( D: Domain ): seq[F] = 
  let N = D.domainSize
  var wvec : seq[F] = newSeq[F]( N )
  let invN     : F = invFr( intToFr(N)  )
  let invOmega : F = invFr( D.domainGen )
  wvec[0] = zeroFr
  for i in 1..<N:
    wvec[i] = invN / ( smallPowFr(invOmega,i) - oneFr )
  return wvec

func calculateWVecBar*( D: Domain ): seq[F] = 
  fftReverseVec( calculateWVec(D) )

func sumOfWVec*( N: int ): F = 
  let fN : F = intToFr( N ) 
  return (oneFr - fN) / (fN + fN)

#-------------------------------------------------------------------------------

# pointwise multiply by `FFT[W]_k = (k + (1-N)/2) / N`
proc inplaceMulByFFTofWVec*( xs: var seq[F] ) = 
  let N    = xs.len
  let fN   = intToFr( N )
  let invN = invFr(  fN )
  let c    = divBy2Fr(oneFr - fN) * invN
  for k in 0..<N:
    let u = c + intToFr(k) * invN
    xs[k] *= u

# pointwise multiply by `FFT[Wbar]_k = Bar[FFT[W]]_k`
proc inplaceMulByFFTofWVecBar*( xs: var seq[F] ) = 
  let N    = xs.len
  let fN   = intToFr( N )
  let invN = invFr(  fN )
  let c    = divBy2Fr(oneFr - fN) * invN
  for k in 0..<N:
    let i = (if k==0: 0 else: N-k)
    let u = c + intToFr(i) * invN
    xs[k] *= u

#---------------------------------------

# pointwise multiply group elements by `FFT[W]_k = (k + (1-N)/2) / N`
proc inplaceScalarMulByFFTofWVec*( gs: var seq[G1] ) = 
  let N    = gs.len
  let fN   = intToFr( N )
  let invN = invFr(  fN )
  let c    = divBy2Fr(oneFr - fN) * invN
  for k in 0..<N:
    let u = c + intToFr(k) * invN
    gs[k] = u ** gs[k]

# pointwise multiply group elements by `FFT[Wbar]_k = Bar[FFT[W]]_k`
proc inplaceScalarMulByFFTofWVecBar*( gs: var seq[G1] ) = 
  let N    = gs.len
  let fN   = intToFr( N )
  let invN = invFr(  fN )
  let c    = divBy2Fr(oneFr - fN) * invN
  for k in 0..<N:
    let i = (if k==0: 0 else: N-k)
    let u = c + intToFr(i) * invN
    gs[k] = u ** gs[k]

#---------------------------------------

proc fieldConvolveWithWVec*( D: Domain, xs: seq[F] ): seq[F] =
  var xsHat = forwardNTT( xs , D )
  inplaceMulByFFTofWVec( xsHat )
  return inverseNTT( xsHat , D)

proc fieldConvolveWithWVecBar*( D: Domain, xs: seq[F] ): seq[F] =
  var xsHat = forwardNTT( xs , D )
  inplaceMulByFFTofWVecBar( xsHat )
  return inverseNTT( xsHat , D)

proc groupConvolveWithWVec*( D: Domain, gs: seq[G1] ): seq[G1] =
  var gsHat = forwardGroupFFT( gs , D )
  inplaceScalarMulByFFTofWVec( gsHat )
  return inverseGroupFFT( gsHat , D)

proc groupConvolveWithWVecBar*( D: Domain, gs: seq[G1] ): seq[G1] =
  var gsHat = forwardGroupFFT( gs , D )
  inplaceScalarMulByFFTofWVecBar( gsHat )
  return inverseGroupFFT( gsHat , D)

#-------------------------------------------------------------------------------

# computes the vectors A*z, B*z (but skips C*z)
func buildOnlyAB*( zkey: ZKey, pwitness: seq[Option[F]] ): OnlyAB =
  let hdr: GrothHeader = zkey.header
  let domSize = hdr.domainSize

  var valuesAz = newSeq[F](domSize)
  var valuesBz = newSeq[F](domSize)
 
  for entry in zkey.coeffs:
    if not isZeroFr(entry.coeff):
      case entry.matrix 
  
        of MatrixA: 
          if isSome(pwitness[entry.col]):
            valuesAz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()
  
        of MatrixB: 
          if isSome(pwitness[entry.col]):
            valuesBz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()
  
        else: raise newException(AssertionDefect, "fatal error")

  return OnlyAB( valuesAz: valuesAz ,
                 valuesBz: valuesBz )

#---------------------------------------

# computes the vectors A*z, B*z (but skips C*z), and also some image masks under A and B
func buildPartialAB*( zkey: ZKey, pwitness: seq[Option[F]] ): PartialAB =
  let hdr: GrothHeader = zkey.header
  let domSize = hdr.domainSize

  var valuesAz = newSeq[F](domSize)
  var valuesBz = newSeq[F](domSize)
 
  # we also compute the image of the complement of the partial witness under A and B
  var complImageA = newSeq[bool](domSize) 
  var complImageB = newSeq[bool](domSize) 
  for i in 0..<domSize: 
    complImageA[i] = false
    complImageB[i] = false

  for entry in zkey.coeffs:
    if not isZeroFr(entry.coeff):
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

  return PartialAB( valuesAz: valuesAz,
                    valuesBz: valuesBz, 
                    complImageA: complImageA,
                    complImageB: complImageB )

#-------------------------------------------------------------------------------
# the phi(x) polynomials and Lagrange product decomposition (for testing purposes)

# evaluates the `phi_{ik}(x)` polynomial at a point `tau`
# we assume that `i != k`
func evalPhiAt*(D: Domain, i: int, k: int, tau: F): F =
  assert( not (i == k) )
  return ( Wcoeff(D , k - i) * evalLagrangePolyAt(D , i , tau) +
           Wcoeff(D , i - k) * evalLagrangePolyAt(D , k , tau) ) 
           
# evaluates the `phi_{ii}(x)` polynomial (as defined in the Dynark parper) at a point `tau`
func evalPhiDiagonalAt*(D: Domain, i:int, tau: F): F =
  let N = D.domainSize
  var s: F = zeroFr
  for k in 0..<N:
    if not (k == i):
      s -= evalPhiAt(D, i, k, tau)
  return s

# the product L_i(tau)L_k(tau) via the key lemma
func evalLagrangeProductViaLemma*(D: Domain, i: int, k: int, tau: F): F =
  let N = D.domainSize
  let ztau = smallPowFr(tau, N) - oneFr
  if i == k:
    return ztau * evalPhiDiagonalAt(D, i, tau) + evalLagrangePolyAt(D, i, tau) 
  else:
    return ztau * evalPhiAt(D, i, k, tau)

#-------------------------------------------------------------------------------
# *** diagonal "phi_ii" points

# the points `delta^{-1} * phi_ii(tau) * (tau^N - 1) * g1` as in the Dynark paper
#
# where
#
# > phi_ii = - sum_j phi_ij = - sum_j ( W[j-i]*L[i] + W[i-j]*L[j] )
# > phi_ij = W[j-i]*L[i] + W[i-j]*L[j] 
#
func calculateDiagPhiFFT1*( wvec: seq[F], deltaLZ: seq[G1], wConvDeltaLZ: seq[G1] ): seq[G1] =
  let N = wvec.len
  assert( N == deltaLZ.len )

  let sumW = sumSeqFr( wvec )

  var hs: seq[G1] = wConvDeltaLZ
  for i in 0..<N:
    hs[i] += (sumW ** deltaLZ[i])
    hs[i] =  negG1(hs[i])

  return hs

func calculateDiagPhiFFT*( wvec: seq[F], deltaLZ: seq[G1] ): seq[G1] =
  return calculateDiagPhiFFT1( wvec, deltaLZ, groupConvolution(wvec , deltaLZ) )

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
# *** cross-term coefficients

#
# Computes the expansion `f(x)` where `f(x) = A(x)B(x) mod (x^N - 1)` in
# terms of `(x^N-1)*L_i(x)`
#
# The inputs `As` and `Bs` are the Lagrange-basis coefficients of A(x) and B(x).
# 
# Note: the remainder modulo `Z(x) := x^N-1` is `sum_k a[k]*b[k]*L_k(x)`
# 
func crossTermCoeffs*(D: Domain, As: seq[F], Bs: seq[F]) : seq[F] =

  let N = D.domainSize
  assert( N == As.len )
  assert( N == Bs.len )

  let ABs = pointwiseProdFr( As, Bs )

  let Aconv  = fieldConvolveWithWVecBar( D , As  )
  let Bconv  = fieldConvolveWithWVecBar( D , Bs  )
  let ABconv = fieldConvolveWithWVecBar( D , ABs )

  let sumW = sumOfWVec( N )

  var output: seq[F] = newSeq[F]( N )
 
  for k in 0..<N:
    output[k] = As[k]*Bconv[k] + Bs[k]*Aconv[k] - ABconv[k] - ABs[k]*sumW

  return output

#---------------------------------------

func crossTermCoeffsSubgroup*(wvec: seq[F], sg: Subgroup, As: seq[F], Bs: seq[F]) : seq[F] =

  let D = sg.smallDomain
  let N = D.domainSize
  assert( N == As.len )
  assert( N == Bs.len )

  let ABs = pointwiseProdFr( As, Bs )

  let wvecBar = selectOnSubgroup( sg , fftReverseVec(wvec) )
 
  let Aconv  = fieldConvolution( wvecBar , As  )
  let Bconv  = fieldConvolution( wvecBar , Bs  )
  let ABconv = fieldConvolution( wvecBar , ABs )

  let sumW = sumOfWVec( N )

  var output: seq[F] = newSeq[F]( N )
 
  for k in 0..<N:
    output[k] = As[k]*Bconv[k] + Bs[k]*Aconv[k]   # - ABconv[k] - ABs[k]*sumW

  return output

#---------------------------------------

proc testCrossTermCoeffs*(N : int): bool = 

  let D = createDomain( N )

  let As = randFrSeq(N)
  let Bs = randFrSeq(N)

  let tau  = randFr()
  let ztau = smallPowFr(tau,N) - oneFr                      # Z(tau) = tau^N - 1

  # reference `A(tau)*B(tau)`
  var Atau: F = zeroFr
  var Btau: F = zeroFr
  for i in 0..<N:
    Atau += As[i] * evalLagrangePolyAt( D , i , tau )       # A(x) = sum_i A_i * L_i(x)
    Btau += Bs[i] * evalLagrangePolyAt( D , i , tau )
  var reference = Atau * Btau

  # correction term (because of the modulo Z(x) behaviour of what we test)
  for k in 0..<N:
    reference -= As[k] * Bs[k] * evalLagrangePolyAt( D , k , tau )

  # the thing we want to test
  let coeffs = crossTermCoeffs( D , As , Bs )
  var smart: F = zeroFr
  for i in 0..<N:
    let lztau = zTau * evalLagrangePolyAt( D , i , tau ) 
    smart += coeffs[i] * lztau

  return (smart === reference)

#-------------------------------------------------------------------------------
