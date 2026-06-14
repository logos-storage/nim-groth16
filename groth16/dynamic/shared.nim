
import std/options

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/poly
#import groth16/math/convolution

import groth16/zkey_types

import groth16/dynamic/types

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

#-------------------------------------------------------------------------------

# computes the vectors A*z, B*z (but skips C*z)
func buildOnlyAB*( zkey: ZKey, pwitness: seq[Option[F]] ): OnlyAB =
  let hdr: GrothHeader = zkey.header
  let domSize = hdr.domainSize

  var valuesAz = newSeq[F](domSize)
  var valuesBz = newSeq[F](domSize)
 
  for entry in zkey.coeffs:
    case entry.matrix 

      of MatrixA: 
        if isSome(pwitness[entry.col]):
          valuesAz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()

      of MatrixB: 
        if isSome(pwitness[entry.col]):
          valuesBz[entry.row] += entry.coeff * pwitness[entry.col].unsafeGet()

      else: raise newException(AssertionDefect, "fatal error")

  return OnlyAB( valuesAz:valuesAz ,
                 valuesBz:valuesBz )

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
