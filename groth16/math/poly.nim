
#
# univariate polynomials over Fr
#
# constantine's implementation is "somewhat lacking", so we have to
# implement these ourselves...
#

import std/sequtils
import std/sugar

import constantine/math/arithmetic
#import constantine/math/io/io_fields
import constantine/named/properties_fields

import groth16/bn128
import groth16/math/domain
import groth16/math/ntt
import groth16/misc

#-------------------------------------------------------------------------------

type 
  Poly* = object
    coeffs* : seq[Fr[BN254_Snarks]]

#-------------------------------------------------------------------------------

func polyDegree*(P: Poly) : int =
  let xs = P.coeffs ; let n = xs.len
  var d : int = n-1
  while isZeroFr(xs[d]) and (d >= 0): d -= 1
  return d

func polyIsZero*(P: Poly) : bool = 
  let xs = P.coeffs ; let n = xs.len
  var b = true
  for i in 0..<n:
    if not isZeroFr(xs[i]):
      b = false
      break
  return b

func polyIsEqual*(P, Q: Poly) : bool = 
  let xs = P.coeffs ; let n = xs.len
  let ys = Q.coeffs ; let m = ys.len
  var b = true
  if n >= m:
    for i in 0..<m: ( if not isEqualFr(xs[i], ys[i]): ( b = false ; break ) )
    for i in m..<n: ( if not isZeroFr( xs[i]       ): ( b = false ; break ) )
  else:
    for i in 0..<n: ( if not isEqualFr(xs[i], ys[i]): ( b = false ; break ) )
    for i in n..<m: ( if not isZeroFr(        ys[i]): ( b = false ; break ) )
  return b

#-------------------------------------------------------------------------------

func polyEvalAt*(P: Poly, x0: Fr[BN254_Snarks]): Fr[BN254_Snarks] =
  let cs = P.coeffs ; let n = cs.len
  var y : Fr[BN254_Snarks] = zeroFr
  var r : Fr[BN254_Snarks] = oneFr
  if n > 0: y = cs[0]
  for i in 1..<n:
    r *= x0
    y += cs[i] * r
  return y

#-------------------------------------------------------------------------------

func polyNeg*(P: Poly) : Poly =
  let zs = map( P.coeffs , negFr )
  return Poly(coeffs: zs)

func polyAdd*(P, Q: Poly) : Poly =
  let xs = P.coeffs ; let n = xs.len
  let ys = Q.coeffs ; let m = ys.len
  var zs = newSeq[Fr[BN254_Snarks]](max(n,m))
  if n >= m:
    for i in 0..<m: zs[i] = ( xs[i] + ys[i] )
    for i in m..<n: zs[i] = ( xs[i]         )
  else:
    for i in 0..<n: zs[i] = ( xs[i] + ys[i] )
    for i in n..<m: zs[i] = (         ys[i] )
  return Poly(coeffs: zs)

func polySub*(P, Q: Poly) : Poly =
  let xs = P.coeffs ; let n = xs.len
  let ys = Q.coeffs ; let m = ys.len
  var zs = newSeq[Fr[BN254_Snarks]](max(n,m))
  if n >= m:
    for i in 0..<m: zs[i] = ( xs[i]  - ys[i] )
    for i in m..<n: zs[i] = ( xs[i]          )
  else:
    for i in 0..<n: zs[i] = ( xs[i]  - ys[i] )
    for i in n..<m: zs[i] = (   negFr( ys[i] ))
  return Poly(coeffs: zs)

#-------------------------------------------------------------------------------

func polyScale*(s: Fr, P: Poly): Poly =
  let zs = map( P.coeffs , proc (x: Fr[BN254_Snarks]): Fr[BN254_Snarks] = s*x )
  return Poly(coeffs: zs)

#-------------------------------------------------------------------------------

func polyMulNaive*(P, Q : Poly): Poly =
  let xs = P.coeffs ; let n1 = xs.len
  let ys = Q.coeffs ; let n2 = ys.len
  let N  = n1 + n2 - 1
  var zs = newSeq[Fr[BN254_Snarks]](N)
  for k in 0..<N:
    # 0 <= i <= min(k , n1-1)
    # 0 <= j <= min(k , n2-1)
    # k = i + j
    # 0 >= i = k - j >= k - min(k , n2-1)
    # 0 >= j = k - i >= k - min(k , n1-1)   
    let A : int = max( 0 , k - min(k , n2-1) )
    let B : int = min( k , n1-1 )
    zs[k] = zeroFr
    for i in A..B:
      let j = k-i
      zs[k] += xs[i] * ys[j]
  return Poly(coeffs: zs)

#-------------------------------------------------------------------------------

# multiply two polynomials using FFT
func polyMulFFT*(P, Q: Poly): Poly = 
  let n1 = P.coeffs.len
  let n2 = Q.coeffs.len

  let log2 : int    = max( ceilingLog2(n1) , ceilingLog2(n2) ) + 1
  let N    : int    = (1 shl log2)
  let D    : Domain = createDomain( N )

  let us = extendAndForwardNTT( P.coeffs, D )
  let vs = extendAndForwardNTT( Q.coeffs, D )
  let zs = collect( newSeq, (for i in 0..<N: us[i]*vs[i] ))
  let ws = inverseNTT( zs, D )

  return Poly(coeffs: ws)

#-------------------------------------------------------------------------------

# WARNING: this is using the naive implementation!
func polyMul*(P, Q : Poly): Poly =
  # return polyMulFFT(P, Q)   
  return polyMulNaive(P, Q)   

#-------------------------------------------------------------------------------

func `==`*(P, Q: Poly): bool = return polyIsEqual(P, Q)

func `+`*(P, Q: Poly): Poly  = return polyAdd(P, Q)
func `-`*(P, Q: Poly): Poly  = return polySub(P, Q)
func `*`*(P, Q: Poly): Poly  = return polyMul(P, Q)

func `*`*(s: Fr  , P: Poly): Poly  = return polyScale(s, P)
func `*`*(P: Poly, s: Fr  ): Poly  = return polyScale(s, P)

#-------------------------------------------------------------------------------

# the generalized vanishing polynomial `(a*x^N - b)`
func generalizedVanishingPoly*(N: int, a: Fr[BN254_Snarks], b: Fr[BN254_Snarks]): Poly =
  assert( N>=1 )
  var cs = newSeq[Fr[BN254_Snarks]]( N+1 )
  cs[0] = negFr(b)
  cs[N] = a
  return Poly(coeffs: cs)

# the vanishing polynomial `(x^N - 1)`
func vanishingPoly*(N: int): Poly = 
  return generalizedVanishingPoly(N, oneFr, oneFr)

func vanishingPoly*(D: Domain): Poly = 
  return vanishingPoly(D.domainSize)

#-------------------------------------------------------------------------------

type
  QuotRem*[T] = object
    quot* : T 
    rem*  : T 

# divide by the vanishing polynomial `(x^N - 1)`
# returns the quotient and remainder
func polyQuotRemByVanishing*(P: Poly, N: int): QuotRem[Poly] = 
  assert( N>=1 )
  let deg  : int     = polyDegree(P)
  let src  = P.coeffs
  var quot = newSeq[Fr[BN254_Snarks]]( max(1, deg - N + 1) )
  var rem  = newSeq[Fr[BN254_Snarks]]( N )

  if deg < N:
    rem = src
  
  else:

    # compute quotient
    for j in countdown(deg-N, 0):
      if j+N <= deg-N:
        quot[j] = src[j+N] + quot[j+N]
      else:
        quot[j] = src[j+N]

    # compute remainder
    for j in 0..<N:
      if j <= deg-N:
        rem[j] = src[j] + quot[j]
      else:
        rem[j] = src[j]

  return QuotRem[Poly]( quot:Poly(coeffs:quot), rem:Poly(coeffs:rem) )

# divide by the vanishing polynomial `(x^N - 1)`
func polyDivideByVanishing*(P: Poly, N: int): Poly = 
  let qr = polyQuotRemByVanishing(P, N)
  assert( polyIsZero(qr.rem) )
  return qr.quot

#-------------------------------------------------------------------------------

# Lagrange basis polynomials
func lagrangePoly*(D: Domain, k: int): Poly =
  let N             = D.domainSize
  let omMinusK = smallPowFr( D.invDomainGen , k )
  let invN     = invFr(intToFr(N))

  var cs  = newSeq[Fr[BN254_Snarks]]( N )
  if k == 0:
    for i in 0..<N: cs[i] = invN
  else:
    var s = invN
    for i in 0..<N: 
      cs[i] = s
      s *= omMinusK

  return Poly(coeffs: cs)

#---------------------------------------

# evaluate a Lagrange basis polynomial at a given point `zeta` (outside the domain)
func evalLagrangePolyAt*(D: Domain, k: int, zeta: Fr[BN254_Snarks]): Fr[BN254_Snarks] =
  let omegaK = smallPowFr(D.domainGen, k)
  let denom  = (zeta - omegaK)
  if bool(isZero(denom)):
    # we are inside the domain
    raise newException(AssertionDefect, "point should be outside the domain")
  else:
    # we are outside the domain
    return omegaK * (smallPowFr(zeta, D.domainSize) - oneFr) * D.invDomainSize * invFr(denom)

#-------------------------------------------------------------------------------

# evaluates a polynomial on an FFT domain
func polyForwardNTT*(P: Poly, D: Domain): seq[Fr[BN254_Snarks]] =
  let n = P.coeffs.len
  assert( n <= D.domainSize , "the domain must be as least as big as the polynomial" )
  let src = P.coeffs
  return forwardNTT(src, D)

#---------------------------------------

# interpolates a polynomial on an FFT domain
func polyInverseNTT*(ys: seq[Fr[BN254_Snarks]], D: Domain): Poly =
  let n = ys.len
  assert( n == D.domainSize , "the domain must be same size as the input" )
  let tgt = inverseNTT(ys, D)
  return Poly(coeffs: tgt)

#-------------------------------------------------------------------------------

