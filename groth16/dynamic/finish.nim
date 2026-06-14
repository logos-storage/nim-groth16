
#
# Finishing a Dynark-style proof
#
# We are implementing the cross-term differently from the paper:
# 
# - in the paper there are only two convolutions (4 field FFT-s), but 2 MSM-s
# - we have three convolutions (6 field FFT-s), but only 1 MSM. 
#
# In practice our version is significantly faster, as MSM-s are much more expensive
# than field FFTs (at least for practical sizes)
#

import std/options

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

import groth16/zkey_types

import groth16/math/domain
import groth16/math/ntt
import groth16/math/poly

import groth16/dynamic/types
import groth16/dynamic/setup
import groth16/dynamic/shared

#-------------------------------------------------------------------------------

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

  var Ahat  = forwardNTT( As  , D )
  var Bhat  = forwardNTT( Bs  , D )
  var ABhat = forwardNTT( ABs , D )

  inplaceMulByFFTofWVecBar( Ahat  )
  inplaceMulByFFTofWVecBar( Bhat  )
  inplaceMulByFFTofWVecBar( ABhat )

  let Aconv  = inverseNTT( Ahat  , D)
  let Bconv  = inverseNTT( Bhat  , D)
  let ABconv = inverseNTT( ABhat , D)

  let sumW = sumOfWVec( N )

  var output: seq[F] = newSeq[F]( N )
 
  for k in 0..<N:
    output[k] = As[k]*Bconv[k] + Bs[k]*Aconv[k] - ABconv[k] - ABs[k]*sumW

  return output

#-------------------------------------------------------------------------------

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
