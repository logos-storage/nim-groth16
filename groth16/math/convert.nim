
#
# convert between powers and Lagrange bases representations
#

{.push raises:[].}

import std/sugar
#import std/sequtils

import constantine/math/arithmetic
#import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/ntt
import groth16/math/group_fft
import groth16/math/poly

#import groth16/zkey_types
import groth16/dynamic/types

#-------------------------------------------------------------------------------
# Lagrange basis version of powers-of-tau

# convert the power basis to Lagrange basis
func powersToLagrange*(inp: PowersOfTau): LagrangeOfTau = 
  let D  = createDomain(inp.elements.len)
  let hs = inverseGroupFFT(inp.elements , D)
  return LagrangeOfTau(domain: D, elements:hs)

#-------------------------------------------------------------------------------
# Convert the "H" points between Jordi's and the Groth16 paper version

#[

the following should hold, because they have both detree 2N-1, 
and they also take the same values on the 2N sized subgroup:

  L_{2i+1}(x)  =  Z(x) * L_i( eta^-1 * x )

from this fact we should be able to convert?

]#

func convertPointsFromJordi*( D: Domain, pointsJordi: seq[G1] ): seq[G1] =
  let N   : int    = D.domainSize
  let D2  : Domain = createDomain( 2*N )
  let eta : F      = D2.domainGen

  # do an FFT
  var arr = forwardGroupFFT( pointsJordi , D )

  # multiply by powers of eta + a constant -2
  var s : F = negFr(twoFr) 
  for i in 0..<N:
    arr[i] = s ** arr[i] 
    s *= eta

  return arr

func convertPointsToJordi*( D: Domain, pointsJens: seq[G1] ): seq[G1] =
  let N   : int    = D.domainSize
  let D2  : Domain = createDomain( 2*N )
  let eta    : F   = D2.domainGen
  let etaInv : F   = invFr(eta)

  # multiply by powers of eta + a constant -1/2
  var arr : seq[G1] = newSeq[G1]( N )
  var s : F = negFr(oneHalfFr) 
  for i in 0..<N:
    arr[i] = s ** pointsJens[i] 
    s *= etaInv

  # do an IFFT
  return inverseGroupFFT( arr , D )

proc testJordiConversion*(N: int, tau: F, delta: F): bool =

  let D  : Domain = createDomain(   N )
  let D2 : Domain = createDomain( 2*N )

  let deltaInv : F  = invFr(delta)
  let ztauG1   : G1 = (smallPowFr(tau,N) - oneFr) ** gen1     # (tau^N - 1)

  #
  # in the original paper, these are the curve points
  #   [ delta^-1 * tau^i * Z(tau) ] 
  #
  let pointsJens: seq[G1] = collect( newSeq , (for i in 0..<N: 
        (deltaInv * smallPowFr(tau,i)) ** ztauG1 ))

  #
  # in the Snarkjs implementation, these are the curve points
  #   [ delta^-1 * L_{2i+1} (tau) ]
  # where L_k are the Lagrange polynomials on the refined domain
  #
  let pointsJordi: seq[G1] = collect( newSeq , (for i in 0..<N: 
        (deltaInv * evalLagrangePolyAt(D2, 2*i+1, tau)) ** gen1 ))

  let jensFromJordi = convertPointsFromJordi( D , pointsJordi )
  let jordiFromJens = convertPointsToJordi(   D , pointsJens  )
  
  let ok1 = isEqualG1Seq( pointsJens  , jensFromJordi )
  let ok2 = isEqualG1Seq( pointsJordi , jordiFromJens )

  # echo "fromJordi = " & ($ok1)
  # echo "toJordi   = " & ($ok2)

  return (ok1 and ok2)

#-------------------------------------------------------------------------------
# for testing purposes, computing the basics from tau

func computePowersOfScalar*(N: int, tau: F): seq[F] =
  var ts: seq[F] = newSeq[F]( N )
  var s: F = oneFr
  for i in 0..<N:
    ts[i] = s
    s *= tau
  return ts

func computePowersOfTau*(N: int, tau: F): PowersOfTau =
  var hs: seq[G1] = newSeq[G1]( N )
  var s: F = oneFr
  for i in 0..<N:
    hs[i] = s ** gen1
    s *= tau
  return PowersOfTau(elements: hs)

# direct computation
func computeLagrangeOfTauV1*(D: Domain, tau: F): LagrangeOfTau =
  let N  = D.domainSize
  var hs: seq[G1] = newSeq[G1]( N )
  for k in 0..<N:
    hs[k] = evalLagrangePolyAt(D, k, tau) ** gen1
  return LagrangeOfTau(domain: D, elements:hs)

# via Fourier transform
func computeLagrangeOfTauV2*(D: Domain, tau: F): LagrangeOfTau =
  let N  = D.domainSize
  let ts = computePowersOfScalar(N, tau) 
  let ls = inverseNTT(ts, D)
  var hs: seq[G1] = newSeq[G1]( N )
  for i in 0..<N:
    hs[i] = ls[i] ** gen1
  return LagrangeOfTau(domain: D, elements:hs)

#-------------------------------------------------------------------------------
