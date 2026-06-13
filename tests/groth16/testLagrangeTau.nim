
{.used.} 

import std/unittest

# import constantine/math/arithmetic
# import constantine/named/properties_fields

#import groth16/bn128/fields
#import groth16/bn128/curves
import groth16/bn128/arrays
import groth16/bn128/rnd
# import groth16/bn128/debug

import groth16/math/domain
import groth16/math/convert
#import groth16/math/ntt
#import groth16/math/group_fft

import groth16/dynamic/types
# import groth16/dynamic/setup

#-------------------------------------------------------------------------------

suite "Lagrange-Tau":

  let N   : int    = 128
  let D   : Domain = createDomain( N )
  let tau : F      = randFr()

  let powers_of_tau = computePowersOfTau(N , tau)
  let lagrange_v0   = powersToLagrange( powers_of_tau )
  let lagrange_v1   = computeLagrangeOfTauV1(D , tau )
  let lagrange_v2   = computeLagrangeOfTauV2(D , tau )

  test "converting from powers of tau == computing directly via evaluation formula":
    check isEqualG1Seq( lagrange_v0.elements, lagrange_v1.elements )

  test "converting from powers of tau == computing via scalar FFT":
    check isEqualG1Seq( lagrange_v0.elements, lagrange_v2.elements )

  test "computing directly via evaluation formulaa == computing via scalar FFT":
    check isEqualG1Seq( lagrange_v1.elements, lagrange_v2.elements )

#-------------------------------------------------------------------------------

suite "to/from Jordi":

  let N = 128
  let tau   = randFr()
  let delta = randFr()

  test "to/from Jordi":  
    check testJordiConversion(N, tau, delta)

#-------------------------------------------------------------------------------
