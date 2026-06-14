
{.push raises:[].}

#import constantine/named/properties_fields

import groth16/bn128
import groth16/math/domain
#import groth16/math/ntt
#import groth16/math/poly

import groth16/partial/types
export types

#-------------------------------------------------------------------------------

type

  # note: already exported from partial/types, and Nim is stupid as usual
  # F* = Fr[BN254_Snarks]

  PowersOfTau* = object
    elements* : seq[G1]

  LagrangeOfTau* = object
    domain*   : Domain
    elements* : seq[G1]

#-------------------------------------------------------------------------------

type 

  # things we can compute at circuit setup time
  DynaSetupV1* = object
    weightVec*     : seq[F]       # the weights `W_k = 1/N/(omega^-k - 1)`
    pointsDeltaLZ* : seq[G1]      # the points `delta^-1 * L_i(tau) * Z(tau) * g1` where `Z(x) = x^N-1`
    wConvDeltaLZ*  : seq[G1]      # the convolution of `W` and `pointsDeltaLZ`

  # things we can compute from the partial witness
  DynaPreprocess* = object
    projA0*       : seq[G1]      # the points `delta^-1 * A0(tau) * L_i(tau) * g1` (well, "modulo Z(tau)")
    projB0*       : seq[G1]      # the same for B0

  DynaPreProof* = object
    partialProof*   : PartialProof
    dyanPreprocess* : DynaPreprocess

  OnlyAB* = object
    valuesAz*    : seq[F]
    valuesBz*    : seq[F]

  PartialAB* = object
    valuesAz*    : seq[F]
    valuesBz*    : seq[F]
    complImageA* : seq[bool]     # image of the complement of the partial witness under A
    complImageB* : seq[bool]     # same for B

#-------------------------------------------------------------------------------
