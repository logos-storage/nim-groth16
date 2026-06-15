
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

  DeltaImages* = object
    imageA*  : seq[bool]          # image of the witness update space under A
    imageB*  : seq[bool]          # image of the witness update space under B
    imageAB* : seq[bool]          # union of the previous two

  # things we can compute at circuit setup time
  DynaSetupV1* = object
    pointsDeltaLZ* : seq[G1]      # the points `delta^-1 * L_i(tau) * Z(tau) * g1` where `Z(x) = x^N-1`
    wConvDeltaLZ*  : seq[G1]      # the convolution of `W` and `pointsDeltaLZ`

  # things we can compute from the partial witness
  DynaPreprocessV1* = object
    projA0*       : seq[G1]      # the points `delta^-1 * A0(tau) * L_i(tau) * g1` (well, "modulo Z(tau)")
    projB0*       : seq[G1]      # the same for B0

  DynaPreProofV1* = object
    partialProof*   : PartialProof       # linear part of the partial proof
    deltaImages*    : DeltaImages       # image of the witness delta under A and B (where can update happen)
    dynaSetup*      : DynaSetupV1        # circuit-setup time calculations
    dynaPreprocess* : DynaPreprocessV1   # "Dynark"-style preprocessing

  # this is used in the "cross-term"
  OnlyAB* = object
    valuesAz*    : seq[F]
    valuesBz*    : seq[F]

  # this one is used in the "projection-term"
  PartialAB* = object
    valuesAz*    : seq[F]
    valuesBz*    : seq[F]
    complImageA* : seq[bool]     # image of the complement of the partial witness under A
    complImageB* : seq[bool]     # same for B

#-------------------------------------------------------------------------------
