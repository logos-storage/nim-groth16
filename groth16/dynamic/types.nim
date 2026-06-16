
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
