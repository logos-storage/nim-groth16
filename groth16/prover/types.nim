
{.push raises:[].}

import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/debug

#-------------------------------------------------------------------------------

type
  Mask* = object
    r*: Fr[BN254_Snarks]              # masking coefficients
    s*: Fr[BN254_Snarks]              # for zero knowledge

#-------------------------------------------------------------------------------
# a Groth16 proof

type
  Proof* = object
    publicIO* : seq[Fr[BN254_Snarks]]
    pi_a*     : G1
    pi_b*     : G2
    pi_c*     : G1
    curve*    : string

func isEqualProof*(prf1, prf2: Proof): bool =
  return (prf1.pi_a === prf2.pi_a) and
         (prf1.pi_b === prf2.pi_b) and
         (prf1.pi_c === prf2.pi_c)

proc debugPrintProof*(prf: Proof) =
  debugPrintG1( "pi_A:" , prf.pi_a )
  debugPrintG2( "pi_B:" , prf.pi_b )
  debugPrintG1( "pi_C:" , prf.pi_c )

#-------------------------------------------------------------------------------
# Az, Bz, Cz column vectors
# 

type
  ABC* = object
    valuesAz* : seq[Fr[BN254_Snarks]]
    valuesBz* : seq[Fr[BN254_Snarks]]
    valuesCz* : seq[Fr[BN254_Snarks]]

#-------------------------------------------------------------------------------
