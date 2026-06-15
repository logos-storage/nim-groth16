
{.push raises:[].}

import std/options

import constantine/named/properties_fields

import groth16/bn128
import groth16/prover/types

export types

#-------------------------------------------------------------------------------

type

  F* = Fr[BN254_Snarks]

  PartialWitness* = object
    values* : seq[Option[F]]

  PartialProof* = object
    partial_mask*  : seq[bool]
    partial_pi_a*  : G1          # = [alpha]_1 + sum z_j*[A_j]_1
    partial_rho*   : G1          # = [beta]_1  + sum z_j*[B_j]_1
    partial_pi_b*  : G2          # = [beta]_2  + sum z_j*[B_j]_2
    partial_pi_c*  : G1          # =             sum z_j*[K_j]_1

#-------------------------------------------------------------------------------

func makePartialWitness*(vals: seq[Option[Fr[BN254_Snarks]]]): PartialWitness =
  return PartialWitness(values: vals)

# true where it's filled
func partialWitnessMask*(pw: PartialWitness): seq[bool] = 
  let N = pw.values.len
  var bs: seq[bool] = newSeq[bool]( N )
  for i in 0..<N:
    bs[i] = isSome(pw.values[i])
  return bs 

#-------------------------------------------------------------------------------

