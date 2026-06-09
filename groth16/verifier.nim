
#
# Groth16 prover
#
# WARNING! 
# the points H in `.zkey` are *NOT* what normal people would think they are
# See <https://geometry.xyz/notebook/the-hidden-little-secret-in-snarkjs>
#

#[
import sugar
import constantine/math/config/curves 
import constantine/math/io/io_fields
import constantine/math/io/io_bigints
import ./zkey
]#

# import constantine/math/arithmetic except Fp, Fr
import constantine/math/io/io_extfields
import constantine/named/properties_fields
import constantine/math/extension_fields/towers

import groth16/bn128
import groth16/zkey_types

from groth16/prover import Proof

#-------------------------------------------------------------------------------
# the verifier
#

proc verifyProof* (vkey: VKey, prf: Proof): bool {.raises: [ValueError].} =

  if prf.curve != "bn128":
    raise newException(ValueError,
      "unsupported curve: expected \"bn128\", got \"" & prf.curve & "\"")

  if not isInSubgroupG1(prf.pi_a): return false
  if not isInSubgroupG2(prf.pi_b): return false
  if not isInSubgroupG1(prf.pi_c): return false

  var pubG1 : G1 = msmG1( prf.publicIO , vkey.vpoints.pointsIC )  

  let lhs   = pairing( negG1(prf.pi_a) , prf.pi_b )          # < -pi_a   , pi_b  >
  let rhs1  = vkey.spec.alphaBeta                            # < alpha  , beta  >
  let rhs2  = pairing( prf.pi_c , vkey.spec.delta2 )         # < pi_c   , delta >
  let rhs3  = pairing( pubG1    , vkey.spec.gamma2 )         # < sum... , gamma >

  var eq : Fp12[BN254_Snarks]
  eq =  lhs  
  eq *= rhs1
  eq *= rhs2
  eq *= rhs3

  return bool(isOne(eq))

#-------------------------------------------------------------------------------
