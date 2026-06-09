
{.used.}

import std/unittest
import std/sequtils

import taskpools

import constantine/math/io/io_fields
import constantine/named/properties_fields

import groth16/prover
import groth16/verifier
import groth16/fake_setup
import groth16/zkey_types
import groth16/files/witness
import groth16/files/r1cs
import groth16/bn128/fields
import groth16/bn128/curves

#-------------------------------------------------------------------------------
# simple hand-crafted arithmetic circuit
#

const myWitnessCfg =
  WitnessConfig( nWires:  8       # including the special variable "1"
               , nPubOut: 1       # public output = input + a*b*c = 1022 + 7*11*13 = 2023
               , nPubIn:  1       # public input  = 1022
               , nPrivIn: 3       # private inputs: 7, 11, 13
               , nLabels: 0 
               )

# 2023 == 1022 + 7*3*11
const myEq1 : Constraint = ( @[] , @[] , @[ (1,minusOneFr) , (2,oneFr) , (7,oneFr) ] )

# 7*11 == 77
const myEq2 : Constraint = ( @[ (3,oneFr) ] , @[ (4,oneFr) ] , @[ (6,oneFr) ] )

# 77*13 == 1001
const myEq3 : Constraint = ( @[ (5,oneFr) ] , @[ (6,oneFr) ] , @[ (7,oneFr) ] )

const myConstraints : seq[Constraint] = @[ myEq1, myEq2, myEq3 ]

const myLabels : seq[int] = @[]

const myR1CS =
  R1CS( r:            primeR
      , cfg:          myWitnessCfg
      , nConstr:      myConstraints.len
      , constraints:  myConstraints
      , wireToLabel:  myLabels
      )

# the equation we want prove is `7*11*13 + 1022 == 2023`
let myWitnessValues = map( @[ 1, 2023, 1022, 7, 11, 13, 7*11, 7*11*13 ] , intToFr )
# wire indices:         ^^^^^^^         0      1    2  3   4   5    6      7

let myWitness = 
  Witness( curve:  "bn128"
         , r:      primeR
         , nvars:  8
         , values: myWitnessValues
         )

#-------------------------------------------------------------------------------

proc testProof(zkey: ZKey, witness: Witness): bool = 
  var pool = Taskpool.new()
  let proof = generateProof( zkey, witness, pool )
  let vkey  = extractVKey( zkey)
  let ok    = verifyProof( vkey, proof )
  pool.shutdown()
  return ok

suite "prover":

  test "prove & verify simple multiplication circuit, `JensGroth` flavour":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=JensGroth )
    check testProof( zkey, myWitness )

  test "prove & verify simple multiplication circuit, `Snarkjs` flavour":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=Snarkjs )
    check testProof( zkey, myWitness )

#-------------------------------------------------------------------------------
# malformed-input handling: verifyProof must return `false`, never crash.

# An on-curve but out-of-subgroup G2 point (same construction as in testCurve.nim).
const badPt2_x1 = fromHex(Fp[BN254_Snarks], "0x2")
const badPt2_xu = fromHex(Fp[BN254_Snarks], "0x0")
const badPt2_y1 = fromHex(Fp[BN254_Snarks], "0x181604d0560080401c08b557815482553e278257d98100d193a011c42782474d")
const badPt2_yu = fromHex(Fp[BN254_Snarks], "0x04f21f9d99cc25f694cf22ff70dc0ac4692e7a721b725dc454a217f04bd03e33")
let badPt2_x   = mkFp2(badPt2_x1, badPt2_xu)
let badPt2_y   = mkFp2(badPt2_y1, badPt2_yu)
let outOfSubgroupG2 = unsafeMkG2(badPt2_x, badPt2_y)

# A G1 value with coordinates (1,1), which does NOT satisfy y^2 = x^3 + 3.
let bogusG1 = unsafeMkG1(oneFp, oneFp)

proc buildGoodProof(): (VKey, Proof) =
  var pool = Taskpool.new()
  let zkey  = createFakeCircuitSetup(myR1cs, flavour=Snarkjs)
  let proof = generateProof(zkey, myWitness, pool)
  let vkey  = extractVKey(zkey)
  pool.shutdown()
  return (vkey, proof)

suite "verifier input validation":

  test "wrong curve string raises ValueError":
    let (vkey, proof) = buildGoodProof()
    var bad = proof
    bad.curve = "bls12-381"
    expect ValueError:
      discard verifyProof(vkey, bad)

  test "out-of-subgroup G2 pi_b returns false":
    let (vkey, proof) = buildGoodProof()
    var bad = proof
    bad.pi_b = outOfSubgroupG2
    check (not verifyProof(vkey, bad))

  test "off-curve G1 pi_a returns false":
    let (vkey, proof) = buildGoodProof()
    var bad = proof
    bad.pi_a = bogusG1
    check (not verifyProof(vkey, bad))

  test "off-curve G1 pi_c returns false":
    let (vkey, proof) = buildGoodProof()
    var bad = proof
    bad.pi_c = bogusG1
    check (not verifyProof(vkey, bad))

  test "known-good proof still verifies":
    let (vkey, proof) = buildGoodProof()
    check verifyProof(vkey, proof)

#-------------------------------------------------------------------------------
# prover input-validation: generateProof must raise ValueError (catchable in
# both --panics modes) for malformed inputs, not AssertionDefect.

suite "prover input validation":

  test "curve mismatch between zkey and witness raises ValueError":
    let zkey = createFakeCircuitSetup(myR1cs, flavour=Snarkjs)
    var badWtns = myWitness
    badWtns.curve = "bls12-381"
    var pool = Taskpool.new()
    expect ValueError:
      discard generateProof(zkey, badWtns, pool)
    pool.shutdown()

  test "wrong witness length raises ValueError":
    let zkey = createFakeCircuitSetup(myR1cs, flavour=Snarkjs)
    var badWtns = myWitness
    badWtns.values = badWtns.values[0 ..< badWtns.values.len - 1]   # drop one wire
    var pool = Taskpool.new()
    expect ValueError:
      discard generateProof(zkey, badWtns, pool)
    pool.shutdown()

#-------------------------------------------------------------------------------
