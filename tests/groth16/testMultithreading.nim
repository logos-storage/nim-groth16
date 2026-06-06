{.used.}

# Multi-threading correctness tests.
#
# Two complementary checks:
#
#   1. Trivial-mask determinism (r=s=0): proof is a pure deterministic function
#      of (zkey, witness), so sweeping the thread count must produce
#      byte-identical proof points. Catches races that produce *different but
#      still valid* proofs across configurations.
#
#   2. Random-mask end-to-end verify: proves with random masking (the
#      production code path) under varied (gc-mode, thread-count) and asserts
#      every resulting proof verifies. Random masks change the MSM coefficient
#      inputs, which exercises the data-dependent (non-constant-time) parts of
#      the MSM where coefficient-magnitude-driven races have historically
#      hidden — invisible under trivial-mask testing.

import std/unittest
import std/sequtils

import taskpools

import groth16/prover
import groth16/prover/groth16 as proverImpl
import groth16/verifier
import groth16/fake_setup
import groth16/zkey_types
import groth16/files/witness
import groth16/files/r1cs
import groth16/bn128/fields

#-------------------------------------------------------------------------------
# Same simple multiplication circuit testProver.nim uses: 7*11*13 + 1022 = 2023.
# Small but exercises the full prover path (4 MSMs + quotient computation).

const myWitnessCfg =
  WitnessConfig( nWires:  8
               , nPubOut: 1
               , nPubIn:  1
               , nPrivIn: 3
               , nLabels: 0
               )

const myEq1 : Constraint = ( @[] , @[] , @[ (1,minusOneFr) , (2,oneFr) , (7,oneFr) ] )
const myEq2 : Constraint = ( @[ (3,oneFr) ] , @[ (4,oneFr) ] , @[ (6,oneFr) ] )
const myEq3 : Constraint = ( @[ (5,oneFr) ] , @[ (6,oneFr) ] , @[ (7,oneFr) ] )

const myConstraints : seq[Constraint] = @[ myEq1, myEq2, myEq3 ]

const myR1CS =
  R1CS( r:            primeR
      , cfg:          myWitnessCfg
      , nConstr:      myConstraints.len
      , constraints:  myConstraints
      , wireToLabel:  @[]
      )

let myWitnessValues = map( @[ 1, 2023, 1022, 7, 11, 13, 7*11, 7*11*13 ] , intToFr )

let myWitness =
  Witness( curve:  "bn128"
         , r:      primeR
         , nvars:  8
         , values: myWitnessValues
         )

const ThreadCounts = [1, 2, 4, 8]

#-------------------------------------------------------------------------------

proc proveWithThreads(zkey: ZKey, witness: Witness, nThreads: int): Proof =
  var pool = Taskpool.new(numThreads = nThreads)
  result = generateProofWithTrivialMask( zkey, witness, pool, printTimings = false )
  pool.shutdown()

proc verifyWith(zkey: ZKey, proof: Proof): bool =
  let vkey = extractVKey(zkey)
  return verifyProof(vkey, proof)

#-------------------------------------------------------------------------------

suite "multithreading":

  test "repeated proofs on the same pool match (no per-call state leak)":
    # Reusing one pool across many proofs must not change the output: rules
    # out residual state in worker-local buffers between invocations.
    let zkey = createFakeCircuitSetup( myR1cs, flavour=Snarkjs )
    var pool = Taskpool.new(numThreads = 4)
    defer: pool.shutdown()
    let first = generateProofWithTrivialMask(zkey, myWitness, pool, false)
    for _ in 0 ..< 4:
      let again = generateProofWithTrivialMask(zkey, myWitness, pool, false)
      check isEqualProof(first, again)
  
  test "trivial-mask proof is deterministic across thread counts (JensGroth)":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=JensGroth )
    let reference = proveWithThreads(zkey, myWitness, ThreadCounts[0])
    check verifyWith(zkey, reference)
    for j in ThreadCounts[1..^1]:
      let proof = proveWithThreads(zkey, myWitness, j)
      check isEqualProof(reference, proof)
      check verifyWith(zkey, proof)

  test "trivial-mask proof is deterministic across thread counts (Snarkjs)":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=Snarkjs )
    let reference = proveWithThreads(zkey, myWitness, ThreadCounts[0])
    check verifyWith(zkey, reference)
    for j in ThreadCounts[1..^1]:
      let proof = proveWithThreads(zkey, myWitness, j)
      check isEqualProof(reference, proof)
      check verifyWith(zkey, proof)

  test "random-mask proofs verify across thread counts (Snarkjs)":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=Snarkjs )
    let vkey = extractVKey(zkey)
    for j in ThreadCounts:
      var pool = Taskpool.new(numThreads = j)
      defer: pool.shutdown()
      for _ in 0 ..< 100:
        let proof = generateProof(zkey, myWitness, pool, false)
        check verifyProof(vkey, proof)

  test "random-mask proofs verify across thread counts (JensGroth)":
    let zkey = createFakeCircuitSetup( myR1cs, flavour=JensGroth )
    let vkey = extractVKey(zkey)
    for j in ThreadCounts:
      var pool = Taskpool.new(numThreads = j)
      defer: pool.shutdown()
      for _ in 0 ..< 100:
        let proof = generateProof(zkey, myWitness, pool, false)
        check verifyProof(vkey, proof)


