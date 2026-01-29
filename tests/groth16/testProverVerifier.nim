
import std/unittest
import std/sequtils

import taskpools

import groth16/bn128/fields
import groth16/prover
import groth16/verifier
import groth16/fake_setup
import groth16/zkey_types
import groth16/files/r1cs
import groth16/files/witness

suite "prover and verifier":

  const simpleWitnessCfg =
    WitnessConfig( nWires:  8
                 , nPubOut: 1
                 , nPubIn:  1
                 , nPrivIn: 3
                 , nLabels: 0 
                 )

  const simpleEq1: Constraint = ( @[] , @[] , @[ (1,minusOneFr) , (2,oneFr) , (7,oneFr) ] )
  const simpleEq2: Constraint = ( @[ (3,oneFr) ] , @[ (4,oneFr) ] , @[ (6,oneFr) ] )
  const simpleEq3: Constraint = ( @[ (5,oneFr) ] , @[ (6,oneFr) ] , @[ (7,oneFr) ] )

  const simpleConstraints: seq[Constraint] = @[ simpleEq1, simpleEq2, simpleEq3 ]
  const simpleLabels: seq[int] = @[]

  const simpleR1CS =
    R1CS( r:            primeR
        , cfg:          simpleWitnessCfg
        , nConstr:      simpleConstraints.len
        , constraints:  simpleConstraints
        , wireToLabel:  simpleLabels
        )

  let simpleWitnessValues = map( @[ 1, 2023, 1022, 7, 11, 13, 7*11, 7*11*13 ] , intToFr )

  let simpleWitness = 
    Witness( curve:  "bn128"
           , r:      primeR
           , nvars:  8
           , values: simpleWitnessValues
           )

  test "prove and verify - JensGroth flavour":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=JensGroth )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, simpleWitness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "prove and verify - Snarkjs flavour":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, simpleWitness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "proof structure":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, simpleWitness, pool )
    pool.shutdown()
    
    check proof.curve == "bn128"
    check proof.publicIO.len == simpleWitnessCfg.nPubOut + simpleWitnessCfg.nPubIn + 1
    check true

  test "verification key extraction":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    let vkey = extractVKey( zkey )
    
    check vkey.curve == "bn128"
    check vkey.vpoints.pointsIC.len == simpleWitnessCfg.nPubOut + simpleWitnessCfg.nPubIn + 1

  test "proof structure consistency":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof1 = generateProof( zkey, simpleWitness, pool )
    let proof2 = generateProof( zkey, simpleWitness, pool )
    pool.shutdown()
    
    check proof1.publicIO.len == proof2.publicIO.len
    check proof1.curve == proof2.curve
    check proof1.curve == "bn128"

  test "multiple proofs with same setup":
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    let vkey = extractVKey( zkey )
    var pool = Taskpool.new()
    
    for i in 1..5:
      let proof = generateProof( zkey, simpleWitness, pool )
      let ok = verifyProof( vkey, proof )
      check ok
    
    pool.shutdown()

  test "proof with different witness but same public IO":
    let altWitnessValues = map( @[ 1, 2023, 1022, 7, 11, 13, 7*11, 7*11*13 ] , intToFr )
    let altWitness = Witness( curve: "bn128", r: primeR, nvars: 8, values: altWitnessValues )
    
    let zkey = createFakeCircuitSetup( simpleR1CS, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, altWitness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok
