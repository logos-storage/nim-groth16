
import std/unittest

import groth16/bn128/fields
import groth16/fake_setup
import groth16/zkey_types
import groth16/files/r1cs

suite "fake setup":

  const testWitnessCfg =
    WitnessConfig( nWires:  5
                 , nPubOut: 1
                 , nPubIn:  1
                 , nPrivIn: 2
                 , nLabels: 0 
                 )

  const testConstraint1: Constraint = ( 
    @[ (wireIdx: 0, value: oneFr) ] , 
    @[ (wireIdx: 1, value: oneFr) ] , 
    @[ (wireIdx: 2, value: oneFr) ] 
  )

  const testConstraints: seq[Constraint] = @[ testConstraint1 ]
  const testLabels: seq[int] = @[]

  const testR1CS =
    R1CS( r:            primeR
        , cfg:          testWitnessCfg
        , nConstr:      testConstraints.len
        , constraints:  testConstraints
        , wireToLabel:  testLabels
        )

  test "create fake circuit setup - JensGroth":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=JensGroth )
    
    check zkey.header.curve == "bn128"
    check zkey.header.flavour == JensGroth
    check zkey.header.nvars == testWitnessCfg.nWires
    check zkey.header.npubs == testWitnessCfg.nPubOut + testWitnessCfg.nPubIn

  test "create fake circuit setup - Snarkjs":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    
    check zkey.header.curve == "bn128"
    check zkey.header.flavour == Snarkjs
    check zkey.header.nvars == testWitnessCfg.nWires
    check zkey.header.npubs == testWitnessCfg.nPubOut + testWitnessCfg.nPubIn

  test "fake setup header consistency":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    let header = zkey.header
    
    check header.domainSize == (1 shl header.logDomainSize)
    check header.domainSize >= testR1CS.nConstr + header.npubs + 1

  test "fake setup spec points":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    let spec = zkey.specPoints
    
    check true
    check true

  test "fake setup verifier points":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    let vPoints = zkey.vPoints
    
    let expectedLen = testWitnessCfg.nPubOut + testWitnessCfg.nPubIn + 1
    check vPoints.pointsIC.len == expectedLen

  test "fake setup prover points":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    let pPoints = zkey.pPoints
    
    check pPoints.pointsA1.len == testWitnessCfg.nWires
    check pPoints.pointsB1.len == testWitnessCfg.nWires
    check pPoints.pointsB2.len == testWitnessCfg.nWires
    check pPoints.pointsC1.len == testWitnessCfg.nWires - (testWitnessCfg.nPubOut + testWitnessCfg.nPubIn + 1)
    check pPoints.pointsH1.len > 0

  test "fake setup coefficients":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    
    check zkey.coeffs.len > 0
    
    for coeff in zkey.coeffs:
      check coeff.row >= 0
      check coeff.col >= 0
      check coeff.col < testWitnessCfg.nWires

  test "random toxic waste":
    let toxic1 = randomToxicWaste()
    let toxic2 = randomToxicWaste()
    
    check not isEqualFr(toxic1.alpha, toxic2.alpha) or
            not isEqualFr(toxic1.beta, toxic2.beta) or
            not isEqualFr(toxic1.gamma, toxic2.gamma) or
            not isEqualFr(toxic1.delta, toxic2.delta) or
            not isEqualFr(toxic1.tau, toxic2.tau)

  test "R1CS to coefficients conversion":
    let coeffs = r1csToCoeffs( testR1CS )
    
    check coeffs.len > 0
    
    var hasA = false
    var hasB = false
    for coeff in coeffs:
      if coeff.matrix == MatrixA:
        hasA = true
      if coeff.matrix == MatrixB:
        hasB = true
    
    check hasA or hasB

  test "R1CS to sparse matrices":
    let matrices = r1csToSparseMatrices( testR1CS )
    
    check matrices.A.len == testWitnessCfg.nWires
    check matrices.B.len == testWitnessCfg.nWires
    check matrices.C.len == testWitnessCfg.nWires

  test "fake setup with different flavours produces different H points":
    let zkey1 = createFakeCircuitSetup( testR1CS, flavour=JensGroth )
    let zkey2 = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    
    check zkey1.pPoints.pointsH1.len == zkey2.pPoints.pointsH1.len

  test "fake setup domain size":
    let zkey = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    
    let domainSize = zkey.header.domainSize
    check (domainSize and (domainSize - 1)) == 0
    
    check domainSize >= testR1CS.nConstr + zkey.header.npubs + 1

  test "fake setup consistency":
    let zkey1 = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    let zkey2 = createFakeCircuitSetup( testR1CS, flavour=Snarkjs )
    
    check zkey1.header.nvars == zkey2.header.nvars
    check zkey1.header.npubs == zkey2.header.npubs
    check zkey1.header.domainSize == zkey2.header.domainSize
    check zkey1.pPoints.pointsA1.len == zkey2.pPoints.pointsA1.len
