
import std/unittest

import taskpools

import groth16/bn128
import groth16/prover
import groth16/verifier
import groth16/fake_setup
import groth16/zkey_types
import groth16/files/r1cs
import groth16/files/witness
import groth16/math/poly
import groth16/math/domain
import groth16/math/ntt

suite "regression tests":

  test "edge case - single constraint circuit":
    const cfg = WitnessConfig( nWires: 3, nPubOut: 0, nPubIn: 1, nPrivIn: 1, nLabels: 0 )
    const constraint: Constraint = ( 
      @[ (wireIdx: 0, value: oneFr) ] , 
      @[ (wireIdx: 1, value: oneFr) ] , 
      @[ (wireIdx: 2, value: oneFr) ] 
    )
    const r1cs = R1CS( r: primeR, cfg: cfg, nConstr: 1, constraints: @[constraint], wireToLabel: @[] )
    let witness = Witness( curve: "bn128", r: primeR, nvars: 3, values: @[oneFr, intToFr(2), intToFr(2)] )
    
    let zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, witness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "edge case - empty public inputs":
    const cfg = WitnessConfig( nWires: 4, nPubOut: 0, nPubIn: 0, nPrivIn: 2, nLabels: 0 )
    const constraint: Constraint = ( 
      @[ (wireIdx: 1, value: oneFr) ] , 
      @[ (wireIdx: 2, value: oneFr) ] , 
      @[ (wireIdx: 3, value: oneFr) ] 
    )
    const r1cs = R1CS( r: primeR, cfg: cfg, nConstr: 1, constraints: @[constraint], wireToLabel: @[] )
    let witness = Witness( curve: "bn128", r: primeR, nvars: 4, values: @[oneFr, intToFr(2), intToFr(2), intToFr(4)] )
    
    let zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, witness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "edge case - large number of constraints":
    let cfg = WitnessConfig( nWires: 10, nPubOut: 1, nPubIn: 1, nPrivIn: 7, nLabels: 0 )
    var constraints: seq[Constraint] = @[]
    
    for i in 0..<20:
      constraints.add((
        @[ (wireIdx: i mod 10, value: oneFr) ],
        @[ (wireIdx: (i+1) mod 10, value: oneFr) ],
        @[ (wireIdx: (i+2) mod 10, value: oneFr) ]
      ))
    
    let witnessValues = @[oneFr, oneFr, oneFr, oneFr, oneFr, oneFr, oneFr, oneFr, oneFr, oneFr]
    let r1cs = R1CS( r: primeR, cfg: cfg, nConstr: constraints.len, constraints: constraints, wireToLabel: @[] )
    let witness = Witness( curve: "bn128", r: primeR, nvars: 10, values: witnessValues )
    
    let zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, witness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "edge case - polynomial with leading zeros":
    let P = Poly(coeffs: @[intToFr(1), intToFr(2), zeroFr, zeroFr])
    check polyDegree(P) == 1
    check not polyIsZero(P)

  test "edge case - domain size 2":
    let D = createDomain(2)
    check D.domainSize == 2
    check D.logDomainSize == 1
    
    let points = enumerateDomain(D)
    check points.len == 2
    check isEqualFr(points[0], oneFr)
    check isEqualFr(points[1], D.domainGen)

  test "edge case - very large domain":
    let D = createDomain(1024)
    check D.domainSize == 1024
    check D.logDomainSize == 10
    
    let points = enumerateDomain(D)
    check points.len == 1024
    let wrapAround = points[1023] * D.domainGen
    check isEqualFr(wrapAround, oneFr)

  test "edge case - field element at prime boundary":
    let minusOne = minusOneFr
    let one = oneFr
    let sum = minusOne + one
    check isZeroFr(sum)

  test "edge case - batch inversion with single element":
    let xs = @[intToFr(5)]
    let ys = batchInverseFr(xs)
    check ys.len == 1
    check isEqualFr(xs[0] * ys[0], oneFr)

  test "edge case - batch inversion with zero element":
    let xs = @[intToFr(1), intToFr(2), intToFr(3)]
    let ys = batchInverseFr(xs)
    for i in 0..<xs.len:
      check isEqualFr(xs[i] * ys[i], oneFr)

  test "regression - polynomial multiplication identity":
    let P = Poly(coeffs: @[intToFr(1), intToFr(2), intToFr(3)])
    let onePoly = Poly(coeffs: @[oneFr])
    let result = P * onePoly
    check polyIsEqual(result, P)

  test "regression - polynomial addition zero":
    let P = Poly(coeffs: @[intToFr(1), intToFr(2)])
    let zeroPoly = Poly(coeffs: @[zeroFr])
    let result = P + zeroPoly
    check polyIsEqual(result, P)

  test "regression - NTT round trip with two elements":
    let D = createDomain(2)
    let input = @[intToFr(1), intToFr(2)]
    let forward = forwardNTT(input, D)
    let inverse = inverseNTT(forward, D)
    check input.len == inverse.len
    check isEqualFr(input[0], inverse[0])
    check isEqualFr(input[1], inverse[1])

  test "regression - vanishing polynomial at domain points":
    let N = 8
    let Z = vanishingPoly(N)
    let D = createDomain(N)
    let points = enumerateDomain(D)
    
    for point in points:
      let eval = polyEvalAt(Z, point)
      check isZeroFr(eval)

  test "regression - Lagrange basis at domain points":
    let N = 4
    let D = createDomain(N)
    let points = enumerateDomain(D)
    
    for k in 0..<N:
      let Lk = lagrangePoly(D, k)
      for i in 0..<N:
        let eval = polyEvalAt(Lk, points[i])
        if i == k:
          check isEqualFr(eval, oneFr)
        else:
          check isZeroFr(eval)

  test "regression - proof generation with minimal circuit":
    const cfg = WitnessConfig( nWires: 2, nPubOut: 0, nPubIn: 0, nPrivIn: 1, nLabels: 0 )
    const constraint: Constraint = ( 
      @[ (wireIdx: 0, value: oneFr) ] , 
      @[ (wireIdx: 1, value: oneFr) ] , 
      @[ (wireIdx: 1, value: oneFr) ] 
    )
    const r1cs = R1CS( r: primeR, cfg: cfg, nConstr: 1, constraints: @[constraint], wireToLabel: @[] )
    let witness = Witness( curve: "bn128", r: primeR, nvars: 2, values: @[oneFr, oneFr] )
    
    let zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
    var pool = Taskpool.new()
    let proof = generateProof( zkey, witness, pool )
    let vkey = extractVKey( zkey )
    let ok = verifyProof( vkey, proof )
    pool.shutdown()
    check ok

  test "regression - multiple proof generations don't interfere":
    const cfg = WitnessConfig( nWires: 3, nPubOut: 0, nPubIn: 1, nPrivIn: 1, nLabels: 0 )
    const constraint: Constraint = ( 
      @[ (wireIdx: 0, value: oneFr) ] , 
      @[ (wireIdx: 1, value: oneFr) ] , 
      @[ (wireIdx: 2, value: oneFr) ] 
    )
    const r1cs = R1CS( r: primeR, cfg: cfg, nConstr: 1, constraints: @[constraint], wireToLabel: @[] )
    let witness = Witness( curve: "bn128", r: primeR, nvars: 3, values: @[oneFr, intToFr(2), intToFr(2)] )
    
    let zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
    let vkey = extractVKey( zkey )
    var pool = Taskpool.new()
    
    for i in 1..10:
      let proof = generateProof( zkey, witness, pool )
      let ok = verifyProof( vkey, proof )
      check ok
    
    pool.shutdown()
