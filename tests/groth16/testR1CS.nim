
import std/unittest

import groth16/bn128/fields
import groth16/files/r1cs

suite "R1CS operations":

  test "R1CS creation":
    let cfg = WitnessConfig(
      nWires: 8,
      nPubOut: 1,
      nPubIn: 1,
      nPrivIn: 3,
      nLabels: 0
    )
    
    let constraint1: Constraint = (
      A: @[(wireIdx: 1, value: minusOneFr), (wireIdx: 2, value: oneFr), (wireIdx: 7, value: oneFr)],
      B: @[],
      C: @[]
    )
    
    let constraint2: Constraint = (
      A: @[(wireIdx: 3, value: oneFr)],
      B: @[(wireIdx: 4, value: oneFr)],
      C: @[(wireIdx: 6, value: oneFr)]
    )
    
    let constraints = @[constraint1, constraint2]
    
    let r1cs = R1CS(
      r: primeR,
      cfg: cfg,
      nConstr: constraints.len,
      constraints: constraints,
      wireToLabel: @[]
    )
    
    check r1cs.cfg.nWires == 8
    check r1cs.nConstr == 2
    check r1cs.constraints.len == 2

  test "R1CS constraint structure":
    let constraint: Constraint = (
      A: @[(wireIdx: 0, value: oneFr), (wireIdx: 1, value: intToFr(2))],
      B: @[(wireIdx: 2, value: intToFr(3))],
      C: @[(wireIdx: 3, value: intToFr(4))]
    )
    
    check constraint.A.len == 2
    check constraint.B.len == 1
    check constraint.C.len == 1
    check isEqualFr(constraint.A[0].value, oneFr)
    check isEqualFr(constraint.A[1].value, intToFr(2))

  test "R1CS witness config":
    let cfg = WitnessConfig(
      nWires: 10,
      nPubOut: 2,
      nPubIn: 3,
      nPrivIn: 4,
      nLabels: 1
    )
    
    check cfg.nWires == 10
    check cfg.nPubOut == 2
    check cfg.nPubIn == 3
    check cfg.nPrivIn == 4
    check cfg.nLabels == 1

  test "R1CS term creation":
    let term1: Term = (wireIdx: 5, value: intToFr(42))
    let term2: Term = (wireIdx: 3, value: minusOneFr)
    
    check term1.wireIdx == 5
    check isEqualFr(term1.value, intToFr(42))
    check term2.wireIdx == 3
    check isEqualFr(term2.value, minusOneFr)

  test "R1CS linear combination":
    let linComb: LinComb = @[
      (wireIdx: 0, value: oneFr),
      (wireIdx: 1, value: intToFr(2)),
      (wireIdx: 2, value: intToFr(3))
    ]
    
    check linComb.len == 3
    check linComb[0].wireIdx == 0
    check isEqualFr(linComb[1].value, intToFr(2))

  test "R1CS constraint evaluation":
    let x = intToFr(5)
    let y = intToFr(7)
    let expected = x * y
    
    let constraint: Constraint = (
      A: @[(wireIdx: 0, value: x)],
      B: @[(wireIdx: 1, value: y)],
      C: @[(wireIdx: 2, value: expected)]
    )
    
    let product = constraint.A[0].value * constraint.B[0].value
    check isEqualFr(product, constraint.C[0].value)

  test "R1CS empty constraint":
    let emptyConstraint: Constraint = (
      A: @[],
      B: @[],
      C: @[]
    )
    
    check emptyConstraint.A.len == 0
    check emptyConstraint.B.len == 0
    check emptyConstraint.C.len == 0

  test "R1CS multiple constraints":
    let constraints = @[
      (A: @[(wireIdx: 0, value: oneFr)], B: @[(wireIdx: 1, value: oneFr)], C: @[(wireIdx: 2, value: oneFr)]),
      (A: @[(wireIdx: 2, value: oneFr)], B: @[(wireIdx: 3, value: oneFr)], C: @[(wireIdx: 4, value: oneFr)]),
      (A: @[(wireIdx: 4, value: oneFr)], B: @[(wireIdx: 5, value: oneFr)], C: @[(wireIdx: 6, value: oneFr)])
    ]
    
    check constraints.len == 3
    for i, constraint in constraints:
      check constraint.A.len > 0 or constraint.B.len > 0 or constraint.C.len > 0

  test "R1CS wire to label mapping":
    let wireToLabel = @[0, 1, 2, 3, 4]
    let r1cs = R1CS(
      r: primeR,
      cfg: WitnessConfig(nWires: 5, nPubOut: 0, nPubIn: 0, nPrivIn: 0, nLabels: 0),
      nConstr: 0,
      constraints: @[],
      wireToLabel: wireToLabel
    )
    
    check r1cs.wireToLabel.len == 5
    check r1cs.wireToLabel[0] == 0
    check r1cs.wireToLabel[4] == 4
