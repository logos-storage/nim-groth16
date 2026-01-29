
import std/unittest
import std/sequtils

import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/fields
import groth16/files/witness

suite "witness operations":

  test "witness creation":
    let values = @[oneFr, intToFr(2), intToFr(3), intToFr(4)]
    let witness = makeWitnessBN254(values)
    
    check witness.curve == "bn128"
    check witness.nvars == 4
    check witness.values.len == 4
    check isEqualFr(witness.values[0], oneFr)
    check isEqualFr(witness.values[1], intToFr(2))

  test "witness structure":
    let values = map(toSeq(1..10), intToFr)
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 10
    check witness.values.len == 10
    check true

  test "witness with zero values":
    let values = @[zeroFr, zeroFr, zeroFr]
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 3
    check witness.values.len == 3
    check isZeroFr(witness.values[0])
    check isZeroFr(witness.values[1])
    check isZeroFr(witness.values[2])

  test "witness with large values":
    let values = @[intToFr(1000000), intToFr(999999), intToFr(123456)]
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 3
    check isEqualFr(witness.values[0], intToFr(1000000))

  test "witness single value":
    let values = @[oneFr]
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 1
    check witness.values.len == 1
    check isEqualFr(witness.values[0], oneFr)

  test "witness empty (should fail or handle gracefully)":
    let values: seq[Fr[BN254_Snarks]] = @[]
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 0
    check witness.values.len == 0

  test "witness field consistency":
    let values = @[oneFr, intToFr(5), intToFr(10)]
    let witness = makeWitnessBN254(values)
    
    # All values should be valid field elements
    check witness.values.len == 3
    check isEqualFr(witness.values[0], oneFr)
    check isEqualFr(witness.values[1], intToFr(5))
    check isEqualFr(witness.values[2], intToFr(10))

  test "witness ordering":
    let publicOut = intToFr(100)
    let publicIn = intToFr(200)
    let privateIn = intToFr(300)
    let secret = intToFr(400)
    
    let values = @[oneFr, publicOut, publicIn, privateIn, secret]
    let witness = makeWitnessBN254(values)
    
    check witness.nvars == 5
    check isEqualFr(witness.values[0], oneFr)
    check isEqualFr(witness.values[1], publicOut)
    check isEqualFr(witness.values[2], publicIn)
    check isEqualFr(witness.values[3], privateIn)
    check isEqualFr(witness.values[4], secret)
