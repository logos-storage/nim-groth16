
{.used.} 

import std/unittest
import std/sequtils
import std/sugar

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128/fields
# import groth16/bn128/rnd
# import groth16/bn128/debug

#-------------------------------------------------------------------------------

suite "field tests":

  test "one-half in Fr":
    let two    = oneFr + oneFr
    let invTwo = oneHalfFr
    check (oneFr === invTwo * two)

  test "sanity check for batch inverse over Fr":
    let xs = map( toSeq(101..137) , intToFr )
    let ys = batchInverseFr( xs )
    let zs = collect( newSeq, (for x in xs: invFr(x)) )
    let n = xs.len
    # for i in 0..<n: echo(i," | batch = ",toDecimalFr(ys[i])," | ref = ",toDecimalFr(zs[i]) )

    var ok: bool = true
    for i in 0..<n:
      ok = ok and bool(ys[i] == zs[i])

    check ok

#-------------------------------------------------------------------------------

