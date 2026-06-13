

{.used.} 

import std/unittest

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128/fields
import groth16/bn128/curves
import groth16/bn128/rnd
import groth16/bn128/debug

import groth16/math/domain
import groth16/math/ntt

#-------------------------------------------------------------------------------

proc randFrSeq(N: int) : seq[Fr[BN254_Snarks]] = 
  var xs : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( N )
  for i in 0..<N:
    xs[i] = randFr()
  return xs

proc scaleFrSeqInPlace(s: Fr[BN254_Snarks], arr: var seq[Fr[BN254_Snarks]] ) = 
  let N = arr.len  
  for i in 0..<N:
    arr[i] *= s

func isEqualFrSeq(xs : seq[Fr[BN254_Snarks]], ys: seq[Fr[BN254_Snarks]]): bool =
  let N = xs.len
  let M = ys.len

  if N != M:
    return false
  else:
    var ok: bool = true
    for i in 0..<N:
      ok = ok and (xs === ys)
    return ok

#---------------------------------------

proc randG1Seq(N: int) : seq[G1] = 
  var xs : seq[G1] = newSeq[G1]( N )
  for i in 0..<N:
    xs[i] = randG1()
  return xs

#-------------------------------------------------------------------------------

suite "FFT checks":

  suite "field NTT":

    let D = createDomain(128)
    let N = D.domainSize
  
    let xs : seq[Fr[BN254_Snarks]] = randFrSeq(N)
  
    test "INTT(NTT(xs) == xs":
      let ys = forwardNTT(xs, D)
      let zs = inverseNTT(ys, D) 
      check isEqualFrSeq( xs , zs )
    
    test "NTT(INTT(xs) == xs":
      let ys = inverseNTT(xs, D)
      let zs = forwardNTT(ys, D) 
      check isEqualFrSeq( xs , zs )

    test "unscaledINTT(NTT(xs) == N * xs":
      let ys = forwardNTT(xs, D)
      let zs = unscaledInverseNTT(ys, D) 
      var ws = xs
      scaleFrSeqInPlace(intToFr(N) , ws)
      check isEqualFrSeq( ws , zs )

#-------------------------------------------------------------------------------
  

