
import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128/fields
import groth16/bn128/curves
import groth16/bn128/rnd

#-------------------------------------------------------------------------------

# proper mathematical modulo operation.
# (all mainstream programming languages are fucked up... Except Haskell, Haskell gets it right)
#
# TODO: maybe move this somewhere else?
func safeMod*(x: int, N: int): int =
  if x >= 0:
    return (x mod N)
  else:
    let d = ((-x) div N)
    let y = x + (d+1)*N
    return (y mod N)

#-------------------------------------------------------------------------------
# generic arrays

func constSeq*[T]( N: int , y: T ): seq[T] =
  var xs : seq[T] = newSeq[T]( N )
  for i in 0..<N:
    xs[i] = y
  return xs

#-------------------------------------------------------------------------------
# Fr arrays

proc randFrSeq*(N: int) : seq[Fr[BN254_Snarks]] = 
  var xs : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( N )
  for i in 0..<N:
    xs[i] = randFr()
  return xs

func isEqualFrSeq*(xs: seq[Fr[BN254_Snarks]], ys: seq[Fr[BN254_Snarks]]): bool =
  let N = xs.len
  let M = ys.len

  if N != M:
    return false
  else:
    var ok: bool = true
    for i in 0..<N:
      ok = ok and (xs[i] === ys[i])
    return ok

proc scaleFrSeqInPlace*(s: Fr[BN254_Snarks], arr: var seq[Fr[BN254_Snarks]] ) = 
  let N = arr.len  
  for i in 0..<N:
    arr[i] *= s

proc pointwiseProdFr*(xs: seq[Fr[BN254_Snarks]], ys: seq[Fr[BN254_Snarks]]): seq[Fr[BN254_Snarks]] =
  let N = xs.len
  assert( N == ys.len )
  var zs : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( N )
  for i in 0..<N:
    zs[i] = xs[i] * ys[i]
  return zs

func dotProdFr*(xs, ys: seq[Fr[BN254_Snarks]]): Fr[BN254_Snarks] =
  let n = xs.len
  assert( n == ys.len, "dotProdFr: incompatible vector lengths" )
  var s : Fr[BN254_Snarks] = zeroFr
  for i in 0..<n:
    s += xs[i] * ys[i]
  return s

func sumSeqFr*(xs : seq[Fr[BN254_Snarks]]): Fr[BN254_Snarks] =
  let n = xs.len
  var s : Fr[BN254_Snarks] = zeroFr
  for i in 0..<n:
    s += xs[i]
  return s

#-------------------------------------------------------------------------------
# G1 arrays

proc randG1Seq*(N: int) : seq[G1] = 
  var xs : seq[G1] = newSeq[G1]( N )
  for i in 0..<N:
    xs[i] = randG1()
  return xs

func isEqualG1Seq*(xs: seq[G1], ys: seq[G1] ): bool =
  let N = xs.len
  let M = ys.len

  if N != M:
    return false
  else:
    var ok: bool = true
    for i in 0..<N:
      ok = ok and (xs[i] === ys[i])
    return ok

proc scaleG1SeqInPlace*(s: Fr[BN254_Snarks], arr: var seq[G1] ) = 
  let N = arr.len  
  for i in 0..<N:
    arr[i] = s ** arr[i]

proc pointwiseScaleG1*(xs: seq[Fr[BN254_Snarks]], gs: seq[G1]): seq[G1] =
  let N = xs.len
  assert( N == gs.len )
  var hs : seq[G1] = newSeq[G1]( N )
  for i in 0..<N:
    hs[i] = xs[i] ** gs[i]
  return hs

#-------------------------------------------------------------------------------
