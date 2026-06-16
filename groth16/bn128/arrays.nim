
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
# bool arrays

func countTrues*( bs: seq[bool] ): int =
  var cnt = 0 
  for b in bs: 
    if b: 
      cnt += 1
  return cnt

func countFalses*( bs: seq[bool] ): int =
  var cnt = 0 
  for b in bs: 
    if not b: 
      cnt += 1
  return cnt

func trueIndices*( bs: seq[bool] ): seq[int] =
  let k = countTrues(bs)
  var idxs: seq[int] = newSeq[int]( k )
  var j = 0
  for (i,b) in bs.pairs:
    if b:
      idxs[j] = i
      j += 1
  return idxs

func falseIndices*( bs: seq[bool] ): seq[int] =
  let k = countFalses(bs)
  var idxs: seq[int] = newSeq[int]( k )
  var j = 0
  for (i,b) in bs.pairs:
    if not b:
      idxs[j] = i
      j += 1
  return idxs

func selectTrues*[T]( bs: seq[bool] , xs: seq[T] ): seq[T] =
  assert( bs.len == xs.len )
  let k = countTrues(bs)
  var ys: seq[T] = newSeq[T]( k )
  var j = 0
  for (i,b) in bs.pairs:
    if b:
      ys[j] = xs[i]
      j += 1
  return ys

func notBoolSeq*( us: seq[bool]): seq[bool] =
  let n = us.len
  var ws: seq[bool] = newSeq[bool]( n )
  for i in 0..<n:
    ws[i] = not us[i]
  return ws

func andBoolSeqs*( us: seq[bool] , vs: seq[bool]): seq[bool] =
  let n = us.len
  assert( n == vs.len )
  var ws: seq[bool] = newSeq[bool]( n )
  for i in 0..<n:
    ws[i] = us[i] and vs[i]
  return ws

func orBoolSeqs*( us: seq[bool] , vs: seq[bool]): seq[bool] =
  let n = us.len
  assert( n == vs.len )
  var ws: seq[bool] = newSeq[bool]( n )
  for i in 0..<n:
    ws[i] = us[i] or vs[i]
  return ws

#-------------------------------------------------------------------------------
# int arrays

func sumIntSeq*(xs : seq[int]): int =
  let n = xs.len
  var s = 0
  for i in 0..<n:
    s += xs[i]
  return s

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

func countNonZerosFr*( xs: seq[Fr[BN254_Snarks]] ): int =
  var cnt = 0 
  for x in xs: 
    if not isZeroFr(x): 
      cnt += 1
  return cnt

# returns a mask and a filtered vector
func selectNonZerosFr*( xs: seq[Fr[BN254_Snarks]] ): (seq[bool] , seq[Fr[BN254_Snarks]]) =
  var cnt = 0 
  var mask: seq[bool] = newSeq[bool]( xs.len ) 
  for (i,x) in xs.pairs: 
    if isZeroFr(x):
      mask[i] = false
    else:
      mask[i] = true
      cnt += 1

  var k = 0
  var short: seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( cnt ) 
  for (i,x) in xs.pairs: 
    if mask[i]:
      short[k] = x
      k += 1

  return (mask,short)

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
