
#
# Ideally we want to the image of update (delta) of the witness under 
# both the A and B matrices to lie in a small subgroup
# 
# This can be always achieved simply by permuting the rows. 
# Here we compute such a permutation
# 

import std/sequtils
import std/algorithm
import std/tables

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/math/domain
import groth16/math/matrix
import groth16/zkey_types
import groth16/misc

#-------------------------------------------------------------------------------

type 

  F = Fr[BN254_Snarks]

  Permutation* = object
    srcIndices* : seq[int]

#-------------------------------------------------------------------------------

func isPermutation*(input: seq[int]): bool =
  let n = input.len
  let tgt: seq[int] = toSeq(0..<n)
  var tmp = input
  tmp.sort()
  return (tmp == tgt)

func isValidPermutation*(perm: Permutation): bool =
  return isPermutation(perm.srcIndices)

func inversePermutation*(P: Permutation): Permutation = 
  let N = P.srcIndices.len
  let ps: seq[int] = P.srcIndices
  var qs: seq[int] = newSeq[int]( N )
  for i in 0..<N:
    qs[ps[i]] = i
  return Permutation(srcIndices: qs)

func permuteSeq*[T](P: Permutation, xs: seq[T]) = 
  let N  = xs.len
  let ps = P.srcIndices
  assert( N == ps.len )

  var ys: seq[T] = newSeq[T]( N )
  for i in 0..<N:
    ys[i] = xs[ps[i]]

#-------------------------------------------------------------------------------

func permuteSparseColumn*[T]( perm: Permutation , column: SparseColumn[T] ): SparseColumn[T] =
  let tgtIndices = inversePermutation(perm).srcIndices
  var output: Table[int,T]
  for (i,y) in column:
    output[tgtIndices[i]] = y
  return output

func permuteMatrixRows*[T]( perm: Permutation , mat: SparseMatrix[T] ): SparseMatrix[T] =
  assert( mat.dims.nrows == perm.srcIndices.len )
  let tgtIndices = inversePermutation(perm).srcIndices
  var outputs: seq[SparseColumn[T]] = newSeq[SparseColumn[T]]( mat.dims.ncols )
  for j in 0..<mat.dims.ncols:
    var column: Table[int,T]
    for (i,y) in mat.columns[j].pairs:
      column[tgtIndices[i]] = y
    outputs[j] = column
  return SparseMatrix[T]( dims: mat.dims, columns: outputs )

#-------------------------------------------------------------------------------

func liesInSubgroup*( sg: Subgroup, mask: seq[bool] ): bool = 
  let N = sg.bigDomain.domainSize
  let K = sg.smallDomain.domainSize
  let J = N div K

  assert( N == mask.len )

  var ok = true
  for (i,b) in mask.pairs:
    if b:
      ok = ok and ((i mod J) == 0)

  return ok

#-------------------------------------------------------------------------------

proc findRowPermutation*(A: SparseMatrix[F], B: SparseMatrix[F], witness_delta_mask: seq[bool]): (int,Permutation) =

  let N = A.dims.nrows
  let M = A.dims.ncols

  assert( isPowerOfTwo(N) )
  assert( A.dims == B.dims )
  assert( witness_delta_mask.len == M )

  let deltaImgA  = sparseMatrixImage( A , witness_delta_mask )
  let deltaImgB  = sparseMatrixImage( B , witness_delta_mask )
  let deltaImgAB = orBoolSeqs( deltaImgA , deltaImgB )

  let k0   = countTrues(deltaImgAB)
  let logK = ceilingLog2(k0)
  let K    = (1 shl logK)           # size of the subgroup

  var sgIndex = N div K             # index of the subgroup

  # echo "domain size    = " & $N
  # echo "image size     = " & $k0
  # echo "subgroup size  = " & $K
  # echo "subgroup index = " & $sgIndex

  var used: seq[bool] = constSeq(N, false)
  var perm: seq[int]  = newSeq[int]( N )

  # first we fill in the subgroup part
  var k: int = 0                 
  for (i,b) in deltaImgAB.pairs:
    if b:
      let j = k * sgIndex
      perm[j] = i
      used[j] = true
      k += 1

  # next we fill in the rest
  var j: int = 0
  for (i,b) in deltaImgAB.pairs:
    if not b:
      while(used[j]): 
        j += 1
      perm[j] = i
      j += 1

  # # debugging
  # echo $trueIndices(deltaImgAB)
  # for j in 0..<K:
  #   echo $perm[sgIndex*j]

  let P = Permutation(srcIndices: perm)

  # sanity checks

  assert( isValidPermutation(P) , "findRowPermutation: the result is not valid permutation :(" )

  let A1 = permuteMatrixRows( P , A )
  let B1 = permuteMatrixRows( P , B )

  let newImgA  = sparseMatrixImage( A1 , witness_delta_mask )
  let newImgB  = sparseMatrixImage( B1 , witness_delta_mask )
  let newImgAB = orBoolSeqs( newImgA , newImgB )
 
  let sg = createSubgroup( createDomain(N) , K )

  assert( liesInSubgroup(sg , newImgAB) , "permuted matrices are not compatible with the subgroup...")

  return (K,P)
  
#-------------------------------------------------------------------------------

proc zkeyFindRowPermutation*(zkey: ZKey, witness_delta_mask: seq[bool]): (int,Permutation) =
  let matrices = zkeyToSparseMatrices(zkey)
  return findRowPermutation( matrices.A , matrices.B , witness_delta_mask )

#-------------------------------------------------------------------------------

