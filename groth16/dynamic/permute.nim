
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
import groth16/files/r1cs
import groth16/zkey_types
import groth16/misc

#-------------------------------------------------------------------------------

type 

  F = Fr[BN254_Snarks]

  Permutation* = object
    srcIndices* : seq[int]

#-------------------------------------------------------------------------------

func permutationSize*(perm: Permutation): int = 
  return perm.srcIndices.len

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

func permuteR1CSRows*( perm: Permutation , r1cs: R1CS ): R1CS =
  let N = perm.srcIndices.len 
  assert( r1cs.nConstr <= N )
  var newConstraints: seq[Constraint] = newSeq[Constraint]( N )
  for i in 0..<N:
    let k = perm.srcIndices[i]
    if k < r1cs.nConstr:
      newConstraints[i] = r1cs.constraints[k]
    else:
      newConstraints[i] = emptyConstraint
  var newR1CS = r1cs
  newR1CS.nConstr     = N
  newR1CS.constraints = newConstraints
  return newR1CS

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

proc r1csPermuteToSubgroup*(r1cs_orig: R1CS, witness_delta_mask: seq[bool], do_add_malleability_eqs: bool = true ): R1CS =
 
  var r1cs: R1CS
  if do_add_malleability_eqs: 
    r1cs = r1csAddMalleabilityEquations( r1cs_orig )
  else:
    r1cs = r1cs_orig

  let log2  = ceilingLog2(r1cs.nConstr)
  let N     = (1 shl log2)
  let npubs = r1cs.cfg.nPubIn + r1cs.cfg.nPubOut + 1
  let N1    = N - npubs
  assert( r1cs.nConstr <= N1  , "snarkjs adds some extra rows for the public inputs" )
  let M    = r1cs.cfg.nWires
  let newDims = MatrixDims( nrows: N, ncols: M )
  let mats = r1csToSparseMatrices( r1cs )
  var A = mats.A ; A.dims = newDims
  var B = mats.B ; B.dims = newDims
  var (K, perm) = findRowPermutation( A , B , witness_delta_mask )

  # snarkjs adds some extra row for the public inputs...
  # so we need to subtract that many to have enough space
  perm.srcIndices = toSeq(perm.srcIndices[0..<N1])
  
  echo "witness update size = " & $countTrues(witness_delta_mask)
  echo "total constraints   = " & $permutationSize(perm)
  echo "subgroup size       = " & $K
  return permuteR1CSRows( perm , r1cs )

proc exportPermutedR1CS*(fname: string, r1cs: R1CS, witness_delta_mask: seq[bool] ) =
  let newR1CS = r1csPermuteToSubgroup( r1cs, witness_delta_mask )
  exportR1CS(fname, newR1CS)

#-------------------------------------------------------------------------------


