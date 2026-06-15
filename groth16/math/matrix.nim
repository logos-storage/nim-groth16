
import std/tables

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

#-------------------------------------------------------------------------------
# dimensions

type
  MatrixDims* = object 
    nrows* : int
    ncols* : int

#-------------------------------------------------------------------------------
# Dense matrices
#
# Note: dense matrices can be very big, this is only feasible for small circuits

type 

  DenseColumn*[T] = seq[T]

  DenseMatrixColumns*[T] = seq[DenseColumn[T]]

  DenseMatrix*[T] = object
    dims*    : MatrixDims
    columns* : seq[DenseColumn[T]]

  DenseMatrices* = object
    A* : DenseMatrix[Fr[BN254_Snarks]]
    B* : DenseMatrix[Fr[BN254_Snarks]]
    C* : DenseMatrix[Fr[BN254_Snarks]]

#-------------------------------------------------------------------------------
# Sparse matrices

type 

  SparseColumn*[T] = Table[int,T]

  SparseMatrixColumns*[T] = seq[SparseColumn[T]]

  SparseMatrix*[T] = object
    dims*    : MatrixDims
    columns* : seq[SparseColumn[T]]

  SparseMatrices* = object
    A* : SparseMatrix[Fr[BN254_Snarks]]
    B* : SparseMatrix[Fr[BN254_Snarks]]
    C* : SparseMatrix[Fr[BN254_Snarks]]

proc columnInsertWithAddFr*( column: var SparseColumn[Fr[BN254_Snarks]] , row: int,  y: Fr[BN254_Snarks] ) =
  var x = getOrDefault( column, row, zeroFr )
  x += y
  column[row] = x

proc sparseDenseDotProdFr*( U: SparseColumn[Fr[BN254_Snarks]], V: DenseColumn[Fr[BN254_Snarks]] ): Fr[BN254_Snarks] =
  var acc : Fr[BN254_Snarks] = zeroFr
  for i,x in U.pairs:
    acc += x * V[i]
  return acc

#-------------------------------------------------------------------------------
# densities

# counts the non-zero elements in each row
func sparseMatrixRowCounts*( A : SparseMatrix[Fr[BN254_Snarks]] ): seq[int] =
  var rowCounts: seq[int] = newSeq[int]( A.dims.nrows )
  for j,column in A.columns.pairs:
    for i,value in column.pairs:
      if not isZeroFr(value):
        rowCounts[i] += 1
  return rowCounts

# counts the non-zero elements in each column
func sparseMatrixColumnCounts*( A : SparseMatrix[Fr[BN254_Snarks]] ): seq[int] =
  var colCounts: seq[int] = newSeq[int]( A.dims.ncols )
  for j,column in A.columns.pairs:
    for i,value in column.pairs:
      if not isZeroFr(value):
        colCounts[j] += 1
  return colCounts

# average count of non-zero elements in the rows
func sparseMatrixAvgRowDensity*( A : SparseMatrix[Fr[BN254_Snarks]] ): float64 =
  let rowCounts = sparseMatrixRowCounts( A )
  var s: float64 = 0
  for x in rowCounts: s += x.float64
  return (s / rowCounts.len.float64)

# average count of non-zero elements in the columns
func sparseMatrixAvgColumnDensity*( A : SparseMatrix[Fr[BN254_Snarks]] ): float64 =
  let colCounts = sparseMatrixColumnCounts( A )
  var s: float64 = 0
  for x in colCounts: s += x.float64
  return (s / colCounts.len.float64)

#-------------------------------------------------------------------------------

# image of subspace (of the witness space) 
func sparseMatrixImage*( A : SparseMatrix[Fr[BN254_Snarks]] , subspace: seq[bool] ): seq[bool] =
  assert( A.dims.ncols == subspace.len )
  var image: seq[bool] = newSeq[bool]( A.dims.nrows )
  for j,column in A.columns.pairs:
    if subspace[j]:
      for i,value in column.pairs:
        if not isZeroFr(value):
          image[i] = true
  return image

#-------------------------------------------------------------------------------
