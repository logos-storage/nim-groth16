
import std/tables

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays

#-------------------------------------------------------------------------------
# dimensions

type
  
  F = Fr[BN254_Snarks]

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
    A* : DenseMatrix[F]
    B* : DenseMatrix[F]
    C* : DenseMatrix[F]

#-------------------------------------------------------------------------------
# Sparse matrices

type 

  SparseColumn*[T] = Table[int,T]

  SparseMatrixColumns*[T] = seq[SparseColumn[T]]

  SparseMatrix*[T] = object
    dims*    : MatrixDims
    columns* : seq[SparseColumn[T]]

  SparseMatrices* = object
    A* : SparseMatrix[F]
    B* : SparseMatrix[F]
    C* : SparseMatrix[F]

proc columnInsertWithAddFr*( column: var SparseColumn[F] , row: int,  y: F ) =
  var x = getOrDefault( column, row, zeroFr )
  x += y
  column[row] = x

proc sparseDenseDotProdFr*( U: SparseColumn[F], V: DenseColumn[F] ): F =
  var acc : F = zeroFr
  for i,x in U.pairs:
    acc += x * V[i]
  return acc

#-------------------------------------------------------------------------------
# sparse matrices represented in the transposed way (we store rows, not columns)

type 

  SparseRow*[T] = Table[int,T]

  SparseMatrixRows*[T] = seq[SparseRow[T]]

  SparseRowMatrix*[T] = object
    dims* : MatrixDims
    rows* : seq[SparseRow[T]]

proc rowInsertWithAddFr*( row: var SparseRow[F] , column: int,  y: F ) =
  var x = getOrDefault( row, column, zeroFr )
  x += y
  row[column] = x

func toSparseRowMatrix*( input: SparseMatrix[F] ): SparseRowMatrix[F] =
  let N = input.dims.nrows
  var rows: seq[SparseRow[F]] = newSeq[SparseRow[F]]( N )
  for (j,col) in input.columns.pairs:
    for (i,val) in col.pairs:
      rowInsertWithAddFr( rows[i] , j , val)
  return SparseRowMatrix[F]( dims: input.dims , rows: rows )

func fromSparseRowMatrix*( input: SparseRowMatrix[F] ): SparseMatrix[F] =
  let M = input.dims.ncols
  var columns: seq[SparseColumn[F]] = newSeq[SparseRow[F]]( M )
  for (i,row) in input.rows.pairs:
    for (j,val) in row.pairs:
      columnInsertWithAddFr( columns[j] , i , val )
  return SparseMatrix[F]( dims: input.dims , columns: columns )

#-----------------------------------------

# render a sparse row as a linear combination
func renderSparseRowLinComb*( row: SparseRow[F] ): string =
  var s:   string = "( "
  var cnt: int    = 0 
  for (j,val) in row.pairs:
    if not isZeroFr(val):
      let (sgn, abs) = renderSignedFr(val) 
      s &= (sgn & " " & abs & " * z" & $j & " ")
      cnt += 1
  if cnt == 0:
    s &= "0 "
  s &= ")"
  return s

# render three sparse rows as an R1CS equation
func renderSparseRowR1CSEq*( rowA, rowB, rowC: SparseRow[F] ): string =
  return
    (renderSparseRowLinComb(rowA) & " * " & 
     renderSparseRowLinComb(rowB) & " + " &        # not sure if the standard conventions is plus or minus here...
     renderSparseRowLinComb(rowC) & " == 0" )

#-------------------------------------------------------------------------------
# densities

# counts the non-zero elements in each row
func sparseMatrixRowCounts*( A : SparseMatrix[F] ): seq[int] =
  var rowCounts: seq[int] = newSeq[int]( A.dims.nrows )
  for j,column in A.columns.pairs:
    for i,value in column.pairs:
      if not isZeroFr(value):
        rowCounts[i] += 1
  return rowCounts

# counts the non-zero elements in each column
func sparseMatrixColumnCounts*( A : SparseMatrix[F] ): seq[int] =
  var colCounts: seq[int] = newSeq[int]( A.dims.ncols )
  for j,column in A.columns.pairs:
    for i,value in column.pairs:
      if not isZeroFr(value):
        colCounts[j] += 1
  return colCounts

# average count of non-zero elements in the rows
func sparseMatrixAvgRowDensity*( A : SparseMatrix[F] ): float64 =
  let rowCounts = sparseMatrixRowCounts( A )
  var s: float64 = 0
  for x in rowCounts: s += x.float64
  return (s / rowCounts.len.float64)

# average count of non-zero elements in the columns
func sparseMatrixAvgColumnDensity*( A : SparseMatrix[F] ): float64 =
  let colCounts = sparseMatrixColumnCounts( A )
  var s: float64 = 0
  for x in colCounts: s += x.float64
  return (s / colCounts.len.float64)

#-------------------------------------------------------------------------------

# image of subspace (of the witness space) 
func sparseMatrixImage*( A : SparseMatrix[F] , subspace: seq[bool] ): seq[bool] =
  assert( A.dims.ncols == subspace.len )
  var image: seq[bool] = newSeq[bool]( A.dims.nrows )
  for j,column in A.columns.pairs:
    if subspace[j]:
      for i,value in column.pairs:
        if not isZeroFr(value):
          image[i] = true
  return image

#-------------------------------------------------------------------------------
