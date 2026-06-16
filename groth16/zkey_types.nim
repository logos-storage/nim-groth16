

import std/tables

import constantine/named/properties_fields
import constantine/math/arithmetic
import constantine/math/extension_fields/towers

import groth16/bn128
import groth16/math/matrix

#-------------------------------------------------------------------------------

type 
 
  Flavour* = enum
    JensGroth          # the version described in the original Groth16 paper
    Snarkjs            # the version implemented by Snarkjs

  GrothHeader* = object
    curve*         : string          # name of the curve, eg. "bn128"
    flavour*       : Flavour         # which variation of the trusted setup
    p*             : BigInt[256]     # size of the base field
    r*             : BigInt[256]     # size of the scalar field
    nvars*         : int             # number of witness variables (including the constant 1)
    npubs*         : int             # number of public input/outputs (**excluding** the constant 1)
    domainSize*    : int             # size of the domain (should be power of two)
    logDomainSize* : int             # base 2 logarithm of the size of the domain

  SpecPoints* = object
    alpha1*      : G1                # = alpha * g1
    beta1*       : G1                # = beta  * g1
    beta2*       : G2                # = beta  * g2
    gamma2*      : G2                # = gamma * g2
    delta1*      : G1                # = delta * g1
    delta2*      : G2                # = delta * g2       
    alphaBeta*   : Fp12[BN254_Snarks] # = <alpha1 , beta2>

  VerifierPoints* = object
    pointsIC*    : seq[G1]           # the points `delta^-1 * ( beta*A_j(tau) + alpha*B_j(tau) + C_j(tau) ) * g1` (for j <= npub)

  ProverPoints* = object
    pointsA1*    : seq[G1]           # the points `A_j(tau) * g1`
    pointsB1*    : seq[G1]           # the points `B_j(tau) * g1`
    pointsB2*    : seq[G2]           # the points `B_j(tau) * g2`
    pointsC1*    : seq[G1]           # the points `delta^-1 * ( beta*A_j(tau) + alpha*B_j(tau) + C_j(tau) ) * g1` (for j > npub)
    pointsH1*    : seq[G1]           # meaning depends on `flavour`

  MatrixSel* = enum
    MatrixA
    MatrixB
    MatrixC

  Coeff* = object
    matrix* : MatrixSel
    row*    : int
    col*    : int
    coeff*  : Fr[BN254_Snarks]

  ZKey* = object
    # sectionMask* : uint32
    header*      : GrothHeader
    specPoints*  : SpecPoints
    vPoints*     : VerifierPoints
    pPoints*     : ProverPoints
    coeffs*      : seq[Coeff]

  VKey* = object 
    curve*   : string
    spec*    : SpecPoints
    vpoints* : VerifierPoints

#-------------------------------------------------------------------------------

# extract the verification key
func extractVKey*(zkey: Zkey): VKey = 
  let curve = zkey.header.curve
  let spec  = zkey.specPoints
  let vpts  = zkey.vPoints
  return VKey(curve:curve, spec:spec, vpoints:vpts)

#-------------------------------------------------------------------------------
# extract the three matrices from the ZKey (TODO: C is actually not present there...)

type F = Fr[BN254_Snarks]

func coeffsToSparseMatrices*( dims: MatrixDims, coeffs: seq[Coeff]): SparseMatrices = 
  var A: SparseMatrixColumns[F] = newSeq[SparseColumn[F]] ( dims.ncols )
  var B: SparseMatrixColumns[F] = newSeq[SparseColumn[F]] ( dims.ncols )
  var C: SparseMatrixColumns[F] = newSeq[SparseColumn[F]] ( dims.ncols )
 
  for cf in coeffs:
    case cf.matrix
      of MatrixA: columnInsertWithAddFr( A[cf.col] , cf.row , cf.coeff )
      of MatrixB: columnInsertWithAddFr( B[cf.col] , cf.row , cf.coeff ) 
      of MatrixC: columnInsertWithAddFr( C[cf.col] , cf.row , cf.coeff )
 
  let mA = SparseMatrix[F]( dims: dims , columns: A )
  let mB = SparseMatrix[F]( dims: dims , columns: B )
  let mC = SparseMatrix[F]( dims: dims , columns: C )
  return SparseMatrices( A: mA, B: mB, C: mC )

func zkeyToSparseMatrices*(zkey: ZKey): SparseMatrices = 
  let dims = MatrixDims( nrows: zkey.header.domainSize , ncols: zkey.header.nvars )
  return coeffsToSparseMatrices( dims , zkey.coeffs )

#-------------------------------------------------------------------------------
# useful for debugging

func doesPublicIOAppears( npubs: int , row: SparseRow[F] , do_include_const1: bool ): bool =
  var yes   = false
  let start = (if do_include_const1: 0 else: 1)
  for j in start..<npubs:
    if row.contains(j):
      yes = true
      break
  return yes

proc zkeyShowPublicIOEquations*(zkey: ZKey, do_include_const1: bool) =
  let N     = zkey.header.domainSize
  let npubs = zkey.header.npubs + 1                  # NOTE: zkey's npubs does NOT include the constant 1 !!!!!!
  let mats  = zkeyToSparseMatrices(zkey)
  let A = toSparseRowMatrix(mats.A)
  let B = toSparseRowMatrix(mats.B)
  let C = toSparseRowMatrix(mats.C)

  echo ""
  for i in 0..<N:
    let rowA = A.rows[i]
    let rowB = B.rows[i]
    let rowC = C.rows[i]

    if ( doesPublicIOAppears( npubs , rowA , do_include_const1 ) or
         doesPublicIOAppears( npubs , rowB , do_include_const1 ) or
         doesPublicIOAppears( npubs , rowC , do_include_const1 ) ):   
      echo "row #" & $i & ": " & renderSparseRowR1CSEq( rowA , rowB , rowC )

  echo ""

#-------------------------------------------------------------------------------

proc printGrothHeader*(hdr: GrothHeader) = 
  echo("")
  echo("curve         = " & ($hdr.curve        ) ) 
  echo("flavour       = " & ($hdr.flavour      ) ) 
  echo("|Fp|          = " & (toDecimalBig(hdr.p)) ) 
  echo("|Fr|          = " & (toDecimalBig(hdr.r)) ) 
  echo("nvars         = " & ($hdr.nvars        ) ) 
  echo("npubs         = " & ($hdr.npubs        ) ) 
  echo("domainSize    = " & ($hdr.domainSize   ) ) 
  echo("logDomainSize = " & ($hdr.logDomainSize) ) 

#-------------------------------------------------------------------------------

proc printMatrixStats*(zkey: ZKey) = 
  let matrices = zkeyToSparseMatrices(zkey)
  echo "average row density of A = " & $(sparseMatrixAvgRowDensity(matrices.A))
  echo "average row density of B = " & $(sparseMatrixAvgRowDensity(matrices.B))
  # echo "average row density of C = " & $(sparseMatrixAvgRowDensity(matrices.C))
  # the matrix C is not included in ZKey

#-------------------------------------------------------------------------------
# debugging

func matrixSelToString(sel: MatrixSel): string = 
  case sel 
    of MatrixA: return "A"
    of MatrixB: return "B"
    of MatrixC: return "C"

proc debugPrintCoeff(cf: Coeff) = 
  echo(    "matrix=", matrixSelToString(cf.matrix)
      , " | i=", cf.row
      , " | j=", cf.col
      , " | val=", signedToDecimalFr(cf.coeff)
      )

proc debugPrintCoeffs*(cfs: seq[Coeff]) = 
  for cf in cfs: debugPrintCoeff(cf)

#-------------------------------------------------------------------------------

