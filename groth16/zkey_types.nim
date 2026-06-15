
import constantine/math/arithmetic
import constantine/named/properties_fields
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
    npubs*         : int             # number of public input/outputs (excluding the constant 1)
    domainSize*    : int             # size of the domain (should be power of two)
    logDomainSize* : int

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
# extract the three matrices from the ZKey

func coeffsToSparseMatrices*( dims: MatrixDims, coeffs: seq[Coeff]): SparseMatrices = 
  var A: SparseMatrixColumns[Fr[BN254_Snarks]] = newSeq[SparseColumn[Fr[BN254_Snarks]]] ( dims.ncols )
  var B: SparseMatrixColumns[Fr[BN254_Snarks]] = newSeq[SparseColumn[Fr[BN254_Snarks]]] ( dims.ncols )
  var C: SparseMatrixColumns[Fr[BN254_Snarks]] = newSeq[SparseColumn[Fr[BN254_Snarks]]] ( dims.ncols )
 
  for cf in coeffs:
    case cf.matrix
      of MatrixA: columnInsertWithAddFr( A[cf.col] , cf.row , cf.coeff )
      of MatrixB: columnInsertWithAddFr( B[cf.col] , cf.row , cf.coeff ) 
      of MatrixC: columnInsertWithAddFr( C[cf.col] , cf.row , cf.coeff )
 
  let mA = SparseMatrix[Fr[BN254_Snarks]]( dims: dims , columns: A )
  let mB = SparseMatrix[Fr[BN254_Snarks]]( dims: dims , columns: B )
  let mC = SparseMatrix[Fr[BN254_Snarks]]( dims: dims , columns: C )
  return SparseMatrices( A: mA, B: mB, C: mC )

func zkeyToSparseMatrices*(zkey: ZKey): SparseMatrices = 
  let dims = MatrixDims( nrows: zkey.header.domainSize , ncols: zkey.header.nvars )
  return coeffsToSparseMatrices( dims , zkey.coeffs )

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

proc printR1csStats*(zkey: ZKey) = 
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

