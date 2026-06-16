
#
# parsing and exporting the `.r1cs` file computed by `circom` witness code genereators
#
# file format
# ===========
# 
# standard iden3 binary container format.
# field elements are in standard representation
#
# sections:
#
# 1: Header
# ---------
#   n8r     : word32    = how many bytes are a field element in Fr
#   r       : n8r bytes = the size of the prime field Fr (the scalar field)
#   nWires  : word32    = number of wires (or witness variables)
#   nPubOut : word32    = number of public outputs
#   nPubIn  : word32    = number of public inputs
#   nPrivIn : word32    = number of private inputs
#   nLabels : word64    = number of labels (variable names in the circom source code)
#   nConstr : word32    = number of constraints
#
# 2: Constraints
# --------------
#   An array of constraints:
#     A : LinComb
#     B : LinComb
#     C : LinComb
#   meaning `A*B=C`, where LinComb looks like this:
#     nTerms : word32     = number of terms
#     <an array of terms>
#   where a term looks like this:
#     idx   : word32      = which witness variable
#     coeff : Fr          = the coefficient
#     
# 3: Wire-to-label mapping
# ------------------------
#   <an array of `nWires` many 64 bit words>
#
# 4: Custom gates list
# --------------------
#   ...
#   ...
#
# 5: Custom gates application
# ---------------------------
#   ...
#   ...
#

import std/streams
import std/sequtils

import constantine/math/arithmetic
import constantine/math/io/io_bigints
import constantine/named/properties_fields

import groth16/bn128
import groth16/math/matrix
import groth16/files/container

#-------------------------------------------------------------------------------

type 
 
  F = Fr[BN254_Snarks]

  WitnessConfig* = object
    nWires*  : int           # total number of wires (or witness variables), including the constant 1 "variable"
    nPubOut* : int           # number of public outputs
    nPubIn*  : int           # number of public inputs
    nPrivIn* : int           # number of private inputs
    nLabels* : int           # number of labels

  Term*       = tuple[ wireIdx: int, value: Fr[BN254_Snarks] ]
  LinComb*    = seq[Term]
  Constraint* = tuple[ A: LinComb, B: LinComb, C: LinComb ]

  R1CS* = object
    r*           : BigInt[256] 
    cfg*         : WitnessConfig
    nConstr*     : int
    constraints* : seq[Constraint]
    wireToLabel* : seq[int]

const emptyLinComb*    : LinComb    = @[]
const emptyConstraint* : Constraint = ( emptyLinComb, emptyLinComb, emptyLinComb )

#-------------------------------------------------------------------------------

proc printWitnessConfig*(cfg: WitnessConfig) = 
  echo "R1CS witness config:"
  echo " - nWires  = " & $cfg.nWires  
  echo " - nPubOut = " & $cfg.nPubOut 
  echo " - nPubIn  = " & $cfg.nPubIn  
  echo " - nPrivIn = " & $cfg.nPrivIn 
  echo " - nLabels = " & $cfg.nLabels 

proc printR1CSMetaData*(r1cs: R1CS) = 
  printWitnessConfig(r1cs.cfg)
  echo "nConstraints = " & $r1cs.nConstr
  assert (r1cs.nConstr == r1cs.constraints.len)

#
# This is to stop a (nontrivial) malleability attack; see:
#
# <https://web.archive.org/web/20240320152158/https://geometry.xyz/notebook/groth16-malleability>
# 
# Note: Snarkjs does it automatically (as does Arkworks), however, it's hidden
# in the ZKey generation proces... And it completely ruins our subgroup image thing
#
# So better do it manually, and patch Snarkjs not to do it if already presents...
#
proc r1csAddMalleabilityEquations*(orig: R1CS): R1CS = 

  # the one is the constant 1 (0th witness element)
  let npubs = 1 + orig.cfg.nPubOut + orig.cfg.nPubIn

  var dummyEqs : seq[Constraint] = newSeq[Constraint]( npubs )
  for i in 0..<npubs:
    let lcA = @[ ( i , oneFr ) ]
    dummyEqs[i] = ( A: lcA , B: emptyLinComb , C: emptyLinComb )

  var r1cs = orig
  r1cs.nConstr     = orig.nConstr + npubs
  r1cs.constraints = concat( orig.constraints , dummyEqs )
  return r1cs

#-------------------------------------------------------------------------------

func r1csToSparseMatrices*( r1cs: R1CS ): SparseMatrices =
  let N = r1cs.nConstr
  let M = r1cs.cfg.nWires
  let dims : MatrixDims = MatrixDims( nrows: N, ncols: M )

  var colsA: seq[SparseColumn[F]] = newSeq[SparseColumn[F]]( M )  
  var colsB: seq[SparseColumn[F]] = newSeq[SparseColumn[F]]( M )  
  var colsC: seq[SparseColumn[F]] = newSeq[SparseColumn[F]]( M )  

  for (i, constraint) in r1cs.constraints.pairs:
    let (lcA, lcB, lcC) = constraint
    for (j,x) in lcA: columnInsertWithAddFr( colsA[j] , i , x )
    for (j,y) in lcB: columnInsertWithAddFr( colsB[j] , i , y )
    for (j,z) in lcC: columnInsertWithAddFr( colsC[j] , i , z )

  let mA : SparseMatrix[F] = SparseMatrix[F]( dims: dims, columns: colsA )  
  let mB : SparseMatrix[F] = SparseMatrix[F]( dims: dims, columns: colsB )  
  let mC : SparseMatrix[F] = SparseMatrix[F]( dims: dims, columns: colsC ) 

  return SparseMatrices( A: mA , B: mB , C: mC )

#-------------------------------------------------------------------------------

proc parseSection1_header( stream: Stream, user: var R1CS, sectionLen: int ) =
  # echo "\nparsing r1cs header"
  
  let (n8r, r) = parsePrimeField( stream )     # size of the scalar field
  user.r = r;

  # echo("r = ",toDecimalBig(r))

  assert( sectionLen == 4 + n8r + 16 + 8 + 4, "unexpected section length")

  assert( bool(r == primeR) , "expecting the alt-bn128 curve" )

  var cfg : WitnessConfig

  cfg.nWires  = int( stream.readUint32() )
  cfg.nPubOut = int( stream.readUint32() )
  cfg.nPubIn  = int( stream.readUint32() )
  cfg.nPrivIn = int( stream.readUint32() )
  cfg.nLabels = int( stream.readUint64() )
  user.cfg = cfg

  let nConstr = int( stream.readUint32() )
  user.nConstr = nConstr

  # echo("witness config = ",cfg)
  # echo("nConstr = ",nConstr)

#-------------------------------------------------------------------------------

proc loadTerm( stream: Stream ): Term = 
  let idx   = int( stream.readUint32() )
  let coeff = loadValueFrStd( stream )
  return (wireIdx:idx, value:coeff)

proc loadLinComb( stream: Stream ): LinComb = 
  let nterms = int( stream.readUint32() )
  var terms : seq[Term]
  for i in 1..nterms:
    terms.add( loadTerm(stream) )
  return terms

proc loadConstraint( stream: Stream ): Constraint = 
  let a = loadLinComb( stream )
  let b = loadLinComb( stream )
  let c = loadLinComb( stream )
  return (A:a, B:b, C:c)

#-------------------------------------------------------------------------------

proc parseSection2_constraints( stream: Stream, user: var R1CS, sectionLen: int ) =
  var constr: seq[Constraint]
  var ncoeffsA, ncoeffsB, ncoeffsC: int
  for i in 1..(user.nConstr):
    let abc = loadConstraint(stream)
    constr.add( abc )
    ncoeffsA += abc.A.len
    ncoeffsB += abc.B.len
    ncoeffsC += abc.C.len
  user.constraints = constr
  # echo( "number of nonzero coefficients in matrix A = ", ncoeffsA )
  # echo( "number of nonzero coefficients in matrix B = ", ncoeffsB )
  # echo( "number of nonzero coefficients in matrix C = ", ncoeffsC )

#-------------------------------------------------------------------------------

proc parseSection3_wireToLabel( stream: Stream, user: var R1CS, sectionLen: int ) =
  assert( sectionLen == 8 * user.cfg.nWires, "unexpected section length")
  var labels: seq[int]
  for i in 1..(user.cfg.nWires):
    let label = int( stream.readUint64() )
    labels.add( label )
  user.wireToLabel = labels

#-------------------------------------------------------------------------------

proc r1csCallback( stream:  Stream
                 , sectId:  int
                 , sectLen: int
                 , user:    var R1CS
                 ) = 
  case sectId
    of 1: parseSection1_header(      stream, user, sectLen )
    of 2: parseSection2_constraints( stream, user, sectLen )
    of 3: parseSection3_wireToLabel( stream, user, sectLen )
    else: discard

proc parseR1CS* (fname: string): R1CS = 
  var r1cs : R1CS
  parseContainer( "r1cs", 1, fname, r1cs, r1csCallback, proc (id: int): bool = id == 1 )
  parseContainer( "r1cs", 1, fname, r1cs, r1csCallback, proc (id: int): bool = id != 1 )
  return r1cs

#-------------------------------------------------------------------------------
# writing R1CS files (required for reordering the rows)

proc writeWitnessConfig(stream: Stream, cfg: WitnessConfig ) =
  write[uint32]( stream , cfg.nWires.uint32  )       # total number of wires (or witness variables)
  write[uint32]( stream , cfg.nPubOut.uint32 )       # number of public outputs
  write[uint32]( stream , cfg.nPubIn.uint32  )       # number of public inputs
  write[uint32]( stream , cfg.nPrivIn.uint32 )       # number of private inputs
  write[uint64]( stream , cfg.nLabels.uint32 )       # number of labels

proc writeR1CSHeader(stream: Stream, r1cs: R1CS ) =
  write[uint32](      stream , 32       )            # n8r
  write(              stream , r1cs.r   )            # the value of r
  writeWitnessConfig( stream , r1cs.cfg )            # witness config
  write[uint32]( stream , r1cs.nConstr.uint32 )      # number of constraints

proc r1csHeaderSection*(r1cs: R1CS): seq[byte] = 
  var stream = newStringStream()
  stream.writeR1CSHeader(r1cs)
  stream.flush()
  stream.setPosition(0)
  let bytes = cast[seq[byte]](stream.readAll())      # WTF nim  
  stream.close()
  return bytes

proc r1csWireToLabelSection*(r1cs: R1CS): seq[byte] = 
  var stream = newStringStream()
  for x in r1cs.wireToLabel:
    write[uint64]( stream , x.uint64 )
  stream.flush()
  stream.setPosition(0)
  let bytes = cast[seq[byte]](stream.readAll())      # WTF nim  
  stream.close()
  return bytes

proc writeFr(stream: Stream, x: Fr[BN254_Snarks]) =
  var big : BigInt[254]  # fucking constantine
  big = x.toBig()
  stream.write(big)

proc writeTerm(stream: Stream, term: Term ) = 
  write[uint32]( stream , term.wireIdx.uint32 )
  writeFr( stream , term.value )

proc writeLinComb(stream: Stream, lc: LinComb ) =
  write[uint32]( stream , lc.len.uint32 )
  for term in lc:
    writeTerm(stream, term)

proc writeConstraint(stream: Stream, con: Constraint ) =
  writeLinComb( stream, con.A )
  writeLinComb( stream, con.B )
  writeLinComb( stream, con.C )

proc r1csConstraintsSection*(r1cs: R1CS): seq[byte] = 
  var stream = newStringStream()
  for c in r1cs.constraints:
    writeConstraint( stream, c )
  stream.flush()
  stream.setPosition(0)
  let bytes = cast[seq[byte]](stream.readAll())      # WTF nim  
  stream.close()
  return bytes

proc exportR1CS*(fname: string, r1cs: R1CS) =
  let section1 = r1csHeaderSection(     r1cs)
  let section2 = r1csConstraintsSection(r1cs)
  let section3 = r1csWireToLabelSection(r1cs)
  var stream = newFileStream(fname, fmWrite)
  writeGlobalHeader(  stream , "r1cs" , 1 , 3 )
  writeSection(       stream , 2 , section2 )
  writeSection(       stream , 1 , section1 )
  writeSection(       stream , 3 , section3 )
  stream.flush()
  stream.close()

#-------------------------------------------------------------------------------
