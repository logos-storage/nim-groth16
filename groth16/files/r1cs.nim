
#
# parsing the `.r1cs` file computed by `circom` witness code genereators
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
#
# 2: Constraints
# --------------
#   nConstr : word32    = number of constraints
#   then an array of constraints:
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

import constantine/math/arithmetic
import constantine/math/io/io_bigints
import constantine/named/properties_fields

import groth16/bn128
import groth16/files/container

#-------------------------------------------------------------------------------

type 
 
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
