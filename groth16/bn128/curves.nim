#
# the `BN254` (aka `alt-bn128`) elliptic curve
#
# See for example <https://hackmd.io/@jpw/bn254>
#
#   p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
#   r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
#
# equation: y^2 = x^3 + 3
#

import std/bitops
import std/options

#import constantine/platforms/abstractions
#import constantine/math/isogenies/frobenius

import constantine/math/io/io_bigints
import constantine/math/arithmetic
import constantine/math/io/io_fields

import constantine/named/properties_fields      as tff
import constantine/math/extension_fields/towers as ext

import constantine/math/elliptic/ec_shortweierstrass_affine     as aff 
import constantine/math/elliptic/ec_shortweierstrass_projective as prj 
import constantine/math/pairings/pairings_bn                    as ate 
import constantine/math/elliptic/ec_scalar_mul_vartime          as scl 

import constantine/math/arithmetic/finite_fields_square_root    as sqrt
import constantine/math/extension_fields/square_root_fp2        as sqrt2

import groth16/bn128/fields

#-------------------------------------------------------------------------------

type G1*   = aff.EC_ShortW_Aff[Fp[BN254_Snarks] , aff.G1]
type G2*   = aff.EC_ShortW_Aff[Fp2[BN254_Snarks], aff.G2]

type ProjG1*  = prj.EC_ShortW_Prj[Fp[BN254_Snarks] , prj.G1]
type ProjG2*  = prj.EC_ShortW_Prj[Fp2[BN254_Snarks], prj.G2]

#-------------------------------------------------------------------------------
# compressed points (supposedly compatible with arkworks-0.5)

type ComprG1* = distinct array[32, byte];
type ComprG2* = distinct array[64, byte];

proc `==` *(a, b: ComprG1): bool {.borrow.}
proc `==` *(a, b: ComprG2): bool {.borrow.}

#-------------------------------------------------------------------------------

func isEqualG1* (x, y: G1 ): bool = bool(x == y)
func isEqualG2* (x, y: G2 ): bool = bool(x == y)

func `===`*(x, y: G1 ): bool = isEqualG1(x,y)
func `===`*(x, y: G2 ): bool = isEqualG2(x,y)

#-------------------------------------------------------------------------------

func unsafeMkG1* ( X, Y: Fp[BN254_Snarks] ) : G1 =
  return aff.EC_ShortW_Aff[Fp[BN254_Snarks], aff.G1](x: X, y: Y)

func unsafeMkG2* ( X, Y: Fp2[BN254_Snarks] ) : G2 =
  return aff.EC_ShortW_Aff[Fp2[BN254_Snarks], aff.G2](x: X, y: Y)

#-------------------------------------------------------------------------------

const infG1* : G1  = unsafeMkG1( zeroFp  , zeroFp  )
const infG2* : G2  = unsafeMkG2( zeroFp2 , zeroFp2 )

func isInfG1*(pt : G1): bool = bool(isNeutral(pt))
func isInfG2*(pt : G2): bool = bool(isNeutral(pt))

func isInfProjG1*(pt : ProjG1): bool = bool(isNeutral(pt))
func isInfProjG2*(pt : ProjG2): bool = bool(isNeutral(pt))

#-------------------------------------------------------------------------------

# y^2 = x^3 + B where B = 3
const theCoeffB = fromHex(Fp[BN254_Snarks], "0x0000000000000000000000000000000000000000000000000000000000000003")

func checkCurveEqG1*( x, y: Fp[BN254_Snarks] ) : bool =
  if bool(isZero(x)) and bool(isZero(y)):
    # the point at infinity is on the curve by definition
    return true
  else:
    var x2 = squareFp(x)
    var y2 = squareFp(y)
    var x3 = x2 * x
    var eq : Fp[BN254_Snarks]
    eq =  x3
    eq += theCoeffB
    eq -= y2
    # echo("eq = ",toDecimalFp(eq))
    return (bool(isZero(eq)))

# note: for BN254, the G1 is the whole curve. This is however not true for other curves like BLS12-381!
func checkSubgroupG1*( x, y: Fp[BN254_Snarks] ) : bool =
  return checkCurveEqG1(x,y)

#---------------------------------------

# y^2 = x^3 + B
# B = b1 + bu*u
# b1 = 19485874751759354771024239261021720505790618469301721065564631296452457478373
# b2 = 266929791119991161246907387137283842545076965332900288569378510910307636690
const twistCoeffB_1 = fromHex(Fp[BN254_Snarks], "0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5")
const twistCoeffB_u = fromHex(Fp[BN254_Snarks], "0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2")
const twistCoeffB   = mkFp2( twistCoeffB_1 , twistCoeffB_u )

func checkCurveEqG2*( x, y: Fp2[BN254_Snarks] ) : bool =
  if isZeroFp2(x) and isZeroFp2(y):
    # the point at infinity is on the curve by definition
    return true
  else:
    var x2 = squareFp2(x)
    var y2 = squareFp2(y)
    var x3 = x2 * x
    var eq : Fp2[BN254_Snarks]
    eq =  x3
    eq += twistCoeffB
    eq -= y2
    return isZeroFp2(eq)

# both just fits into 254 bits
const G2_cofactor: BigInt[254] = fromHex( BigInt[254] , "0x30644e72e131a029b85045b68181585e06ceecda572a2489345f2299c0f9fa8d" , bigEndian )
const GroupOrder : BigInt[254] = fromHex( BigInt[254] , "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001" , bigEndian )

# checks that the point is in the right subgroup
func checkSubgroupG2*( x, y: Fp2[BN254_Snarks] ) : bool =

  if not checkCurveEqG2(x,y):
    # reject if not on the curve
    return false
  else:
    # multiply by the order. Note: this is slow
    # TODO: implement <https://hackmd.io/@jpw/bn254#mathbb-G_2-membership-check-using-efficient-endomorphism>
    let point = unsafeMkG2(x,y)
    var q : ProjG2
    prj.fromAffine( q , point )
    scl.scalarMul_vartime(  q , GroupOrder )
    return bool(isNeutral(q))
    # var r : G2
    # prj.affine( r, q )

#-------------------------------------------------------------------------------

func mkG1*( x, y: Fp[BN254_Snarks] ) : G1 =
  if isZeroFp(x) and isZeroFp(y):
    return infG1
  else:
    assert( checkCurveEqG1(x,y) , "mkG1: not a G1 group point (in case of BN254, this is the same a curve point)" )
    return unsafeMkG1(x,y)

func mkCurve2*( x, y: Fp2[BN254_Snarks] ) : G2 =
  if isZeroFp2(x) and isZeroFp2(y):
    return infG2
  else:
    assert( checkCurveEqG2(x,y) , "mkCurve2: not a curve point on the curve over the extended field" )
    return unsafeMkG2(x,y)

func mkG2*( x, y: Fp2[BN254_Snarks] ) : G2 =
  if isZeroFp2(x) and isZeroFp2(y):
    return infG2
  else:
    assert( checkSubgroupG2(x,y) , "mkG2: not a G2 group point" )
    return unsafeMkG2(x,y)

#-------------------------------------------------------------------------------
# group generators

const gen1_x  = fromHex(Fp[BN254_Snarks], "0x01")
const gen1_y  = fromHex(Fp[BN254_Snarks], "0x02")

const gen2_xi  = fromHex(Fp[BN254_Snarks], "0x1adcd0ed10df9cb87040f46655e3808f98aa68a570acf5b0bde23fab1f149701")
const gen2_xu  = fromHex(Fp[BN254_Snarks], "0x09e847e9f05a6082c3cd2a1d0a3a82e6fbfbe620f7f31269fa15d21c1c13b23b")
const gen2_yi  = fromHex(Fp[BN254_Snarks], "0x056c01168a5319461f7ca7aa19d4fcfd1c7cdf52dbfc4cbee6f915250b7f6fc8")
const gen2_yu  = fromHex(Fp[BN254_Snarks], "0x0efe500a2d02dd77f5f401329f30895df553b878fc3c0dadaaa86456a623235c")

const gen2_x   = mkFp2( gen2_xi, gen2_xu )
const gen2_y   = mkFp2( gen2_yi, gen2_yu )

const gen1* : G1 = unsafeMkG1( gen1_x, gen1_y )
const gen2* : G2 = unsafeMkG2( gen2_x, gen2_y )

#-------------------------------------------------------------------------------

func isOnCurve1* ( p: G1 ) : bool =
  return checkCurveEqG1( p.x, p.y )

func isOnCurve2* ( p: G2 ) : bool =
  return checkCurveEqG2( p.x, p.y )

func isInSubgroupG1* ( p: G1 ) : bool =
  return checkSubgroupG1( p.x, p.y )

func isInSubgroupG2* ( p: G2 ) : bool =
  return checkSubgroupG2( p.x, p.y )

#-------------------------------------------------------------------------------

func unwrapComprG1*( c1: ComprG1 ): array[32,byte] = 
  return array[32,byte](c1) 

func unwrapComprG2*( c2: ComprG2 ): array[64,byte] = 
  return array[64,byte](c2)

func bigInt256_to_254(inp: BigInt[256]): BigInt[254] =
  var res: BigInt[254]
  res.copyTruncatedFrom(inp)
  return res

const halfPrime256 : BigInt[256] = fromHex( B, "0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3", bigEndian )
const halfPrime254 : BigInt[254] = bigInt256_to_254( halfPrime256 )
const thePrime254  : BigInt[254] = bigInt256_to_254( primeP       )

# little-endian encoding of the X coord, with bit 255 set if `Y > P/2`
func compressG1*( pt : G1 ) : ComprG1 = 
  var xbig : BigInt[254] 
  var ybig : BigInt[254]
  xbig.fromField( pt.x )
  ybig.fromField( pt.y )
  let flag : bool = bool(ybig > halfPrime254)
  var buf  : array[32,byte]
  buf.marshal(xbig, littleEndian)
  if (flag):
    buf[31] = bitor( buf[31] , 0x80 );
  return ComprG1(buf)

func uncompressG1*( compr1 : ComprG1 ) : Option[G1] = 
  var buf  : array[32,byte] = unwrapComprG1(compr1)
  let flag : bool = (buf[31] >= 0x80)
  buf[31] = bitand( buf[31] , 0x7f )
  var xbig : BigInt[254]
  xbig.unmarshal(buf, littleEndian)
  if bool(xbig >= thePrime254):
    return none(G1)
  else:
    var x : Fp[BN254_Snarks]
    var y : Fp[BN254_Snarks]
    x.fromBig(xbig)
    y = x*x*x + theCoeffB
    let ok = bool( sqrt.sqrt_if_square_vartime(y) )
    if ok:
      var ybig : BigInt[254]
      ybig.fromField( y )
      let switch = bool(ybig > halfPrime254) xor flag
      if switch:
        y.neg()
      let g1 = unsafeMkG1(x,y)
      return some(g1)
    else:
      return none(G1)

#---------------------------------------

# little-endian encoding of the X coord (real and imaginary components, in that order), 
# with the last, 511-th bit set if `Y_imag > P/2`
func compressG2*( pt : G2 ) : ComprG2 = 
  var x_real_big : BigInt[254] 
  var x_imag_big : BigInt[254] 
  var y_imag_big : BigInt[254]
  x_real_big.fromField( pt.x.coords[0] )
  x_imag_big.fromField( pt.x.coords[1] )
  y_imag_big.fromField( pt.y.coords[1] )
  let flag : bool = bool(y_imag_big > halfPrime254)
  var buf_real : array[32,byte]
  var buf_imag : array[32,byte]
  marshal(buf_real , x_real_big , littleEndian)
  marshal(buf_imag , x_imag_big , littleEndian)
  var buf: array[64,byte]
  buf[ 0..31] = buf_real
  buf[32..63] = buf_imag
  if (flag):
    buf[63] = bitor( buf[63] , 0x80 );
  return ComprG2(buf)

func uncompressG2*( compr2 : ComprG2 ) : Option[G2] = 
  var buf  : array[64,byte] = unwrapComprG2(compr2)
  let flag : bool = (buf[63] >= 0x80)
  buf[63] = bitand( buf[63] , 0x7f )
  var x_big_real : BigInt[254]
  var x_big_imag : BigInt[254]
  unmarshal(x_big_real , buf         , littleEndian)
  unmarshal(x_big_imag , buf[32..63] , littleEndian)
  if bool(x_big_real >= thePrime254) or bool(x_big_imag >= thePrime254):
    return none(G2)
  else:
    var x_real : Fp[BN254_Snarks]
    var x_imag : Fp[BN254_Snarks]
    var y      : Fp2[BN254_Snarks]
    x_real.fromBig(x_big_real)
    x_imag.fromBig(x_big_imag)
    let x: Fp2[BN254_Snarks] = mkFp2( x_real, x_imag )
    y = x*x*x + twistCoeffB
    let ok = bool( sqrt2.sqrt_if_square(y) )
    if ok:
      var y_big_imag : BigInt[254]
      y_big_imag.fromField( y.coords[1] )
      let switch = bool(y_big_imag > halfPrime254) xor flag
      if switch:
        y.neg()
      let g2 = unsafeMkG2(x,y)
      return some(g2)
    else:
      return none(G2)

#===============================================================================

func addG1*(p,q: G1): G1 =
  var r, x, y : ProjG1
  prj.fromAffine(x, p)
  prj.fromAffine(y, q)
  prj.sum(r, x, y)
  var s : G1
  prj.affine(s, r)
  return s

#---------------------------------------

func addG2*(p,q: G2): G2 =
  var r, x, y : ProjG2
  prj.fromAffine(x, p)
  prj.fromAffine(y, q)
  prj.sum(r, x, y)
  var s : G2
  prj.affine(s, r)
  return s

func negG1*(p: G1): G1 =
  var r : G1 = p
  neg(r)
  return r

func negG2*(p: G2): G2 =
  var r : G2 = p
  neg(r)
  return r

#---------------------------------------

func `+`*(p,q: G1): G1 = addG1(p,q)
func `+`*(p,q: G2): G2 = addG2(p,q)

func `+=`*(p: var G1, q: G1) =    p = addG1(p,q)
func `+=`*(p: var G2, q: G2) =    p = addG2(p,q)

func `-=`*(p: var G1, q: G1) =    p = addG1(p,negG1(q))
func `-=`*(p: var G2, q: G2) =    p = addG2(p,negG2(q))

#-------------------------------------------------------------------------------
#
# (affine) scalar multiplication
#

func `**`*( coeff: Fr , point: G1 ) : G1 =
  var q : ProjG1
  prj.fromAffine( q , point )
  scl.scalarMul_vartime(  q , coeff.toBig() )
  var r : G1
  prj.affine( r, q )
  return r

func `**`*( coeff: Fr , point: G2 ) : G2 =
  var q : ProjG2
  prj.fromAffine( q , point )
  scl.scalarMul_vartime(  q , coeff.toBig() )
  var r : G2
  prj.affine( r, q )
  return r

#-------------------

func `**`*( coeff: BigInt , point: G1 ) : G1 =
  var q : ProjG1
  prj.fromAffine( q , point )
  scl.scalarMul_vartime(  q , coeff )
  var r : G1
  prj.affine( r, q )
  return r

func `**`*( coeff: BigInt , point: G2 ) : G2 =
  var q : ProjG2
  prj.fromAffine( q , point )
  scl.scalarMul_vartime(  q , coeff )
  var r : G2
  prj.affine( r, q )
  return r

#-------------------------------------------------------------------------------

func pairing* (p: G1, q: G2) : Fp12[BN254_Snarks] =
  var t : Fp12[BN254_Snarks]
  ate.pairing_bn( t, p, q )
  return t

#-------------------------------------------------------------------------------

#[
proc sanityCheckGroupGens*() =
  echo( "gen1 on the curve        = ",  checkCurveEqG1(gen1.x,gen1.y) )
  echo( "gen2 on the curve        = ",  checkCurveEqG2(gen2.x,gen2.y) )
  echo( "gen2 is in the subgroup  = ", checkSubgroupG2(gen2.x,gen2.y) )

  let primeR : BigInt[254] = fromHex( BigInt[254], "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", bigEndian )
  echo( "order of gen1 is R  = ", (not bool(isNeutral(gen1))) and bool(isNeutral(primeR ** gen1)) )
  echo( "order of gen2 is R  = ", (not bool(isNeutral(gen2))) and bool(isNeutral(primeR ** gen2)) )

#
# the point (computed via Sage)
#
#   (2 : 2237046587054574173616397632856518880513033439888792180868262182050662989363*u + 10894412225134874879786325788974416805327887441035008073952212076423500941133 : 1)
#
# should be on the curve but not in the subgroup
#
proc sanityCheckInSubgroupG2*() = 
  let pt2_x1  = fromHex(Fp[BN254_Snarks], "0x2")
  let pt2_xu  = fromHex(Fp[BN254_Snarks], "0x0")
  let pt2_y1  = fromHex(Fp[BN254_Snarks], "0x181604d0560080401c08b557815482553e278257d98100d193a011c42782474d")
  let pt2_yu  = fromHex(Fp[BN254_Snarks], "0x04f21f9d99cc25f694cf22ff70dc0ac4692e7a721b725dc454a217f04bd03e33")
  let pt2_x   = mkFp2( pt2_x1, pt2_xu )
  let pt2_y   = mkFp2( pt2_y1, pt2_yu )
  echo("pt2 is on the curve    (should be true )   = " ,  checkCurveEqG2(pt2_x, pt2_y) )
  echo("pt2 is in the subgroup (should be false)   = " , checkSubgroupG2(pt2_x, pt2_y) )

]#

#-------------------------------------------------------------------------------
