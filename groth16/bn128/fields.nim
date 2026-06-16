
#
# the prime fields Fp and Fr with sizes
#
#   p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
#   r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
#

import sugar
import std/bitops
import std/options
# import std/sequtils

import constantine/platforms/abstractions
import constantine/math/arithmetic
import constantine/math/io/io_fields
import constantine/math/io/io_bigints
import constantine/named/properties_fields      as tff
import constantine/math/extension_fields/towers as ext

#-------------------------------------------------------------------------------

type B*    = BigInt[256]
type B254* = BigInt[254]       # some more recent constantine idiosyncrasies.....

func mkFp2* (i: Fp[BN254_Snarks], u: Fp[BN254_Snarks]) : Fp2[BN254_Snarks] =
  let c : array[2, Fp[BN254_Snarks]] = [i,u]
  return ext.QuadraticExt[Fp[BN254_Snarks]]( coords: c )

#-------------------------------------------------------------------------------

const primeR*  : B = fromHex( B, "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", bigEndian )
const primeP*  : B = fromHex( B, "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47", bigEndian )

const primeR_B254* : BigInt[254] = fromHex( BigInt[254] , "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", bigEndian )

#-------------------------------------------------------------------------------

const zeroFp*   = fromHex( Fp[BN254_Snarks], "0x00" )
const zeroFr*   = fromHex( Fr[BN254_Snarks], "0x00" )

const oneFp*    = fromHex( Fp[BN254_Snarks], "0x01" )
const oneFr*    = fromHex( Fr[BN254_Snarks], "0x01" )

const twoFp*    = fromHex( Fp[BN254_Snarks], "0x02" )
const twoFr*    = fromHex( Fr[BN254_Snarks], "0x02" )

const zeroFp2* = mkFp2( zeroFp, zeroFp )
const oneFp2*  = mkFp2( oneFp , zeroFp )

const minusOneFp* = fromHex( Fp[BN254_Snarks], "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46" )
const minusOneFr* = fromHex( Fr[BN254_Snarks], "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000" )

const oneHalfFr*  = fromHex( Fr[BN254_Snarks], "0x183227397098d014dc2822db40c0ac2e9419f4243cdcb848a1f0fac9f8000001" )

#-------------------------------------------------------------------------------

func intToB*(a: uint): B =
  var y : B
  y.setUint(a)
  return y

func intToB254*(a: uint): B254 =
  var y : B254
  y.setUint(a)
  return y

func intToFp*(a: int): Fp[BN254_Snarks] =
  var y : Fp[BN254_Snarks]
  y.fromInt(a)
  return y

func intToFr*(a: int): Fr[BN254_Snarks] =
  var y : Fr[BN254_Snarks]
  y.fromInt(a)
  return y

# -- seriously... 
# func frToB*(x: Fr[BN254_Snarks]) : B = 
#   var b : B
#   b.fromField(x)
#   return b

func frToB254*(x: Fr[BN254_Snarks]) : B254 = 
  var b : BigInt[254]
  b.fromField(x)
  return b

#-------------------------------------------------------------------------------

# if the field element is a small number (either positive or negative!) then we return that
# useful for printing equation coefficients
#
func mbIsSmallFr*(x: Fr[BN254_Snarks]): Option[int] =
  let bpos  : B254 = frToB254(x)
  let thres : B254 = intToB254(1000000000)
  if bool(bpos < thres):
    return some(bpos.limbs[0].int)
  else:
    var y = zeroFr
    y -= x
    let bneg : B254 = frToB254(y)
    if bool(bneg < thres):
      return some(-bneg.limbs[0].int)
    else:
      return none(int)

#-------------------------------------------------------------------------------

func isZeroFp* (x: Fp ): bool = bool(isZero(x))
func isZeroFp2*(x: Fp2): bool = bool(isZero(x))
func isZeroFr* (x: Fr ): bool = bool(isZero(x))

func isEqualFp* (x, y: Fp ): bool = bool(x == y)
func isEqualFp2*(x, y: Fp2): bool = bool(x == y)
func isEqualFr* (x, y: Fr ): bool = bool(x == y)

func `===`*(x, y: Fp ): bool = isEqualFp(x,y)
func `===`*(x, y: Fp2): bool = isEqualFp2(x,y)
func `===`*(x, y: Fr ): bool = isEqualFr(x,y)

#-------------------

func isEqualFpSeq*(xs, ys: seq[Fp]): bool =
  let n = xs.len
  assert( n == ys.len )
  var b = true
  for i in 0..<n:
    if not bool(xs[i] == ys[i]):
      b = false
      break
  return b

func isEqualFrSeq*(xs, ys: seq[Fr]): bool =
  let n = xs.len
  assert( n == ys.len )
  var b = true
  for i in 0..<n:
    if not bool(xs[i] == ys[i]):
      b = false
      break
  return b

func `===`*(xs, ys: seq[Fp]): bool = isEqualFpSeq(xs,ys)
func `===`*(xs, ys: seq[Fr]): bool = isEqualFrSeq(xs,ys)

#-------------------------------------------------------------------------------

func `+`*[n](x, y: BigInt[n] ): BigInt[n] = ( var z : BigInt[n] = x ; z += y ; return z )
func `+`*(x, y: Fp ): Fp  =  ( var z : Fp  = x ; z += y ; return z )
func `+`*(x, y: Fp2): Fp2 =  ( var z : Fp2 = x ; z += y ; return z )
func `+`*(x, y: Fr ): Fr  =  ( var z : Fr  = x ; z += y ; return z )

func `-`*[n](x, y: BigInt[n] ): BigInt[n] = ( var z : BigInt[n] = x ; z -= y ; return z )
func `-`*(x, y: Fp ): Fp  =  ( var z : Fp  = x ; z -= y ; return z )
func `-`*(x, y: Fp2): Fp2 =  ( var z : Fp2 = x ; z -= y ; return z )
func `-`*(x, y: Fr ): Fr  =  ( var z : Fr  = x ; z -= y ; return z )

func `*`*(x, y: Fp ): Fp  =  ( var z : Fp  = x ; z *= y ; return z )
func `*`*(x, y: Fp2): Fp2 =  ( var z : Fp2 = x ; z *= y ; return z )
func `*`*(x, y: Fr ): Fr  =  ( var z : Fr  = x ; z *= y ; return z )

func negFp* (y: Fp ): Fp  =  ( var z : Fp  = zeroFp  ; z -= y ; return z )
func negFp2*(y: Fp2): Fp2 =  ( var z : Fp2 = zeroFp2 ; z -= y ; return z )
func negFr* (y: Fr ): Fr  =  ( var z : Fr  = zeroFr  ; z -= y ; return z )

func invFp*(y: Fp): Fp =  ( var z : Fp = y ; inv(z) ; return z )
func invFr*(y: Fr): Fr =  ( var z : Fr = y ; inv(z) ; return z )

func squareFp* (y: Fp):  Fp  =  ( var z : Fp  = y ; square(z) ; return z )
func squareFp2*(y: Fp2): Fp2 =  ( var z : Fp2 = y ; square(z) ; return z )
func squareFr* (y: Fr):  Fr  =  ( var z : Fr  = y ; square(z) ; return z )

func `/`*(x, y: Fp ): Fp  =  ( var z : Fp  = x ; z *= invFr(y) ; return z )
func `/`*(x, y: Fp2): Fp2 =  ( var z : Fp2 = x ; z *= invFr(y) ; return z )
func `/`*(x, y: Fr ): Fr  =  ( var z : Fr  = x ; z *= invFr(y) ; return z )

func divBy2Fp* (x: Fp ): Fp = ( var z : Fp  = x ; div2(z) ; return z )
func divBy2Fp2*(x: Fp2): Fp2= ( var z : Fp2 = x ; div2(z) ; return z )
func divBy2Fr* (x: Fr ): Fr = ( var z : Fr  = x ; div2(z) ; return z )

# template/generic instantiation of `pow_vartime` from here
# /Users/bkomuves/.nimble/pkgs/constantine-0.0.1/constantine/math/arithmetic/finite_fields.nim(389, 7) template/generic instantiation of `fieldMod` from here
# /Users/bkomuves/.nimble/pkgs/constantine-0.0.1/constantine/math/config/curves_prop_field_derived.nim(67, 5) Error: undeclared identifier: 'getCurveOrder'
# ...
func smallPowFr*(base: Fr, expo: uint): Fr =
  var a : Fr   = oneFr
  var s : Fr   = base
  var e : uint = expo
  while (e > 0):
    if bitand(e,1) > 0: a *= s
    e = (e shr 1)
    square(s)
  return a

func smallPowFr*(base: Fr, expo: int): Fr =
  if expo >= 0:
    return smallPowFr( base, uint(expo) )
  else:
    return smallPowFr( invFr(base) , uint(-expo) )

#-------------------------------------------------------------------------------

func deltaFr*(i, j: int) : Fr[BN254_Snarks] =
  return (if (i == j): oneFr else: zeroFr)

#-------------------------------------------------------------------------------

# Montgomery batch inversion
func batchInverseFr*( xs: seq[Fr[BN254_Snarks]] ) : seq[Fr[BN254_Snarks]] =
  let n = xs.len
  assert(n>0)
  var us : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]](n+1)
  var a = xs[0]
  us[0] = oneFr
  us[1] = a
  for i in 1..<n: ( a *= xs[i] ; us[i+1] = a )
  var vs : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]](n)
  vs[n-1] = invFr( us[n] )
  for i in countdown(n-2,0): vs[i] = vs[i+1] * xs[i+1]
  return collect( newSeq, (for i in 0..<n: us[i]*vs[i] ) )


#-------------------------------------------------------------------------------

