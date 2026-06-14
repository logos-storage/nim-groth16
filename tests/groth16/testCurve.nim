
{.used.} 

import std/unittest

import constantine/math/io/io_bigints
import constantine/math/arithmetic
import constantine/math/io/io_fields
import constantine/named/properties_fields 
import constantine/math/extension_fields/towers 

import groth16/bn128
import groth16/bn128/fields
import groth16/bn128/curves
import groth16/bn128/rnd

#-------------------------------------------------------------------------------

suite "curve arithmetic":

  test "G1 add/sub":
    let g = randG1()
    let h = randG1()
    let lhs = g - h
    let rhs = g + negG1(h)
    var x = g
    x -= h
    check ((lhs === rhs) and (x === lhs))

  test "G2 add/sub":
    let g = randG2()
    let h = randG2()
    let lhs = g - h
    let rhs = g + negG2(h)
    var x = g
    x -= h
    check ((lhs === rhs) and (x === lhs))

#-------------------------------------------------------------------------------

#
# the point (computed via Sage)
#
#   pt2 = (2 : 2237046587054574173616397632856518880513033439888792180868262182050662989363*u + 10894412225134874879786325788974416805327887441035008073952212076423500941133 : 1)
#
# should be on the curve but not in the subgroup
#

const pt2_x1  = fromHex(Fp[BN254_Snarks], "0x2")
const pt2_xu  = fromHex(Fp[BN254_Snarks], "0x0")
const pt2_y1  = fromHex(Fp[BN254_Snarks], "0x181604d0560080401c08b557815482553e278257d98100d193a011c42782474d")
const pt2_yu  = fromHex(Fp[BN254_Snarks], "0x04f21f9d99cc25f694cf22ff70dc0ac4692e7a721b725dc454a217f04bd03e33")
const pt2_x   = mkFp2( pt2_x1, pt2_xu )
const pt2_y   = mkFp2( pt2_y1, pt2_yu )

suite "curve and subgroup checks":

  test "gen1 is on the curve":
    check checkCurveEqG1(gen1.x,gen1.y) 

  test "gen1 is in the subgroup G1":
    check checkSubgroupG1(gen1.x,gen1.y) 

  test "gen2 is on the curve over Fp2":
    check checkCurveEqG2(gen2.x,gen2.y)
  
  test "gen2 is in the subgroup G2":
    check checkSubgroupG2(gen2.x,gen2.y)

  let prime254 : BigInt[254] = fromHex( BigInt[254], "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", bigEndian )

  test "order of gen1 equals to R":
    check ( (not bool(isInfG1(gen1))) and bool(isInfG1(prime254 ** gen1)) )

  test "order of gen2 equals to R":
    check ( (not bool(isInfG2(gen2))) and bool(isInfG2(prime254 ** gen2)) )

  test "pt2 is on the curve over Fp2":
    check checkCurveEqG2(pt2_x, pt2_y) 

  test "pt2 is NOT in the subgroup G2":
    check (not checkSubgroupG2(pt2_x, pt2_y))

#-------------------------------------------------------------------------------
