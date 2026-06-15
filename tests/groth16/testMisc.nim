
{.used.} 

import std/unittest
import std/math

import groth16/misc

#-------------------------------------------------------------------------------

suite "misc":

  test "floorLog2":

    var ok: bool = true
    for i in 1..100:
      let x = float64(i)
      let lhs = floorLog2(i) 
      let rhs = floor(log2(x)).int
      # echo $lhs & " = " & $rhs 
      ok = ok and ( lhs == rhs )
    check ok

  test "ceilingLog2":

    var ok: bool = true
    for i in 1..100:
      let x = float64(i)
      let lhs = ceilingLog2(i) 
      let rhs = ceil(log2(x)).int
      # echo $lhs & " = " & $rhs 
      ok = ok and ( lhs == rhs )
    check ok

#-------------------------------------------------------------------------------
