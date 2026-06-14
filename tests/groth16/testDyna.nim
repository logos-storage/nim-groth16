
{.used.} 

import std/unittest

# import constantine/math/arithmetic
# import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/rnd

import groth16/dynamic/types
import groth16/dynamic/preprocess
#import groth16/dynamic/setup

#-------------------------------------------------------------------------------

suite "dynamic proof tests":

  let N : int = 128

  let tau   : F = randFr()
  let delta : F = randFr()

  test "projection term V1 calculation":

    check testProjectionElementsV1( N , tau , delta )

#-------------------------------------------------------------------------------
