
#
# power-of-two sized multiplicative FFT domains in the scalar field
#

import constantine/math/arithmetic
import constantine/math/io/io_fields
#import constantine/math/io/io_bigints
import constantine/named/properties_fields

import groth16/bn128
import groth16/misc

#-------------------------------------------------------------------------------

type 

  Domain* = object
    domainSize*    : int                # `N = 2^n`
    logDomainSize* : int                # `n = log2(N)`
    domainGen*     : Fr[BN254_Snarks]   # `g`
    invDomainGen*  : Fr[BN254_Snarks]   # `g^-1`
    invDomainSize* : Fr[BN254_Snarks]   # `1/n`

  Subgroup* = object
    bigDomain*   : Domain
    smallDomain* : Domain

#-------------------------------------------------------------------------------

# the generator of the multiplicative subgroup with size `2^28`
const gen28 = fromHex( Fr[BN254_Snarks], "0x2a3c09f0a58a7e8500e0a7eb8ef62abc402d111e41112ed49bd61b6e725b19f0" )

func createDomain*(size: int): Domain = 
  let log2 = ceilingLog2(size)
  assert( (1 shl log2) == size , "domain must have a power-of-two size" )

  let expo : uint = 1'u shl (28 - log2)
  let gen  = smallPowFr(gen28, expo)

  let halfSize = size div 2
  let a = smallPowFr(gen, size    )
  let b = smallPowFr(gen, halfSize)
  assert(     bool(a == oneFr) , "domain generator sanity check /A" )
  assert( not bool(b == oneFr) , "domain generator sanity check /B" )

  return Domain( domainSize:    size
               , logDomainSize: log2
               , domainGen:     gen
               , invDomainGen:  invFr(gen)
               , invDomainSize: invFr(intToFr(size)) 
               )

#-------------------------------------------------------------------------------

func enumerateDomain*(D: Domain): seq[Fr[BN254_Snarks]] =
  var xs = newSeq[Fr[BN254_Snarks]](D.domainSize)
  var g = oneFr
  for i in 0..<D.domainSize:
    xs[i] = g
    g *= D.domainGen
  return xs

#-------------------------------------------------------------------------------

func subgroupIndex*(sg: Subgroup): int = 
  return (sg.bigDomain.domainSize div sg.smallDomain.domainSize)

func createSubgroup*(D: Domain, K: int): Subgroup =
  let N = D.domainSize
  assert( K >= 1 )
  assert( K <= N )
  assert( (N mod K) == 0 )
  return Subgroup( bigDomain: D, smallDomain: createDomain(K) )

func selectOnSubgroup*[T]( sg: Subgroup , xs: seq[T] ): seq[T] = 
  let N = sg.bigDomain.domainSize
  let K = sg.smallDomain.domainSize
  assert( xs.len == N )
  var ys : seq[T] = newSeq[T]( K )
  let J = N div K
  for i in 0..<K:
    ys[i] = xs[i*J]
  return ys

func liesInSubgroup*( sg: Subgroup, mask: seq[bool] ): bool = 
  let N = sg.bigDomain.domainSize
  let K = sg.smallDomain.domainSize
  let J = N div K

  assert( N == mask.len )

  var ok = true
  for (i,b) in mask.pairs:
    if b:
      ok = ok and ((i mod J) == 0)

  return ok

#-------------------------------------------------------------------------------
