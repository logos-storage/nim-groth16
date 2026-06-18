
import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/curves
import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/ntt
import groth16/math/group_fft

#-------------------------------------------------------------------------------

type F = Fr[BN254_Snarks]

func fieldConvolution*(xs: seq[F], ys: seq[F]): seq[F] = 
  let N = xs.len
  assert( N == ys.len )
  let D = createDomain( N )

  let us = forwardNTT(xs, D) 
  let vs = forwardNTT(ys, D)
  let ws = pointwiseProdFr(us, vs)

  return inverseNTT(ws, D)

#-------------------------------------------------------------------------------

func groupConvolution*(xs: seq[F], gs: seq[G1]): seq[G1] =
  let N = xs.len
  assert( N == gs.len )
  let D = createDomain( N )

  let us =       forwardNTT(xs, D) 
  let hs =  forwardGroupFFT(gs, D)
  let rs = pointwiseScaleG1(us, hs)

  return inverseGroupFFT(rs, D)

#-------------------------------------------------------------------------------

func fieldConvolutionOnSubgroup*(sg: Subgroup, xs: seq[F], ys: seq[F]): seq[F] = 
  let N = xs.len
  assert( N == ys.len )
  assert( N == sg.bigDomain.domainSize )
  let D = sg.bigDomain
  let K = sg.smallDomain.domainSize
  let ell = N div K
  let fell    : F = intToFr(ell)
  let invfell : F = invFr(fell)

  let us = forwardNTT(xs, D) 
  let vs = forwardNTT(ys, D)

  var ws: seq[F] = newSeq[F]( K )
  for k in 0..<K:
    var sum: F = zeroFr
    for i in 0..<ell:
      sum += us[k + i*K] * vs[k + i*K]
    ws[k] = sum * invfell

  return inverseNTT(ws, sg.smallDomain)

#---------------------------------------

# only computes the result in a subgroup
func groupConvolutionOnSubgroup*(sg: Subgroup, xs: seq[F], gs: seq[G1]): seq[G1] =
  let N = xs.len
  assert( N == gs.len )
  assert( N == sg.bigDomain.domainSize )
  let D = sg.bigDomain
  let K = sg.smallDomain.domainSize
  let ell = N div K
  let fell    : F = intToFr(ell)
  let invfell : F = invFr(fell)

  let us =      forwardNTT(xs, D) 
  let hs = forwardGroupFFT(gs, D)
  
  var rs      : seq[G1] = newSeq[G1]( K   )
  var miniPts : seq[G1] = newSeq[G1]( ell )
  var miniCfs : seq[F]  = newSeq[F] ( ell )
  for k in 0..<K:
    for i in 0..<ell:
      miniCfs[i] = us[k + K*i] * invfell
      miniPts[i] = hs[k + K*i]
    rs[k] = msmConstantineG1( miniCfs , miniPts )

  return inverseGroupFFT(rs, sg.smallDomain)

#---------------------------------------

proc testFieldConvolutionOnSubgroup*(N: int, K: int): bool = 
  let sg = createSubgroup( createDomain(N) , K )

  let xs = randFrSeq(N)
  let ys = randFrSeq(N)

  let target = selectOnSubgroup( sg , fieldConvolution(xs, ys) )
  let smart  = fieldConvolutionOnSubgroup( sg, xs, ys )
  return isEqualFrSeq(target , smart)


proc testGroupConvolutionOnSubgroup*(N: int, K: int): bool = 
  let sg = createSubgroup( createDomain(N) , K )

  let xs = randFrSeq(N)
  let ps = randG1Seq(N)

  let target = selectOnSubgroup( sg , groupConvolution(xs, ps) )
  let smart  = groupConvolutionOnSubgroup( sg, xs, ps )
  return isEqualG1Seq(target , smart)

#-------------------------------------------------------------------------------
# naive reference implementations for testing

func naiveFieldConvolution*(xs: seq[F], ys: seq[F]): seq[F] = 
  let N = xs.len
  assert( N == ys.len )

  var zs: seq[F] = newSeq[F]( N )
  for k in 0..<N:
    var acc: F = zeroFr 
    for i in 0..<N:
      acc += xs[ (k-i+N) mod N ] * ys[ i ]
    zs[k] = acc

  return zs  

func naiveGroupConvolution*(xs: seq[F], gs: seq[G1]): seq[G1] = 
  let N = xs.len
  assert( N == gs.len )

  var hs: seq[G1] = newSeq[G1]( N )
  for k in 0..<N:
    var acc: G1 = infG1 
    for i in 0..<N:
      acc += xs[ (k-i+N) mod N ] ** gs[ i ]
    hs[k] = acc

  return hs

#-------------------------------------------------------------------------------
