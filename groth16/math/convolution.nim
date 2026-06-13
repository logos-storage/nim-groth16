
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

  let us =      forwardNTT(xs, D) 
  let hs = forwardGroupFFT(gs, D)
  let rs = pointwiseScaleG1(us, hs)

  return inverseGroupFFT(rs, D)

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
