
#
# Number-theoretic transform 
# (that is, FFT for polynomials over finite fields)
#

#-------------------------------------------------------------------------------

import constantine/math/arithmetic
import constantine/math/io/io_fields
import constantine/named/properties_fields

import groth16/bn128
import groth16/math/domain

#-------------------------------------------------------------------------------

func forwardNTT_worker( m: int
                      , srcStride: int
                      , gpows: seq[Fr[BN254_Snarks]]
                      , src:     seq[Fr[BN254_Snarks]] , srcOfs: int
                      , buf: var seq[Fr[BN254_Snarks]] , bufOfs: int
                      , tgt: var seq[Fr[BN254_Snarks]] , tgtOfs: int ) =
  case m 

    of 0: 
      tgt[tgtOfs] = src[srcOfs]
  
    of 1:
      tgt[tgtOfs  ] = src[srcOfs] + src[srcOfs+srcStride]
      tgt[tgtOfs+1] = src[srcOfs] - src[srcOfs+srcStride]

    else:
      let N     : int =  1 shl  m  
      let halfN : int =  1 shl (m-1)  
      forwardNTT_worker( m-1
                       , srcStride shl 1
                       , gpows
                       , src , srcOfs
                       , buf , bufOfs + N
                       , buf , bufOfs )
      forwardNTT_worker( m-1
                       , srcStride shl 1
                       , gpows
                       , src , srcOfs + srcStride
                       , buf , bufOfs + N
                       , buf , bufOfs + halfN )
      for j in 0..<halfN:
        let y = gpows[j*srcStride] * buf[bufOfs+j+halFN]
        tgt[tgtOfs+j      ] = buf[bufOfs+j] + y
        tgt[tgtOfs+j+halfN] = buf[bufOfs+j] - y

#---------------------------------------

# forward number-theoretical transform (corresponds to polynomial evaluation)
func forwardNTT*(src: seq[Fr[BN254_Snarks]], D: Domain): seq[Fr[BN254_Snarks]] =
  assert( D.domainSize == (1 shl D.logDomainSize) , "domain must have a power-of-two size" )
  assert( D.domainSize == src.len , "input must have the same size as the domain" )
  var buf = newSeq[Fr[BN254_Snarks] ]( 2 * D.domainSize )
  var tgt = newSeq[Fr[BN254_Snarks]](     D.domainSize )

  # precalc powers of gen
  let N     = D.domainSize
  let halFN = N div 2
  var gpows = newSeq[Fr[BN254_Snarks]]( halFN )
  var x     = oneFr
  let gen   = D.domainGen
  for i in 0..<halfN:
    gpows[i] = x
    x *= gen

  forwardNTT_worker( D.logDomainSize
                   , 1
                   , gpows
                   , src , 0
                   , buf , 0
                   , tgt , 0 )
  return tgt

# pads the input with zeros to get a pwoer of two size
# TODO: optimize the FFT so that it doesn't do the multiplications with zeros 
func extendAndForwardNTT*(src: seq[Fr[BN254_Snarks]], D: Domain): seq[Fr[BN254_Snarks]] =
  let n = src.len
  let N = D.domainSize 
  assert( n <= N )
  if n == N:
    return forwardNTT(src, D)
  else:
    var padded = newSeq[Fr[BN254_Snarks]]( N )
    for i in 0..<n: padded[i] = src[i]
    # for i in n..<N: padded[i] = zeroFr 
    return forwardNTT(padded, D)

#-------------------------------------------------------------------------------

# unscaled!
func inverseNTT_worker( m: int
                      , tgtStride: int
                      , gpows: seq[Fr[BN254_Snarks]]
                      , src:     seq[Fr[BN254_Snarks]] , srcOfs: int
                      , buf: var seq[Fr[BN254_Snarks]] , bufOfs: int
                      , tgt: var seq[Fr[BN254_Snarks]] , tgtOfs: int ) =
  case m 

    of 0: 
      tgt[tgtOfs] = src[srcOfs]
  
    of 1:
      tgt[tgtOfs          ] = ( src[srcOfs] + src[srcOfs+1] ) 
      tgt[tgtOfs+tgtStride] = ( src[srcOfs] - src[srcOfs+1] ) 

    else:
      let N     : int =  1 shl  m  
      let halfN : int =  1 shl (m-1)  

      for j in 0..<halfN:
        buf[bufOfs+j      ] = ( src[srcOfs+j] + src[srcOfs+j+halfN] ) 
        buf[bufOfs+j+halfN] = ( src[srcOfs+j] - src[srcOfs+j+halfN] ) * gpows[ j*tgtStride ]

      inverseNTT_worker( m-1
                       , tgtStride shl 1
                       , gpows
                       , buf , bufOfs
                       , buf , bufOfs + N
                       , tgt , tgtOfs )
      inverseNTT_worker( m-1
                       , tgtStride shl 1
                       , gpows
                       , buf , bufOfs + halfN
                       , buf , bufOfs + N
                       , tgt , tgtOfs + tgtStride )

#---------------------------------------

# inverse number-theoretical transform (corresponds to polynomial interpolation)
func toggableInverseNTT(src: seq[Fr[BN254_Snarks]], D: Domain, do_rescale: bool): seq[Fr[BN254_Snarks]] =
  assert( D.domainSize == (1 shl D.logDomainSize) , "domain must have a power-of-two size" )
  assert( D.domainSize == src.len , "input must have the same size as the domain" )
  var buf = newSeq[Fr[BN254_Snarks]]( 2 * D.domainSize )
  var tgt = newSeq[Fr[BN254_Snarks]](     D.domainSize )

  # precalc 1/2 times powers of gen^-1
  let N     = D.domainSize
  let halFN = N div 2
  var gpows = newSeq[Fr[BN254_Snarks]]( halFN )
  var x     = oneFr
  let ginv  = invFr( D.domainGen )
  for i in 0..<halfN:
    gpows[i] = x
    x *= ginv

  inverseNTT_worker( D.logDomainSize
                   , 1
                   , gpows
                   , src , 0
                   , buf , 0
                   , tgt , 0 )

  if do_rescale:
    var invN : Fr[BN254_Snarks]
    invN.fromInt(N)
    invN.inv() 
    for i in 0..<N:
      tgt[i] *= invN

  return tgt

#-------------------------------------------------------------------------------

# inverse number-theoretical transform (corresponds to polynomial interpolation)
func inverseNTT*(src: seq[Fr[BN254_Snarks]], D: Domain): seq[Fr[BN254_Snarks]] =
  toggableInverseNTT(src, D, true)

# inverse number-theoretical transform, without the 1/N rescaling
func unscaledInverseNTT*(src: seq[Fr[BN254_Snarks]], D: Domain): seq[Fr[BN254_Snarks]] =
  toggableInverseNTT(src, D, false)

#-------------------------------------------------------------------------------
