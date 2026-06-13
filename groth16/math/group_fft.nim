
#
# FFT for groups elements
#

#-------------------------------------------------------------------------------

import constantine/math/arithmetic
import constantine/math/io/io_fields
import constantine/named/properties_fields

# import constantine/math/elliptic/ec_shortweierstrass_affine     as aff 
import constantine/math/elliptic/ec_shortweierstrass_projective as prj 
import constantine/math/elliptic/ec_scalar_mul_vartime          as scl 

import groth16/bn128
import groth16/bn128/curves
import groth16/math/domain

#-------------------------------------------------------------------------------

func forward_FFT_worker( m: int
                       , srcStride: int
                       , gpows: seq[Fr[BN254_Snarks]]
                       , src:     seq[ProjG1] , srcOfs: int
                       , buf: var seq[ProjG1] , bufOfs: int
                       , tgt: var seq[ProjG1] , tgtOfs: int ) =
  case m 

    of 0: 
      tgt[tgtOfs] = src[srcOfs]
  
    of 1:
      tgt[tgtOfs  ] = src[srcOfs] ; tgt[tgtOfs  ] += src[srcOfs+srcStride]
      tgt[tgtOfs+1] = src[srcOfs] ; tgt[tgtOfs+1] -= src[srcOfs+srcStride]

    else:
      let N     : int =  1 shl  m  
      let halfN : int =  1 shl (m-1)  
      forward_FFT_worker( m-1
                        , srcStride shl 1
                        , gpows
                        , src , srcOfs
                        , buf , bufOfs + N
                        , buf , bufOfs )
      forward_FFT_worker( m-1
                        , srcStride shl 1
                        , gpows
                        , src , srcOfs + srcStride
                        , buf , bufOfs + N
                        , buf , bufOfs + halfN )
      for j in 0..<halfN:
        var y = buf[bufOfs+j+halFN]
        y.scalarMul_vartime( gpows[j*srcStride] ) 
        tgt[tgtOfs+j      ] = buf[bufOfs+j] ; tgt[tgtOfs+j      ] += y
        tgt[tgtOfs+j+halfN] = buf[bufOfs+j] ; tgt[tgtOfs+j+halfN] -= y

#---------------------------------------

# forward number-theoretical transform (corresponds to polynomial evaluation)
func forwardGroupFFT*(affSrc: seq[AffG1], D: Domain): seq[AffG1] =
  let N = D.domainSize
  assert( N == (1 shl D.logDomainSize) , "domain must have a power-of-two size"        )
  assert( N == affSrc.len              , "input must have the same size as the domain" )

  var src = newSeq[ProjG1](     N )
  var buf = newSeq[ProjG1]( 2 * N )
  var tgt = newSeq[ProjG1](     N )

  for i in 0..<N:
    src[i].fromAffine(affSrc[i])

  # precalc powers of gen
  let halFN = N div 2
  var gpows = newSeq[Fr[BN254_Snarks]]( halFN )
  var x     = oneFr
  let gen   = D.domainGen
  for i in 0..<halfN:
    gpows[i] = x
    x *= gen

  forward_FFT_worker( D.logDomainSize
                    , 1
                    , gpows
                    , src , 0
                    , buf , 0
                    , tgt , 0 )

  var affTgt = newSeq[AffG1]( N )
  for i in 0..<N:
    affTgt[i].affine(tgt[i])

  return affTgt

#-------------------------------------------------------------------------------

# unscaled!
func inverse_FFT_worker( m: int
                       , tgtStride: int
                       , gpows: seq[Fr[BN254_Snarks]]
                       , src:     seq[ProjG1] , srcOfs: int
                       , buf: var seq[ProjG1] , bufOfs: int
                       , tgt: var seq[ProjG1] , tgtOfs: int ) =
  case m 

    of 0: 
      tgt[tgtOfs] = src[srcOfs]
  
    of 1:
      tgt[tgtOfs          ] = src[srcOfs] ; tgt[tgtOfs          ] += src[srcOfs+1] 
      tgt[tgtOfs+tgtStride] = src[srcOfs] ; tgt[tgtOfs+tgtStride] -= src[srcOfs+1] 

    else:
      let N     : int =  1 shl  m  
      let halfN : int =  1 shl (m-1)  

      for j in 0..<halfN:
        buf[bufOfs+j      ] = src[srcOfs+j] ; buf[bufOfs+j      ] += src[srcOfs+j+halfN] 
        buf[bufOfs+j+halfN] = src[srcOfs+j] ; buf[bufOfs+j+halfN] -= src[srcOfs+j+halfN] 
        buf[bufOfs+j+halfN].scalarMul_vartime( gpows[ j*tgtStride ] ) 

      inverse_FFT_worker( m-1
                        , tgtStride shl 1
                        , gpows
                        , buf , bufOfs
                        , buf , bufOfs + N
                        , tgt , tgtOfs )
      inverse_FFT_worker( m-1
                        , tgtStride shl 1
                        , gpows
                        , buf , bufOfs + halfN
                        , buf , bufOfs + N
                        , tgt , tgtOfs + tgtStride )

#---------------------------------------

# inverse number-theoretical transform (with the 1/N rescaling toggable)
func toggableInverseGroupFFT(affSrc: seq[AffG1], D: Domain, do_rescale: bool): seq[AffG1] =
  let N = D.domainSize
  assert( N == (1 shl D.logDomainSize) , "domain must have a power-of-two size"        )
  assert( N == affSrc.len              , "input must have the same size as the domain" )

  var src = newSeq[ProjG1](     N )
  var buf = newSeq[ProjG1]( 2 * N )
  var tgt = newSeq[ProjG1](     N )

  for i in 0..<N:
    src[i].fromAffine(affSrc[i])

  # precalc times powers of gen^-1
  let halFN = N div 2
  var gpows = newSeq[Fr[BN254_Snarks]]( halFN )
  var x     = oneFr
  let ginv  = invFr( D.domainGen )
  for i in 0..<halfN:
    gpows[i] = x
    x *= ginv

  inverse_FFT_worker( D.logDomainSize
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
      tgt[i].scalarMul_vartime( invN )

  var affTgt = newSeq[AffG1]( N )
  for i in 0..<N:
    affTgt[i].affine(tgt[i])

  return affTgt

#---------------------------------------

func inverseGroupFFT*(src: seq[AffG1], D: Domain): seq[AffG1] =
  toggableInverseGroupFFT(src , D , true)

# inverse number-theoretical transform (without the 1/N rescaling)
func unscaledInverseGroupFFT*(src: seq[AffG1], D: Domain): seq[AffG1] =
  toggableInverseGroupFFT(src , D , false)

#-------------------------------------------------------------------------------
