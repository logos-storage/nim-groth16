

{.used.} 

import std/unittest

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128/fields
import groth16/bn128/curves
import groth16/bn128/arrays
# import groth16/bn128/rnd
# import groth16/bn128/debug

import groth16/math/domain
import groth16/math/ntt
import groth16/math/group_fft
import groth16/math/convolution

#-------------------------------------------------------------------------------

suite "field NTT checks":

  let D = createDomain(128)
  let N = D.domainSize
  let xs : seq[Fr[BN254_Snarks]] = randFrSeq(N)

  test "INTT(NTT(xs)) == xs":
    let ys = forwardNTT(xs, D)
    let zs = inverseNTT(ys, D) 
    check isEqualFrSeq( xs , zs )
  
  test "NTT(INTT(xs)) == xs":
    let ys = inverseNTT(xs, D)
    let zs = forwardNTT(ys, D) 
    check isEqualFrSeq( xs , zs )

  test "unscaledINTT(NTT(xs)) == N * xs":
    let ys = forwardNTT(xs, D)
    let zs = unscaledInverseNTT(ys, D) 
    var ws = xs
    scaleFrSeqInPlace(intToFr(N) , ws)
    check isEqualFrSeq( ws , zs )

#-------------------------------------

suite "group FFT checks (for the group G1)":

  let D = createDomain(64)
  let N = D.domainSize
  let gs : seq[G1] = randG1Seq(N)

  test "IFFT(FFT(gs)) == gs":
    let hs = forwardGroupFFT(gs, D)
    let rs = inverseGroupFFT(hs, D) 
    check isEqualG1Seq( gs , rs )
  
  test "FFT(IFFT(xs)) == xs":
    let hs = inverseGroupFFT(gs, D)
    let rs = forwardGroupFFT(hs, D) 
    check isEqualG1Seq( gs , rs )

  test "FFT(unscaledIFFT(gs)) == N * gs":
    let hs = unscaledInverseGroupFFT(gs, D)
    let rs = forwardGroupFFT(hs, D) 
    var qs = gs
    scaleG1SeqInPlace(intToFr(N) , qs)
    check isEqualG1Seq( qs , rs )

#-------------------------------------------------------------------------------
  
suite "convolution checks":

  let N = 64

  let xs = randFrSeq(N)
  let ys = randFrSeq(N)
  let gs = randG1Seq(N)

  test "FFT field convolutions vs. slow reference":
    check isEqualFrSeq( fieldConvolution(xs,ys) , naiveFieldConvolution(xs,ys) )

  test "FFT group convolutions vs. slow reference":
    check isEqualG1Seq( groupConvolution(xs,gs) , naiveGroupConvolution(xs,gs) )

#-------------------------------------------------------------------------------
