
{.used.} 

import std/unittest

# import constantine/math/arithmetic
# import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/rnd
import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/ntt
import groth16/math/poly
import groth16/math/convolution

import groth16/dynamic/types
import groth16/dynamic/shared
import groth16/dynamic/v1/preprocess
#import groth16/dynamic/v1/setup
#import groth16/dynamic/v1/finish

#-------------------------------------------------------------------------------

func ones(N : int): seq[F] = constSeq(N, oneFr)

#-------------------------------------------------------------------------------

suite "dynamic proof (Dynark) tests":

  let N : int = 128
  let K : int =  32
  let D  = createDomain( N )
  let sg = createSubgroup( D , K )
  let tau   : F = randFr()
  let delta : F = randFr()

  test "Wvec vs Wcoeff":
    let xs = calculateWVec(D) 
    var ys: seq[F] = newSeq[F]( N )
    for i in 0..<N:
      ys[i] = Wcoeff( D , i )
    check isEqualFrSeq( xs , ys )

  test "WvecBar vs Wcoeff":
    let xs = calculateWVecBar(D) 
    var ys: seq[F] = newSeq[F]( N )
    for i in 0..<N:
      ys[i] = Wcoeff( D , -i )
    check isEqualFrSeq( xs , ys )

  test "sum of Wvec":
    check sumOfWVec(N) === sumSeqFr( calculateWVec(D) )

  test "FFT of Wvec":
    var xs = ones(N)
    inplaceMulByFFTofWVec( xs ) 
    let ys = forwardNTT( calculateWVec(D) , D ) 
    check isEqualFrSeq( xs , ys )

  test "FFT of WvecBar":
    var xs = ones(N)
    inplaceMulByFFTofWVecBar( xs ) 
    let ys = forwardNTT( calculateWVecBar(D) , D ) 
    check isEqualFrSeq( xs , ys)

  test "field convolution with Wvec":
    let xs  = randFrSeq( N )
    let lhs = fieldConvolveWithWVec( D , xs ) 
    let rhs = fieldConvolution( xs , calculateWVec(D) )
    check isEqualFrSeq( lhs , rhs )  

  test "field convolution with WvecBar":
    let xs  = randFrSeq( N )
    let lhs = fieldConvolveWithWVecBar( D , xs ) 
    let rhs = fieldConvolution( xs , calculateWVecBar(D) )
    check isEqualFrSeq( lhs , rhs )  

  test "field convolution with WvecBar on a subgroup":
    let xs  = randFrSeq( N )
    let lhs = fieldConvolveWithWVecBarOnSubgroup( sg , xs ) 
    let rhs = selectOnSubgroup( sg , fieldConvolution( xs , calculateWVecBar(D) ) )
    check isEqualFrSeq( lhs , rhs )  

  test "group convolution with Wvec":
    let gs  = randG1Seq( N )
    let lhs = groupConvolveWithWVec( D , gs ) 
    let rhs = groupConvolution( calculateWVec(D) , gs )
    check isEqualG1Seq( lhs , rhs )  

  test "group convolution with WvecBar":
    let gs  = randG1Seq( N )
    let lhs = groupConvolveWithWVecBar( D , gs ) 
    let rhs = groupConvolution( calculateWVecBar(D) , gs )
    check isEqualG1Seq( lhs , rhs )  

  test "group convolution with Wvec on a subgroup":
    let gs  = randG1Seq( N )
    let lhs = groupConvolveWithWVecOnSubgroup( sg , gs ) 
    let rhs = selectOnSubgroup( sg , groupConvolution( calculateWVec(D) , gs ) )
    check isEqualG1Seq( lhs , rhs )  

  test "expansion of the product of Lagrange polynomials":
    var ok = true
    for i in 0..<N:
      for k in 0..<N:
        let lhs = evalLagrangeProductViaLemma(D, i, k, tau)
        let rhs = evalLagrangePolyAt(D, i, tau) * evalLagrangePolyAt(D, k, tau)
        let local_ok = (lhs === rhs) 
        ok = ok and local_ok
    check ok

  test "cross-term coefficients":
    check testCrossTermCoeffs(N)

  test "projection term V1 calculation":
    check testProjectionElementsV1( N , tau , delta )

#-------------------------------------------------------------------------------
