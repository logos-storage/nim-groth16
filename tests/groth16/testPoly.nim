
{.used.} 

import std/unittest
import std/sequtils
import std/sugar

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128/fields
import groth16/bn128/arrays
# import groth16/bn128/rnd

import groth16/math/domain
import groth16/math/poly
# import groth16/math/ntt

#-------------------------------------------------------------------------------

type F = Fr[BN254_Snarks]

suite "polynomial tests":

  test "sanity check for (division by) vanishing polynomial":

    var js : seq[int] = toSeq(101..112)
    let cs : seq[F]   = map( js, intToFr )
    let P  : Poly     = Poly( coeffs:cs )

    let ok1 = polyDegree(P) == 11

    let n  : int = 5
    let QR = polyQuotRemByVanishing(P, n)
    let Q  = QR.quot
    let R  = QR.rem

    let ok2 = (polyDegree(R) < 5) and (polyDegree(Q) == 11-5)

    let Z : Poly = vanishingPoly(n)
    let S : Poly = Q * Z + R

    let ok3 = polyIsEqual(P,S) 

    check (ok1 and ok2 and ok3)

  #-----------------------------------------------------------------------------

  test "sanity check for NTT":

    var js : seq[int] = toSeq(101..116)
    let cs : seq[F]   = map( js, intToFr )
    let P  : Poly     = Poly( coeffs:cs )
    let D  : Domain   = createDomain(16)
    let xs : seq[F]   = D.enumerateDomain()
    let ys : seq[F]   = collect( newSeq, (for x in xs: polyEvalAt(P,x)) ) 
    let zs : seq[F]   = polyForwardNTT(P ,D)
    let Q  : Poly     = polyInverseNTT(zs,D)

    let ok1 = isEqualFrSeq( ys , zs )         # naive evaluations = NTT
    let ok2 = isEqualFrSeq( Q.coeffs , cs)    # inverse NTT = starting

    check (ok1 and ok2)

  #-----------------------------------------------------------------------------

  test "sanity check for FFT polynomial multiplication":

    var js : seq[int] = toSeq(101..110)
    let cs : seq[F]   = map( js, intToFr )
    let P  : Poly     = Poly( coeffs:cs )
  
    var ks : seq[int] = toSeq(1001..1020)
    let ds : seq[F]   = map( ks, intToFr )
    let Q  : Poly     = Poly( coeffs:ds )
  
    let R1 : Poly = polyMulNaive( P , Q )
    let R2 : Poly = polyMulFFT(   P , Q )
  
    # debugPrintFrSeq("naive coeffs", R1.coeffs)
    # debugPrintFrSeq("fft coeffs",   R2.coeffs)
  
    check  polyIsEqual(R1,R2) 

  #-----------------------------------------------------------------------------

  test "sanity check for Lagrange bases":

    let n  = 8
    let D  = createDomain(n)
  
    let L : seq[Poly] = collect( newSeq, (for k in 0..<n: lagrangePoly(D,k) ))
  
    let xs = enumerateDomain(D)
    let ys0 : seq[F]  = collect( newSeq, (for x in xs: polyEvalAt(L[0],x) )) 
    let ys1 : seq[F]  = collect( newSeq, (for x in xs: polyEvalAt(L[1],x) )) 
    let ys5 : seq[F]  = collect( newSeq, (for x in xs: polyEvalAt(L[5],x) )) 
    let zs0 : seq[F]  = collect( newSeq, (for i in 0..<n: deltaFr(0,i)  )) 
    let zs1 : seq[F]  = collect( newSeq, (for i in 0..<n: deltaFr(1,i)  )) 
    let zs5 : seq[F]  = collect( newSeq, (for i in 0..<n: deltaFr(5,i)  )) 
     
    let zeta = intToFr(123457)
    let us : seq[F]  = collect( newSeq, (for i in 0..<n: polyEvalAt(L[i],zeta)) ) 
    let vs : seq[F]  = collect( newSeq, (for i in 0..<n: evalLagrangePolyAt(D,i,zeta)) ) 
  
    check ( ys0===zs0 and ys1===zs1 and ys5===zs5 and us===vs )

#-------------------------------------------------------------------------------
