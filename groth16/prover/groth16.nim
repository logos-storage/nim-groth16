
#
# Groth16 prover
#
# WARNING! 
# the points H in `.zkey` are *NOT* what normal people would think they are
# See <https://geometry.xyz/notebook/the-hidden-little-secret-in-snarkjs>
#

{.push raises:[].}

import std/times
import system
import taskpools
import constantine/math/arithmetic
import constantine/named/properties_fields

#import constantine/math/io/io_extfields except Fp12

import groth16/bn128
#import groth16/math/domain
import groth16/math/poly
import groth16/zkey_types
import groth16/files/witness
import groth16/misc

import groth16/prover/types
import groth16/prover/shared

#-------------------------------------------------------------------------------
# the prover
#

proc generateProofWithMask*( zkey: ZKey, wtns: Witness, mask: Mask, pool: Taskpool, printTimings: bool): Proof {.raises: [ValueError].} =

  if zkey.header.curve != wtns.curve:
    raise newException(ValueError,
      "curve mismatch: zkey.header.curve=" & $zkey.header.curve &
      " wtns.curve=" & $wtns.curve)

  let witness = wtns.values

  let hdr  : GrothHeader  = zkey.header
  let spec : SpecPoints   = zkey.specPoints
  let pts  : ProverPoints = zkey.pPoints

  let nvars = hdr.nvars
  let npubs = hdr.npubs

  if nvars != witness.len:
    raise newException(ValueError,
      "wrong witness length: expected " & $nvars &
      " (zkey.header.nvars), got " & $witness.len)

  # remark: with the special variable "1" we actuall have (npub+1) public IO variables
  var pubIO = newSeq[Fr[BN254_Snarks]]( npubs + 1)
  for i in 0..npubs: pubIO[i] = witness[i]             

  var abc : ABC 
  withMeasureTime(printTimings,"building 'ABC'"):
    abc = buildABC( zkey, witness )

  var qs : seq[Fr[BN254_Snarks]]
  withMeasureTime(printTimings,"computing the quotient (FFTs)"):
    case zkey.header.flavour

      # the points H are [delta^-1 * tau^i * Z(tau)]
      of JensGroth:
        let polyQ = computeQuotientPointwise( abc, pool )
        qs = polyQ.coeffs
  
      # the points H are `[delta^-1 * L_{2i+1}(tau)]_1`
      # where L_i are Lagrange basis polynomials on the double-sized domain
      of Snarkjs:
        qs = computeSnarkjsScalarCoeffs( abc, pool )

  var zs = newSeq[Fr[BN254_Snarks]]( nvars - npubs - 1 )
  for j in npubs+1..<nvars:
    zs[j-npubs-1] = witness[j]

  # masking coeffs
  let r = mask.r
  let s = mask.s

  assert( witness.len == pts.pointsA1.len )
  assert( witness.len == pts.pointsB1.len )
  assert( witness.len == pts.pointsB2.len )
  assert( hdr.domainSize    == qs.len           )
  assert( hdr.domainSize    == pts.pointsH1.len )
  assert( nvars - npubs - 1 == zs.len           )
  assert( nvars - npubs - 1 == pts.pointsC1.len )

  var pi_a : G1 
  withMeasureTime(printTimings,"computing pi_A (G1 MSM)"):
    pi_a =  spec.alpha1
    pi_a += r ** spec.delta1
    pi_a += msmMultiThreadedG1( witness , pts.pointsA1, pool )

  var rho : G1 
  withMeasureTime(printTimings,"computing rho (G1 MSM)"):
    rho =  spec.beta1
    rho += s ** spec.delta1
    rho += msmMultiThreadedG1(  witness , pts.pointsB1, pool )

  var pi_b : G2
  withMeasureTime(printTimings,"computing pi_B (G2 MSM)"):
    pi_b =  spec.beta2
    pi_b += s ** spec.delta2
    pi_b += msmMultiThreadedG2( witness , pts.pointsB2, pool )

  var pi_c : G1
  withMeasureTime(printTimings,"computing pi_C (2x G1 MSM)"):
    pi_c =  s ** pi_a
    pi_c += r ** rho
    pi_c += negFr(r*s) ** spec.delta1
    pi_c += msmMultiThreadedG1( qs , pts.pointsH1, pool )
    pi_c += msmMultiThreadedG1( zs , pts.pointsC1, pool )

  return Proof( curve:"bn128", publicIO:pubIO, pi_a:pi_a, pi_b:pi_b, pi_c:pi_c )

#-------------------------------------------------------------------------------

proc generateProofWithTrivialMask*( zkey: ZKey, wtns: Witness, pool: Taskpool, printTimings: bool ): Proof {.raises: [ValueError].} =
  let mask = Mask( r: zeroFr , s: zeroFr )
  return generateProofWithMask( zkey, wtns, mask, pool, printTimings )

proc generateProof*( zkey: ZKey, wtns: Witness, pool: Taskpool, printTimings = false ): Proof {.raises: [ValueError].} =
  let mask = randomMask()
  return generateProofWithMask( zkey, wtns, mask, pool, printTimings )

#-------------------------------------------------------------------------------
