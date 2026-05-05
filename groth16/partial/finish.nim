
#
# Finish a "partial proof"
#

{.push raises:[].}

import std/times
import system
import taskpools
import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
#import groth16/math/domain
import groth16/math/poly
import groth16/zkey_types
import groth16/files/witness
import groth16/misc

import groth16/partial/types
import groth16/prover/types
import groth16/prover/shared

#-------------------------------------------------------------------------------
# the finishing prover
#

proc finishPartialProofWithMask*( zkey: ZKey, wtns: Witness, partialProof: PartialProof, mask: Mask, pool: Taskpool, printTimings: bool): Proof =

  # if (zkey.header.curve != wtns.curve):
  #   echo( "zkey.header.curve = " & ($zkey.header.curve) )
  #   echo( "wtns.curve        = " & ($wtns.curve       ) )

  assert( zkey.header.curve == wtns.curve )

  let witness = wtns.values

  let hdr  : GrothHeader  = zkey.header
  let spec : SpecPoints   = zkey.specPoints
  let pts  : ProverPoints = zkey.pPoints     

  let nvars = hdr.nvars
  let npubs = hdr.npubs

  assert( nvars == witness.len , "wrong witness length" )

  let partial_mask = partialProof.partial_mask

  # remark: with the special variable "1" we actuall have (npub+1) public IO variables
  var pubIO = newSeq[Fr[BN254_Snarks]]( npubs + 1)
  for i in 0..npubs: pubIO[i] = witness[i]             

  # note: we have to ignore public inputs (plus 1 for the special first entry)
  let startIdx : int = npubs + 1
  var count:  int = 0       # count the number of NEW witness elements 
  var zs_cnt: int = 0       # count the number of NEW witness elements, which are NOT public IO
  for i in 0..<startidx:
    if (not partial_mask[i]):
      count  += 1
  for i in startIdx..<nvars:
    if (not partial_mask[i]):
      count  += 1
      zs_cnt += 1

  # compactify the stuff so we can call normal vector MSM
  var compact_witness:  seq[F]  = newSeq[F](count)
  var compact_pointsA1: seq[G1] = newSeq[G1](count)
  var compact_pointsB1: seq[G1] = newSeq[G1](count)
  var compact_pointsB2: seq[G2] = newSeq[G2](count)
  var compact_zs:       seq[F]  = newSeq[F](zs_cnt)
  var compact_pointsC1: seq[G1] = newSeq[G1](zs_cnt)
  block:
    var j: int = 0;
    var k: int = 0;
    for i in 0..<nvars:
      if (not partial_mask[i]):
        compact_witness[j ] = witness[i]
        compact_pointsA1[j] = pts.pointsA1[i]
        compact_pointsB1[j] = pts.pointsB1[i]
        compact_pointsB2[j] = pts.pointsB2[i]
        j += 1
        if i >= startIdx:
          compact_zs[k]       = witness[i]
          compact_pointsC1[k] = pts.pointsC1[ i - startIdx ]
          k += 1

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

  # masking coeffs
  let r = mask.r
  let s = mask.s

  assert( witness.len == pts.pointsA1.len )
  assert( witness.len == pts.pointsB1.len )
  assert( witness.len == pts.pointsB2.len )
  assert( hdr.domainSize    == qs.len           )
  assert( hdr.domainSize    == pts.pointsH1.len )
  assert( nvars - npubs - 1 == pts.pointsC1.len )

  var pi_a : G1 = partialProof.partial_pi_a
  withMeasureTime(printTimings,"computing pi_A (G1 MSM)"):
    pi_a += r ** spec.delta1
    pi_a += msmMultiThreadedG1( compact_witness , compact_pointsA1, pool )

  var rho : G1 = partialProof.partial_rho
  withMeasureTime(printTimings,"computing rho (G1 MSM)"):
    rho += s ** spec.delta1
    rho += msmMultiThreadedG1( compact_witness , compact_pointsB1, pool )

  var pi_b : G2 = partialProof.partial_pi_b
  withMeasureTime(printTimings,"computing pi_B (G2 MSM)"):
    pi_b += s ** spec.delta2
    pi_b += msmMultiThreadedG2( compact_witness , compact_pointsB2, pool )

  var pi_c : G1 = partialProof.partial_pi_c
  withMeasureTime(printTimings,"computing pi_C (2x G1 MSM)"):
    pi_c += s ** pi_a
    pi_c += r ** rho
    pi_c += negFr(r*s) ** spec.delta1
    pi_c += msmMultiThreadedG1( qs , pts.pointsH1, pool )
    pi_c += msmMultiThreadedG1( compact_zs , compact_pointsC1, pool )

  return Proof( curve:"bn128", publicIO:pubIO, pi_a:pi_a, pi_b:pi_b, pi_c:pi_c )

#-------------------------------------------------------------------------------

proc finishPartialProofWithTrivialMask*( zkey: ZKey, wtns: Witness, partial: PartialProof, pool: Taskpool, printTimings: bool ): Proof =
  let mask = Mask( r: zeroFr , s: zeroFr )
  return finishPartialProofWithMask( zkey, wtns, partial, mask, pool, printTimings )

proc finishPartialProof*( zkey: ZKey, wtns: Witness, partial: PartialProof, pool: Taskpool, printTimings = false ): Proof =
  let mask = randomMask()
  return finishPartialProofWithMask( zkey, wtns, partial, mask, pool, printTimings )

#-------------------------------------------------------------------------------
