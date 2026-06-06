
#
# "a partial prover": precalculate stuff based on a partial witness
#
# this is useful when a significant portion of the witness does not changes 
# between proofs - then such precalculation can result in a potentially big speedup.
#
# an example use is RLN proofs, where the circuit is dominated by the Merkle 
#vinclusion proof check, which is expected to change relatively rarely.
#

{.push raises:[].}

import std/times
import std/options
import system

import taskpools

import groth16/bn128
import groth16/zkey_types
import groth16/misc

import groth16/partial/types

#import groth16/math/domain
#import groth16/math/poly
#import groth16/prover/shared

#-------------------------------------------------------------------------------
# the prover
#

proc generatePartialProof*( zkey: ZKey, pwtns: PartialWitness, pool: Taskpool, printTimings: bool): PartialProof =

  # assert( zkey.header.curve == wtns.curve )

  let partial_witness = pwtns.values

  let hdr  : GrothHeader  = zkey.header
  let spec : SpecPoints   = zkey.specPoints
  let pts  : ProverPoints = zkey.pPoints     

  let nvars = hdr.nvars
  let npubs = hdr.npubs

  assert( nvars == partial_witness.len , "wrong witness length" )

  assert( partial_witness.len == pts.pointsA1.len )
  assert( partial_witness.len == pts.pointsB1.len )
  assert( partial_witness.len == pts.pointsB2.len )

  # note: in "zs" we have to ignore public inputs (plus 1 for the special first entry)
  var partial_mask : seq[bool] = newSeq[bool](nvars)
  for i in 0..<nvars:
    partial_mask[i] = false

  let startIdx: int = npubs + 1  # in "zs" we have to skip the public IO
  var count:    int = 0          # count the number of existing witness elements 
  var zs_cnt:   int = 0          # count the number of existing witness elements, which are NOT public IO
  for i in 0..<startidx:
    if isSome(partial_witness[i]):
      partial_mask[i] = true
      count  += 1
  for i in startIdx..<nvars:
    if isSome(partial_witness[i]):
      partial_mask[i] = true
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
      if isSome(partial_witness[i]):
        let wtns_value = partial_witness[i].unsafeGet()
        compact_witness[j]  = wtns_value
        compact_pointsA1[j] = pts.pointsA1[ i ]
        compact_pointsB1[j] = pts.pointsB1[ i ]
        compact_pointsB2[j] = pts.pointsB2[ i ]
        j += 1
        if i >= startIdx:
          compact_zs[k]       = wtns_value
          compact_pointsC1[k] = pts.pointsC1[ i - startIdx ]
          k += 1
      
#[
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

]#

  assert( hdr.domainSize    == pts.pointsH1.len )
  assert( nvars - npubs - 1 == pts.pointsC1.len )

  var pi_a : G1 
  withMeasureTime(printTimings,"computing partial pi_A (G1 MSM)"):
    pi_a =  spec.alpha1
    pi_a += msmMultiThreadedG1( compact_witness , compact_pointsA1 , pool )

  var rho : G1 
  withMeasureTime(printTimings,"computing partial rho (G1 MSM)"):
    rho =  spec.beta1
    rho += msmMultiThreadedG1(  compact_witness , compact_pointsB1 , pool )

  var pi_b : G2
  withMeasureTime(printTimings,"computing partial pi_B (G2 MSM)"):
    pi_b =  spec.beta2
    pi_b += msmMultiThreadedG2( compact_witness , compact_pointsB2 , pool )

  var pi_c : G1
  withMeasureTime(printTimings,"partial computing pi_C (G1 MSM)"):
    pi_c =  msmMultiThreadedG1( compact_zs , compact_pointsC1 , pool )

  return PartialProof(
      partial_mask  : partial_mask
    , partial_pi_a  : pi_a
    , partial_rho   : rho
    , partial_pi_b  : pi_b
    , partial_pi_c  : pi_c
    )

#-------------------------------------------------------------------------------

