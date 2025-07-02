
#
# Groth16 prover
#
# WARNING! 
# the points H in `.zkey` are *NOT* what normal people would think they are
# See <https://geometry.xyz/notebook/the-hidden-little-secret-in-snarkjs>
#

{.push raises:[].}

#[
import sugar
import constantine/math/config/curves 
import constantine/math/io/io_fields
import constantine/math/io/io_bigints
import ./zkey
]#

import std/os
import std/times
import std/cpuinfo
import system
import taskpools
import constantine/math/arithmetic
import constantine/named/properties_fields
import constantine/math/extension_fields/towers

#import constantine/math/io/io_extfields except Fp12

import groth16/bn128
import groth16/math/domain
import groth16/math/poly
import groth16/zkey_types
import groth16/files/witness
import groth16/misc

#-------------------------------------------------------------------------------

type
  Proof* = object
    publicIO* : seq[Fr[BN254_Snarks]]
    pi_a*     : G1
    pi_b*     : G2
    pi_c*     : G1
    curve*    : string

#-------------------------------------------------------------------------------
# Az, Bz, Cz column vectors
# 

type
  ABC = object
    valuesAz : seq[Fr[BN254_Snarks]]
    valuesBz : seq[Fr[BN254_Snarks]]
    valuesCz : seq[Fr[BN254_Snarks]]

# computes the vectors A*z, B*z, C*z where z is the witness
func buildABC( zkey: ZKey, witness: seq[Fr[BN254_Snarks]] ): ABC =
  let hdr: GrothHeader = zkey.header
  let domSize = hdr.domainSize

  var valuesAz = newSeq[Fr[BN254_Snarks]](domSize)
  var valuesBz = newSeq[Fr[BN254_Snarks]](domSize)

  for entry in zkey.coeffs:
    case entry.matrix 
      of MatrixA: valuesAz[entry.row] += entry.coeff * witness[entry.col]
      of MatrixB: valuesBz[entry.row] += entry.coeff * witness[entry.col]
      else: raise newException(AssertionDefect, "fatal error")

  var valuesCz = newSeq[Fr[BN254_Snarks]](domSize)
  for i in 0..<domSize:
    valuesCz[i] = valuesAz[i] * valuesBz[i]

  return ABC( valuesAz:valuesAz, valuesBz:valuesBz, valuesCz:valuesCz )

#-------------------------------------------------------------------------------
# quotient poly
#

# interpolates A,B,C, and computes the quotient polynomial Q = (A*B - C) / Z
func computeQuotientNaive( abc: ABC ): Poly=
  let n = abc.valuesAz.len
  assert( abc.valuesBz.len == n )
  assert( abc.valuesCz.len == n )
  let D = createDomain(n)
  let polyA : Poly = polyInverseNTT( abc.valuesAz , D )
  let polyB : Poly = polyInverseNTT( abc.valuesBz , D )
  let polyC : Poly = polyInverseNTT( abc.valuesCz , D )
  let polyBig = polyMulFFT( polyA , polyB ) - polyC
  var polyQ   = polyDivideByVanishing(polyBig, D.domainSize)
  polyQ.coeffs.add( zeroFr )    # make it a power of two
  return polyQ

#---------------------------------------

# returns [ eta^i * xs[i] | i<-[0..n-1] ]
func multiplyByPowers( xs: seq[Fr[BN254_Snarks]], eta: Fr[BN254_Snarks] ): seq[Fr[BN254_Snarks]] =
  let n = xs.len
  assert(n >= 1)
  var ys = newSeq[Fr[BN254_Snarks]](n)
  ys[0] = xs[0]
  if n >= 1: ys[1] = eta * xs[1]
  var spow = eta
  for i in 2..<n: 
    spow *= eta
    ys[i] = spow * xs[i]
  return ys

# interpolates a polynomial, shift the variable by `eta`, and compute the shifted values
func shiftEvalDomain(
  values: seq[Fr[BN254_Snarks]],
  D: Domain,
  eta: Fr[BN254_Snarks] ): seq[Fr[BN254_Snarks]] =
  let poly : Poly = polyInverseNTT( values , D )
  let cs : seq[Fr[BN254_Snarks]] = poly.coeffs
  var ds : seq[Fr[BN254_Snarks]] = multiplyByPowers( cs, eta )
  return polyForwardNTT( Poly(coeffs:ds), D )

# Wraps shiftEvalDomain such that it can be called by Taskpool.spawn. The result
# is written to the output parameter. Has an unused return type because
# Taskpool.spawn cannot handle a void return type.
func shiftEvalDomainTask(
  values: seq[Fr[BN254_Snarks]],
  D: Domain,
  eta: Fr[BN254_Snarks],
  output: ptr Isolated[seq[Fr[BN254_Snarks]]]): bool =

  output[] = isolate shiftEvalDomain(values, D, eta)

# computes the quotient polynomial Q = (A*B - C) / Z
# by computing the values on a shifted domain, and interpolating the result
# remark: Q has degree `n-2`, so it's enough to use a domain of size n
proc computeQuotientPointwise( abc: ABC, pool: TaskPool ): Poly =
  let n    = abc.valuesAz.len
  assert( abc.valuesBz.len == n )
  assert( abc.valuesCz.len == n )

  let D    = createDomain(n)
  
  # (eta*omega^j)^n - 1 = eta^n - 1 
  # 1 / [ (eta*omega^j)^n - 1] = 1/(eta^n - 1)
  let eta   = createDomain(2*n).domainGen
  let invZ1 = invFr( smallPowFr(eta,n) - oneFr )

  var outputA1, outputB1, outputC1: Isolated[seq[Fr[BN254_Snarks]]]

  var taskA1 = pool.spawn shiftEvalDomainTask( abc.valuesAz, D, eta, addr outputA1 )
  var taskB1 = pool.spawn shiftEvalDomainTask( abc.valuesBz, D, eta, addr outputB1 )
  var taskC1 = pool.spawn shiftEvalDomainTask( abc.valuesCz, D, eta, addr outputC1 )

  discard sync taskA1
  discard sync taskB1
  discard sync taskC1

  let A1 = outputA1.extract()
  let B1 = outputB1.extract()
  let C1 = outputC1.extract()

  var ys : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( n )
  for j in 0..<n: ys[j] = ( A1[j]*B1[j] - C1[j] ) * invZ1
  let Q1 = polyInverseNTT( ys, D )
  let cs = multiplyByPowers( Q1.coeffs, invFr(eta) )

  return Poly(coeffs: cs)

#---------------------------------------

# Snarkjs does something different, not actually computing the quotient poly
# they can get away with this, because during the trusted setup, they
# replace the points encoding the values `delta^-1 * tau^i * Z(tau)` by 
# (shifted) Lagrange bases.
# see <https://geometry.xyz/notebook/the-hidden-little-secret-in-snarkjs>
#
proc computeSnarkjsScalarCoeffs( abc: ABC, pool: TaskPool ): seq[Fr[BN254_Snarks]] =
  let n    = abc.valuesAz.len
  assert( abc.valuesBz.len == n )
  assert( abc.valuesCz.len == n )
  let D    = createDomain(n)
  let eta  = createDomain(2*n).domainGen

  var outputA1, outputB1, outputC1: Isolated[seq[Fr[BN254_Snarks]]]

  var taskA1 = pool.spawn shiftEvalDomainTask( abc.valuesAz, D, eta, addr outputA1 )
  var taskB1 = pool.spawn shiftEvalDomainTask( abc.valuesBz, D, eta, addr outputB1 )
  var taskC1 = pool.spawn shiftEvalDomainTask( abc.valuesCz, D, eta, addr outputC1 )

  discard sync taskA1
  discard sync taskB1
  discard sync taskC1

  let A1 = outputA1.extract()
  let B1 = outputB1.extract()
  let C1 = outputC1.extract()

  var ys : seq[Fr[BN254_Snarks]] = newSeq[Fr[BN254_Snarks]]( n )
  for j in 0..<n: ys[j] = ( A1[j] * B1[j] - C1[j] ) 

  return ys

#[

proc computeSnarkjsScalarCoeffs_st( abc: ABC ): seq[Fr] =
  let n    = abc.valuesAz.len
  assert( abc.valuesBz.len == n )
  assert( abc.valuesCz.len == n )
  let D    = createDomain(n)
  let eta  = createDomain(2*n).domainGen
  let A1 : seq[Fr] = shiftEvalDomain( abc.valuesAz, D, eta )
  let B1 : seq[Fr] = shiftEvalDomain( abc.valuesBz, D, eta )
  let C1 : seq[Fr] = shiftEvalDomain( abc.valuesCz, D, eta )
  var ys : seq[Fr] = newSeq[Fr]( n )
  for j in 0..<n: ys[j] = ( A1[j] * B1[j] - C1[j] ) 
  return ys

proc computeSnarkjsScalarCoeffs( nthreads: int, abc: ABC ): seq[Fr] =
  if nthreads <= 1:
    computeSnarkjsScalarCoeffs_st( abc )
  else:
    computeSnarkjsScalarCoeffs_mt( nthreads, abc )

]#

#-------------------------------------------------------------------------------
# the prover
#

type
  Mask* = object
    r*: Fr[BN254_Snarks]              # masking coefficients
    s*: Fr[BN254_Snarks]              # for zero knowledge

proc generateProofWithMask*( zkey: ZKey, wtns: Witness, mask: Mask, pool: Taskpool, printTimings: bool): Proof =

  when not (defined(gcArc) or defined(gcOrc) or defined(gcAtomicArc)):
    {.fatal: "Compile with arc/orc!".}

  # if (zkey.header.curve != wtns.curve):
  #   echo( "zkey.header.curve = " & ($zkey.header.curve) )
  #   echo( "wtns.curve        = " & ($wtns.curve       ) )

  assert( zkey.header.curve == wtns.curve )
  var start : float = 0

  let witness = wtns.values

  let hdr  : GrothHeader  = zkey.header
  let spec : SpecPoints   = zkey.specPoints
  let pts  : ProverPoints = zkey.pPoints     

  let nvars = hdr.nvars
  let npubs = hdr.npubs

  assert( nvars == witness.len , "wrong witness length" )

  # remark: with the special variable "1" we actuall have (npub+1) public IO variables
  var pubIO = newSeq[Fr[BN254_Snarks]]( npubs + 1)
  for i in 0..npubs: pubIO[i] = witness[i]             

  start = cpuTime()
  var abc : ABC 
  withMeasureTime(printTimings,"building 'ABC'"):
    abc = buildABC( zkey, witness )

  start = cpuTime()
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

proc generateProofWithTrivialMask*( zkey: ZKey, wtns: Witness, pool: Taskpool, printTimings: bool ): Proof =
  let mask = Mask( r: zeroFr , s: zeroFr )
  return generateProofWithMask( zkey, wtns, mask, pool, printTimings )

proc generateProof*( zkey: ZKey, wtns: Witness, pool: Taskpool, printTimings = false ): Proof =

  # masking coeffs
  let r = randFr()
  let s = randFr()
  let mask = Mask(r: r, s: s)

  return generateProofWithMask( zkey, wtns, mask, pool, printTimings )
