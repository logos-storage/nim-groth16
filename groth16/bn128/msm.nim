
#
# Multi-Scalar Multiplication (MSM)
# 

import system
import std/cpuinfo
import taskpools

# import constantine/curves_primitives except Fp, Fp2, Fr
 
import constantine/platforms/abstractions   except Subgroup
import constantine/math/endomorphisms/frobenius except Subgroup

import constantine/math/io/io_bigints
import constantine/named/properties_fields except Subgroup
import constantine/math/arithmetic
import constantine/math/io/io_fields

import constantine/math/extension_fields/towers                 as ext
import constantine/math/elliptic/ec_shortweierstrass_affine     as aff except Subgroup
import constantine/math/elliptic/ec_shortweierstrass_projective as prj except Subgroup
import constantine/math/elliptic/ec_scalar_mul_vartime          as scl except Subgroup
import constantine/math/elliptic/ec_multi_scalar_mul            as msm except Subgroup

import groth16/bn128/fields
import groth16/bn128/curves as mycurves

import groth16/misc    # TEMP DEBUGGING
import std/times

#-------------------------------------------------------------------------------

proc msmConstantineG1*( coeffs: openArray[Fr[BN254_Snarks]] , points: openArray[G1] ): G1 =

  # let start = cpuTime()

  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )

  var bigcfs : seq[BigInt[254]]
  for x in coeffs:
    bigcfs.add( x.toBig() )

  var r : ProjG1

  # [Fp,aff.G1]
  msm.multiScalarMul_vartime( r,
    toOpenArray(bigcfs, 0, N-1),
    toOpenArray(points, 0, N-1) )

  var rAff: G1
  prj.affine(rAff, r)

  # let elapsed = cpuTime() - start
  # echo("computing an MSM of size " & ($N) & " took " & seconds(elapsed))

  return rAff

#---------------------------------------

func msmConstantineG2*( coeffs: openArray[Fr[BN254_Snarks]] , points: openArray[G2] ): G2 =

  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )

  var bigcfs : seq[BigInt[254]]
  for x in coeffs:
    bigcfs.add( x.toBig() )

  var r : ProjG2

  # note: at the moment of writing this, `multiScalarMul_vartime` is buggy.
  # however, the "reference" one is _much_ slower.
  msm.multiScalarMul_reference_vartime( r,
    toOpenArray(bigcfs, 0, N-1),
    toOpenArray(points, 0, N-1) )

  var rAff: G2
  prj.affine(rAff, r)

  return rAff

#-------------------------------------------------------------------------------

const task_multiplier : int = 1

proc msmMultiThreadedG1*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G1], pool: Taskpool ): G1 =

  # for N <= 255 , we use 1 thread
  # for N == 256 , we use 2 threads
  # for N == 512 , we use 4 threads 
  # for N >= 1024, we use 8+ threads 

  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )
  let nthreads_target = min( pool.numThreads, 256 )
  let nthreads = max( 1 , min( N div 128 , nthreads_target ) )
  let ntasks   = if nthreads>1: (nthreads*task_multiplier) else: 1

  var pending : seq[FlowVar[mycurves.G1]] = newSeq[FlowVar[mycurves.G1]](ntasks)

  var a : int = 0
  var b : int
  for k in 0..<ntasks:
    if k < ntasks-1:
      b = (N*(k+1)) div ntasks
    else:
      b = N
    let cs = coeffs[a..<b]
    let ps = points[a..<b]
    pending[k] = pool.spawn msmConstantineG1( cs, ps );
    a = b

  var res : G1 = infG1
  for k in 0..<ntasks:
    res += sync pending[k]

  return res

#---------------------------------------

proc msmMultiThreadedG2*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G2], pool: Taskpool ): G2 =

  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )
  let nthreads_target = min( pool.numThreads, 256 )
  let nthreads = max( 1 , min( N div 128 , nthreads_target ) )
  let ntasks   = if nthreads>1: (nthreads*task_multiplier) else: 1

  var pending : seq[FlowVar[mycurves.G2]] = newSeq[FlowVar[mycurves.G2]](ntasks)

  var a : int = 0
  var b : int
  for k in 0..<ntasks:
    if k < ntasks-1:
      b = (N*(k+1)) div ntasks
    else:
      b = N
    let cs = coeffs[a..<b]
    let ps = points[a..<b]
    pending[k] = pool.spawn msmConstantineG2( cs, ps );
    a = b

  var res : G2 = infG2
  for k in 0..<ntasks:
    res += sync pending[k]

  return res

#-------------------------------------------------------------------------------

func msmNaiveG1*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G1] ): G1 =
  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )

  var s : ProjG1
  s.setNeutral()

  for i in 0..<N:
    var t : ProjG1
    prj.fromAffine( t, points[i] )
    scl.scalarMul_vartime( t , coeffs[i].toBig() )
    s += t

  var r : G1
  prj.affine( r, s )

  return r

#---------------------------------------

func msmNaiveG2*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G2] ): G2 =
  let N = coeffs.len
  assert( N == points.len, "incompatible sequence lengths" )

  var s : ProjG2
  s.setNeutral()

  for i in 0..<N:
    var t : ProjG2
    prj.fromAffine( t, points[i] )
    scl.scalarMul_vartime( t , coeffs[i].toBig() )
    s += t

  var r : G2
  prj.affine( r, s)

  return r

#-------------------------------------------------------------------------------

proc msmG1*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G1] ): G1 = msmConstantineG1(coeffs, points)
proc msmG2*( coeffs: seq[Fr[BN254_Snarks]] , points: seq[G2] ): G2 = msmConstantineG2(coeffs, points)

#-------------------------------------------------------------------------------

