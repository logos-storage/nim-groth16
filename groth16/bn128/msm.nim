
#
# Multi-Scalar Multiplication (MSM)
# 

import system
import taskpools

import constantine/platforms/abstractions       except Subgroup
# import constantine/math/endomorphisms/frobenius except Subgroup

import constantine/math/io/io_bigints
import constantine/named/properties_fields except Subgroup
import constantine/math/arithmetic
import constantine/math/io/io_fields

import constantine/math/extension_fields/towers                 as ext
import constantine/math/elliptic/ec_shortweierstrass_affine     as aff except Subgroup
import constantine/math/elliptic/ec_shortweierstrass_projective as prj except Subgroup
import constantine/math/elliptic/ec_scalar_mul_vartime          as scl except Subgroup
import constantine/math/elliptic/ec_multi_scalar_mul            as msm except Subgroup

#import groth16/bn128/fields
import groth16/bn128/curves as mycurves
import groth16/sharedbuf

#import groth16/misc    # TEMP DEBUGGING
#import std/cpuinfo
#import std/times

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
# spawnable wrappers: take SharedBuf views, delegate to the core.
# These are what `pool.spawn` calls, so they carry {.gcsafe, raises: [].}.
#
# Local aliases `AffG1`/`AffG2` are required because taskpools' `spawn` macro
# does `getImpl().replaceSymsByIdents()`, which strips qualifications. With
# `SharedBuf[mycurves.G1]` the bare ident `G1` then re-resolves to the enum
# value `aff.G1` (of type `Subgroup`), not the type alias. Renaming dodges
# the collision.

type
  AffG1 = mycurves.G1
  AffG2 = mycurves.G2
  FrBN  = Fr[BN254_Snarks]

func msmConstantineG1Range( coeffs: SharedBuf[FrBN] ,
                            points: SharedBuf[AffG1] ): AffG1
                            {.gcsafe, raises: [].} =
  msmConstantineG1( toOpenArray(coeffs.payload, 0, coeffs.len - 1),
                    toOpenArray(points.payload, 0, points.len - 1) )

func msmConstantineG2Range( coeffs: SharedBuf[FrBN] ,
                            points: SharedBuf[AffG2] ): AffG2
                            {.gcsafe, raises: [].} =
  msmConstantineG2( toOpenArray(coeffs.payload, 0, coeffs.len - 1),
                    toOpenArray(points.payload, 0, points.len - 1) )

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
    let cs = SharedBuf.view(toOpenArray(coeffs, a, b - 1))
    let ps = SharedBuf.view(toOpenArray(points, a, b - 1))
    pending[k] = pool.spawn msmConstantineG1Range(cs, ps)
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
    let cs = SharedBuf.view(toOpenArray(coeffs, a, b - 1))
    let ps = SharedBuf.view(toOpenArray(points, a, b - 1))
    pending[k] = pool.spawn msmConstantineG2Range(cs, ps)
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

