
#
# nimble build -d:release
#

import strformat

import taskpools

import constantine/math/arithmetic
import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/msm
import groth16/bn128/arrays

import groth16/math/domain
import groth16/math/ntt
import groth16/dynamic/group_fft

import shared

#-------------------------------------------------------------------------------

when isMainModule:

  echo "quick & dirty FFT benchmarks"
  
  let nthreads: int = 8
  let N:        int = 8192

  let D = createDomain(N)
  let xs : seq[Fr[BN254_Snarks]] = randFrSeq(N)
  let gs : seq[AffG1]            = randG1Seq(N)

  var ys, zs: seq[Fr[BN254_Snarks]]
  var hs, rs: seq[AffG1]

  var pool = Taskpool.new(nthreads)

  #-----------------------------------------------------------------------------

  echo "\n----------------------------------------"
  echo "*** scalar multiplications\n"

  hs = newSeq[AffG1](N)
  withMeasureTime(true , fmt"{N} individual scalar multiplications      "):
    for i in 0..<N: 
      hs[i] = xs[i] ** gs[i]

  var sum_naive: AffG1
  withMeasureTime(true , fmt"naive simulated MSM of size {N}            "):
    for i in 0..<N: 
      hs[i] = xs[i] ** gs[i]
      if i==0:
        sum_naive = hs[i]
      else:
        sum_naive += hs[i]

  var sum_msm: AffG1
  withMeasureTime(true , fmt"proper MSM of size {N}                     "):
    sum_msm = msmConstantineG1( xs, gs )

  var sum_multi: AffG1
  withMeasureTime(true , fmt"multithreaded MSM of size {N} ({nthreads} threads)  "):
    let sum_multi = msmMultiThreadedG1( xs , gs , pool )

  echo "naive == msm : " & $(sum_naive === sum_msm)

  #-----------------------------------------------------------------------------

  echo "\n----------------------------------------"
  echo "*** field FFTs\n"

  withMeasureTime(true , fmt"field NTT of size {N}           "): 
    ys = forwardNTT(xs, D)

  withMeasureTime(true , fmt"unscaled field INTT of size {N} "): 
    zs = unscaledInverseNTT(ys, D) 

  withMeasureTime(true , fmt"field INTT of size {N}          "): 
    zs = inverseNTT(ys, D) 

  echo "xs == zs : " & $isEqualFrSeq(xs , zs)
 
  #-----------------------------------------------------------------------------

  echo "\n----------------------------------------"
  echo "*** group FFTs\n"

  withMeasureTime(true , fmt"group FFT of size {N}           "): 
    hs = forwardGroupFFT(gs, D)

  withMeasureTime(true , fmt"unscaled group IFFT of size {N} "): 
    rs = unscaledInverseGroupFFT(hs, D) 

  withMeasureTime(true , fmt"group IFFT of size {N}          "): 
    rs = inverseGroupFFT(hs, D) 

  echo "gs == rs : " & $isEqualG1Seq(gs , rs)

  echo "\n"

#-------------------------------------------------------------------------------
