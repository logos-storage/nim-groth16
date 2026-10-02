
version     = "0.0.1"
author      = "Balazs Komuves"
description = "FFT benchmarks"
license     = "MIT OR Apache-2.0"

binDir      = "build"
bin         = @["bench_fft"]

requires "nim >= 2.2.0"
requires "https://github.com/status-im/nim-taskpools >= 0.0.5"
requires "https://github.com/mratsim/constantine"
requires "groth16 >= 0.1.2"
