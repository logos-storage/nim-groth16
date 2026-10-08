
pragma circom 2.1.1;

include "parametric_bench.circom";

component main { public [staticInp,updateInp] } = Benchmark(1234,567);