
Groth16 prover written in Nim
-----------------------------

This is a Groth16 prover implementation in Nim, using the 
[`constantine`](https://github.com/mratsim/constantine)
library as an arithmetic / curve backend.

The implementation is compatible with the `circom` + `snarkjs` ecosystem.

At the moment only the `BN254` (aka. `alt-bn128`) curve is supported.

### Partial proofs

Groth16 has an interesting internal structure in that most of the proof
computation is _linear_ in the witness.

This implies that if you have to generate a lot of proofs where a large part
of the witness is unchanged, then you can precalculate quite a lot, making
subsequent proofs faster. A situation where this happens is when you have
a long-term identity, which however you don't want to reveal (eg. anonymous 
credentials).

As splitting the proof into two has negligible overhead in practice, this
is an essentially free win in such situations. However the speedup potential
is limited, because the single nonlinear term depends on the full witness.
In practice we've seen 2x-3x speedups when 10% of the witness changes.

This is implemented under `groth16/partial/`.

### Semi-dynamic proofs

In a recent paper ["Dynark: Making Groth16 Dynamic"](https://eprint.iacr.org/2025/1897) 
(2025), a solution was proposed for also speeding up the computation of the nonlinear 
term. This results in an asymptotic speedup, however, unlike the previous 
(essentially trivial) solution, this one is not _for free_: There is a somewhat 
heavy precalculation when the full witness changes.

Whether that's worth it depends on the particular use-case.

A version of this idea is implemented under `groth16/dynamic/`.

### License

Licensed and distributed under either of the
[MIT license](http://opensource.org/licenses/MIT) or
[Apache License, v2.0](http://www.apache.org/licenses/LICENSE-2.0),
at your choice. 

### TODO

- [ ] clean up the code
- [ ] more comprehensive testing
- [ ] compare `.r1cs` to the "coeffs" section of `.zkey`
- [x] generate fake circuit-specific setup ourselves
- [x] make a CLI interface
- [x] multithreading support (MSM, and possibly also FFT)
- [ ] add Groth16 explanatory notes
- [ ] document the `snarkjs` circuit-specific setup `H` points convention
- [x] precalculate stuff for "partial" proofs
- [ ] implement Dynark style "dynamic proofs" too
- [ ] benchmarks
- [ ] make it work for different curves
