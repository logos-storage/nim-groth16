
Groth16 prover written in Nim
-----------------------------

This is Groth16 prover implementation in Nim, using the 
[`constantine`](https://github.com/mratsim/constantine)
library as an arithmetic / curve backend.

The implementation is compatible with the `circom` + `snarkjs` ecosystem.

At the moment only the `BN254` (aka. `alt-bn128`) curve is supported.

### License

Licensed and distributed under either of the
[MIT license](http://opensource.org/licenses/MIT) or
[Apache License, v2.0](http://www.apache.org/licenses/LICENSE-2.0),
at your choice. 

### TODO

- [x] find and fix the _second_ totally surreal bug
- [ ] clean up the code
- [x] make it compatible with the latest constantine and also Nim 2.0.x
- [x] make it a nimble package
- [ ] compare `.r1cs` to the "coeffs" section of `.zkey`
- [x] generate fake circuit-specific setup ourselves
- [x] make a CLI interface
- [x] multithreading support (MSM, and possibly also FFT)
- [ ] add Groth16 notes
- [ ] document the `snarkjs` circuit-specific setup `H` points convention
- [x] precalculate stuff for "partial" proofs
- [ ] implement Dynark style "dynamic proofs" too
- [ ] make it work for different curves
