
"Semi-dynamic" Groth16 proofs
-----------------------------

This is loosely based on the [Dynark paper](https://eprint.iacr.org/2025/1897):

- _"Dynark: Making Groth16 Dynamic"_ by Tianyu Zhang, Yupeng Ouyang and Yupeng Zhang

See also [this write-up](https://hackmd.io/@bkomuves/HyBL5V5xMe) (mostly following the paper).

Here some details are different from the paper though, and our implementation 
is noticeably more efficient in practice.

### Different versions

We provide several slightly different versions with different tradeoffs

- `dynamic/v1`: This is a version mostly following the Dynark paper (TODO: sparse convolution algorithm)
- `dyanmic/v2`: This introduces some of our optimizations: 
     - the projection term updates become essentially free
     - we reorder the circuit so that the changes are in a subgroup (TODO)
