
"Semi-dynamic" Groth16 proofs
-----------------------------

This is loosely based on the [Dynark paper](https://eprint.iacr.org/2025/1897):

- _"Dynark: Making Groth16 Dynamic"_ by Tianyu Zhang, Yupeng Ouyang and Yupeng Zhang

See also [this write-up](https://hackmd.io/@bkomuves/HyBL5V5xMe) (mostly following the paper).

Here some details are different from the paper though, and our implementation 
is noticeably more efficient in practice.

### Different versions

We provide several slightly different versions with different tradeoffs

- `dynamic/v1`: This is a version mostly following the Dynark paper 
     - TODO: sparse convolution algorithm
- `dyanmic/v2`: This introduces some of our optimizations: 
     - the projection term updates become essentially free
     - micro-optim: concatenate all the `pi_C` MSMs into one
- `dyanmic/v3`: 
     - we reorder the circuit so that the changes are in a subgroup
     - so cross-terms becomes `O(D*log(D)`, and also just 1 MSM instead of 2
     - so the updates now are essentially optimal
     - TODO: half of the very slow group FFTs (during the slow preprocessing step) 
       also becomes `O(D*log)`, so that's essenitally a 2x speedup

