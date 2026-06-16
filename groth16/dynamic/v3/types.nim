
{.push raises:[].}

#import constantine/named/properties_fields

import groth16/bn128
import groth16/math/domain
#import groth16/math/ntt
#import groth16/math/poly

import groth16/partial/types as ptypes
import groth16/dynamic/types as dtypes

export ptypes
export dtypes

#-------------------------------------------------------------------------------

type 

  # things we can compute at circuit setup time
  DynaSetupV3* = object
    imageSubgroup*  : Subgroup     # the subgroup the changes fall on
    weightVec*      : seq[F]       # the values `W_k` (it's surprisingly slow to re-compute in the finish part!)
    pointsDeltaLZ*  : seq[G1]      # the points `delta^-1 * L_i(tau) * Z(tau) * g1` where `Z(x) = x^N-1`
    wConvDeltaLZ*   : seq[G1]      # the convolution of `W` and `pointsDeltaLZ`
    miniDiagPoints* : seq[G1]      # corrections for subgroup convolution...

  # things we can compute from the partial witness
  DynaPreprocessV3* = object
    unified*      : seq[G1]       # the points `X_j := sum_k ( B_{kj} * U_k + A_{kj} * V_k )

  # "pre-proof" of the Dynark paper, which we can finish 
  DynaPreProofV3* = object
    partialProof*   : PartialProof       # linear part of the partial proof
    deltaImages*    : DeltaImages        # image of the witness delta under A and B (where can update happen)
    dynaSetup*      : DynaSetupV3        # circuit-setup time calculations
    dynaPreprocess* : DynaPreprocessV3   # "Dynark"-style preprocessing
    
#-------------------------------------------------------------------------------
