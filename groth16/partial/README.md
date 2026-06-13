
Partial proofs
--------------

If a large part of your witness is constant (or changing infrequently),
you can precalculate a significant portion of the proof.

See [this forum post](https://forum.research.logos.co/t/speeding-up-rln-proofs/663) 
for more details.

This is essentially free (the overhead is minimal, so already when doing 2 proofs
with same shared part of the witness, it's worth do it).

A more advanced version, but also with a more heave trade-off, is developed in
[`groth16/dynamic`](../dynamic/).

