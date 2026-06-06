

import std/[times,os]
import strformat

import taskpools

import groth16/prover
import groth16/verifier
import groth16/files/witness
import groth16/files/r1cs
import groth16/files/zkey
import groth16/zkey_types
import groth16/fake_setup

func seconds(x: float): string = fmt"{x:.4f} seconds"

#-------------------------------------------------------------------------------

proc testProveAndVerify*( zkey_fname, wtns_fname: string): (VKey,Proof) = 

  echo("parsing witness & zkey files...")
  let witness = parseWitness( wtns_fname)
  let zkey    = parseZKey( zkey_fname)

  echo("generating proof...")
  let start = cpuTime()
  var pool = Taskpool.new()
  let proof = generateProof( zkey, witness ,pool)
  let elapsed = cpuTime() - start
  echo("proving took ",seconds(elapsed))

  echo("verifying the proof...")
  let vkey = extractVKey( zkey)
  let ok   = verifyProof( vkey, proof )
  echo("verification succeeded = ",ok)

  return (vkey,proof)

#-------------------------------------------------------------------------------

proc testFakeSetupAndVerify*( r1cs_fname, wtns_fname: string, flavour=Snarkjs): (VKey,Proof) = 
  echo("trusted setup flavour = ",flavour)

  echo("parsing witness & r1cs files...")
  let witness = parseWitness( wtns_fname)
  let r1cs    = parseR1CS( r1cs_fname)

  echo("performing fake trusted setup...")
  let start1 = cpuTime()
  let zkey = createFakeCircuitSetup( r1cs, flavour=flavour )
  let elapsed1 = cpuTime() - start1
  echo("fake setup took ",seconds(elapsed1))

  # printCoeffs(zkey.coeffs)

  echo("generating proof...")
  let vkey  = extractVKey( zkey)

  let start = cpuTime()
  var pool = Taskpool.new()
  let proof = generateProof( zkey, witness ,pool)
  let elapsed = cpuTime() - start
  echo("proving took ",seconds(elapsed))

  echo("verifying the proof...")
  let ok = verifyProof( vkey, proof )
  echo("verification succeeded = ",ok)

  return (vkey,proof)
