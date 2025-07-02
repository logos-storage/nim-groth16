
import sugar
import std/strutils
import std/sequtils
import std/os
import std/cpuinfo
import std/parseopt
import std/times
import std/options
# import strformat

import taskpools

import groth16/prover
import groth16/verifier
import groth16/files/witness
import groth16/files/r1cs
import groth16/files/zkey
import groth16/files/export_json
import groth16/zkey_types
import groth16/fake_setup
import groth16/misc

#-------------------------------------------------------------------------------

proc printHelp() =
  echo "usage:"
  echo "$ nim-groth16 [options] --output=proof_input.json --circom=proof_main.circom"
  echo ""
  echo "available options:"
  echo " -h, --help                      : print this help"
  echo " -v, --verbose                   : verbose output"
  echo " -d, --debug                     : debug output"
  echo " -j, --nthreads                  : number of CPU threads"
  echo " -t, --time                      : print time measurements"
  echo " -p, --prove                     : create a proof"
  echo " -y, --verify                    : verify a proof"
  echo " -u, --setup                     : perform (fake) trusted setup"
  echo " -n, --nomask                    : don't use random masking for full ZK"
  echo " -z, --zkey   = <circuit.zkey>   : the `.zkey` file"
  echo " -w, --wtns   = <circuit.wtns>   : the `.wtns` file" 
  echo " -r, --r1cs   = <circuit.r1cs>   : the `.r1cs` file" 
  echo " -o, --output = <proof.json>     : the proof file"
  echo " -i, --io     = <public.json>    : the public input/output file" 

#-------------------------------------------------------------------------------

type Config = object
  zkey_file:    string
  r1cs_file:    string
  wtns_file:    string
  output_file:  string
  io_file:      string
  verbose:      bool
  debug:        bool
  measure_time: bool
  do_prove:     bool
  do_verify:    bool
  do_setup:     bool
  no_masking:   bool
  nthreads:     int

const dummyConfig = 
  Config( zkey_file:    ""
        , r1cs_file:    ""
        , wtns_file:    ""
        , output_file:  ""
        , io_file:      ""
        , verbose:      false
        , measure_time: false
        , do_prove:     false
        , do_verify:    false
        , do_setup:     false
        , no_masking:   false
        , nthreads:     0
        )

proc printConfig(cfg: Config) =
  echo ""
  echo "zkey file        = " & ($cfg.zkey_file)
  echo "witness file     = " & ($cfg.wtns_file)
  echo "r1cs file        = " & ($cfg.r1cs_file)
  echo "proof file       = " & ($cfg.output_file)
  echo "public i/o file  = " & ($cfg.io_file)

#-------------------------------------------------------------------------------

proc parseCliOptions(): Config =

  var argCtr: int = 0
  var swCtr:  int = 0

  var cfg: Config = dummyConfig

  for kind, key, value in getOpt():
    case kind

    # Positional arguments
    of cmdArgument:
      # echo ("arg #" & $argCtr & " = " & key)
      argCtr += 1

    # Switches
    of cmdLongOption, cmdShortOption:
      swCtr += 1
      case key

      of "h", "help"             : printHelp()
      of "v", "verbose"          : cfg.verbose        = true
      of "d", "debug"            : cfg.debug          = true
      of "j", "nthreads"         : cfg.nthreads       = parseInt(value)
      of "t", "time"             : cfg.measure_time   = true
      of "p", "prove"            : cfg.do_prove       = true
      of "y", "verify"           : cfg.do_verify      = true
      of "u", "setup"            : cfg.do_setup       = true
      of "n", "nomask"           : cfg.no_masking     = true
      of "o", "output"           : cfg.output_file    = value
      of "r", "r1cs"             : cfg.r1cs_file      = value
      of "z", "zkey"             : cfg.zkey_file      = value
      of "w", "wtns", "witness"  : cfg.wtns_file      = value
      of "i", "io", "input"      : cfg.io_file        = value
      else:
        echo "Unknown option: ", key
        echo "use --help to get a list of options"
        quit()

    of cmdEnd:
      discard  

  if swCtr==0 and argCtr==0:
    printHelp()
    quit()

  if cfg.nthreads <= 0: 
    cfg.nthreads = countProcessors()

  return cfg

#-------------------------------------------------------------------------------

#[
proc testProveAndVerify*( zkey_fname, wtns_fname: string): (VKey,Proof) = 

  echo("parsing witness & zkey files...")
  let witness = parseWitness( wtns_fname)
  let zkey    = parseZKey( zkey_fname)

  echo("generating proof...")
  let start = cpuTime()
  let proof = generateProof( zkey, witness )
  let elapsed = cpuTime() - start
  echo("proving took ",seconds(elapsed))

  echo("verifying the proof...")
  let vkey = extractVKey( zkey)
  let ok   = verifyProof( vkey, proof )
  echo("verification succeeded = ",ok)

  return (vkey,proof)
]#

#-------------------------------------------------------------------------------

proc cliMain(cfg: Config) =

  var wtns:  Witness
  var zkey:  ZKey
  var r1cs:  R1CS
  var proof: Proof

  if not (cfg.wtns_file == ""):
    echo("\nparsing witness file " & quoted(cfg.wtns_file))
    withMeasureTime(cfg.measure_time,"parsing the witness"):
      wtns = parseWitness(cfg.wtns_file)
  
  if not (cfg.zkey_file == ""):
    echo("\nparsing zkey file " & quoted(cfg.zkey_file))
    withMeasureTime(cfg.measure_time,"parsing the zkey"):
      zkey = parseZKey(cfg.zkey_file)
  
  if not (cfg.r1cs_file == ""):
    echo("\nparsing r1cs file " & quoted(cfg.r1cs_file))
    withMeasureTime(cfg.measure_time,"parsing the r1cs"):
      r1cs = parseR1CS(cfg.r1cs_file)
  
  if cfg.do_setup:
    if not (cfg.zkey_file == ""):
      echo("\nwe are doing a fake trusted setup, don't specify the zkey file!")   
      quit()
    if (cfg.r1cs_file == ""):
      echo("\nerror: r1cs file is required for the fake setup!")
      quit()
    echo("\nperforming fake trusted setup...")
    withMeasureTime(cfg.measure_time,"fake setup"):
      zkey = createFakeCircuitSetup( r1cs, flavour=Snarkjs )
  
  if cfg.debug:
    printGrothHeader(zkey.header)
    # debugPrintCoeffs(zkey.coeffs)

  if cfg.do_prove:
    if (cfg.wtns_file=="") or (cfg.zkey_file=="" and cfg.do_setup==false):
      echo("cannot prove: missing witness and/or zkey file!")      
      quit()
    else:
      echo("generating proof...")
      let print_timings = cfg.measure_time and cfg.verbose
      var pool = Taskpool.new(cfg.nthreads)
      withMeasureTime(cfg.measure_time,"proving"):
        if cfg.no_masking:
          proof = generateProofWithTrivialMask(zkey, wtns, pool, print_timings)
        else:
          proof = generateProof(zkey, wtns, pool, print_timings)
    
      if not (cfg.output_file == ""):
        echo("exporting the proof to " & quoted(cfg.output_file))
        exportProof( cfg.output_file, proof )
      if not (cfg.io_file == ""):
        echo("exporting the public IO to " & quoted(cfg.io_file))
        exportPublicIO( cfg.io_file, proof )

  if cfg.do_verify:
    if (cfg.zkey_file == "" and cfg.do_setup==false):
      echo("cannot verify: missing vkey (well, zkey)")      
      quit()
    else:
      let vkey = extractVKey( zkey)
      echo("\nverifying the proof...")
      var ok : bool
      withMeasureTime(cfg.measure_time,"verifying"):
        ok = verifyProof( vkey, proof )
        echo("verification succeeded = ",ok)
    
  echo("")

#-------------------------------------------------------------------------------

when isMainModule:
  let cfg = parseCliOptions()
  if cfg.verbose: printConfig(cfg)
  cliMain(cfg)
