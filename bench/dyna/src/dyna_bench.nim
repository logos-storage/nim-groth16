
{. warning[UnusedImport]:off .}

import sugar
import std/strutils
import std/sequtils
import std/tables
import std/options
import std/random
import std/parseopt
import std/strformat

import system
import std/os
import std/osproc
import std/streams
import std/paths
import std/dirs
import std/times

import marshal

import taskpools

import constantine/named/properties_fields

import groth16/bn128
import groth16/bn128/arrays
import groth16/files/zkey
import groth16/files/witness
import groth16/files/r1cs

import groth16/prover
import groth16/prover/shared
import groth16/prover/types
import groth16/verifier

import groth16/partial/types as ptypes
import groth16/partial/precalc
import groth16/partial/finish

import groth16/dynamic/permute
import groth16/dynamic/v1/types
import groth16/dynamic/v1/setup
import groth16/dynamic/v1/preprocess
import groth16/dynamic/v1/finish

import groth16/dynamic/permute
import groth16/dynamic/v3/types
import groth16/dynamic/v3/setup
import groth16/dynamic/v3/preprocess
import groth16/dynamic/v3/finish

import groth16/zkey_types
import groth16/misc as gmisc

import circom_witnessgen/types as wgentypes
import circom_witnessgen/input_json
import circom_witnessgen/graph
import circom_witnessgen/load
import circom_witnessgen/dependencies
import circom_witnessgen/witness
import circom_witnessgen/partial
import circom_witnessgen/export_wtns

#-------------------------------------------------------------------------------

const circomExe    = "circom"  
const snarkjsExe   = "snarkjs" 
const malleableExe = "snarkjs-malleable"
const graphExe     = "build-circuit"
const ptauPath     = "~/zk/ptau"

#-------------------------------------------------------------------------------

type 

  Params = object
    fullSize       : int
    witnessUpdate  : int
    imageUpdate    : int

  Groth16Timings = object
    proofTime      : float64
    verifyTime     : float64

  PartialTimings = object
    partialTime    : float64
    finishTime     : float64

  DynaV1Timings = object
    setupTime      : float64
    preprocessTime : float64
    finishTime     : float64

  DynaV3Timings = object
    setupTime      : float64
    preprocessTime : float64
    finishTime     : float64

  Timings = object
    tgtParams   : Params
    realParams  : Params
    witnessgen  : float64
    groth16     : Groth16Timings
    partial     : PartialTimings
    dynaV1      : DynaV1Timings
    dynaV3      : DynaV3Timings

#-------------------------------------------------------------------------------

proc echoPhase(text: string) =
  echo "\n------------------------------------------------------------"
  echo (text & "\n")

#-------------------------------------------------------------------------------

template measuringTime*(doPrint: bool, text: string, code: untyped): float64 =
  block:
    let t0 = epochTime()
    code
    let elapsed = epochTime() - t0
    if doPrint:
      let elapsedStr = elapsed.formatFloat(format = ffDecimal, precision = 4)
      echo ( text & " took " & elapsedStr & " seconds" )
    elapsed

#-------------------------------------------------------------------------------
#
# fucking stupid Nim shitstorm
# <https://github.com/nim-lang/Nim/issues/9953>
#

proc nimExecProcess(command: string, args: openArray[string]) = 
  let opts = { poUsePath, poEchoCmd, poStdErrToStdOut, poParentStreams } 
  let process = osproc.startProcess( command = command, args = args, env = nil, options = opts )
  discard process.waitForExit()

proc nimExecCmd(command: string) = 
  let opts = { poUsePath, poEchoCmd, poStdErrToStdOut }   # , poParentStreams } 
  discard osproc.execCmdEx( command = command, options = opts )

#-------------------------------------------------------------------------------

proc exportJsonFr(fname: string, fields: seq[(string,F)]) = 

  let f = open(fname, fmWrite)
  defer: f.close()
  for i,pair in fields.pairs():
    let prefix = if i == 0: "{ " else: ", "
    let ln = prefix & "\"" & pair[0] & "\": \"" & toDecimalFr(pair[1]) & "\""
    f.writeLine(ln)
  f.writeLine("}")

func inputSeqToTable( fields: seq[(string,Fr)] ) : Table[string,seq[Fr]] =
  var table : Table[string,seq[Fr]] = initTable[string,seq[Fr]]()
  for pair in fields:
    let list: seq[Fr] = @[ pair[1] ]
    table[pair[0]] = list
  return table

#-------------------------------------------------------------------------------

func nextNumber( what: int, available: seq[int] ): int = 
  var res: int = 0
  for y in available:
    if what <= y:
      res = y
      break
  return res

func guessPtauFile( nconstr: int ): string = 
  let log2 = ceilingLog2( nConstr )
  let available = @[12,14,16,18,20,21,22]
  let next = nextNumber( log2 , available )
  let fname = "powersOfTau28_hez_final_" & ($next) & ".ptau"
  return (ptauPath & "/" & fname)

#-------------------------------------------------------------------------------

proc benchmarkTarget(pool: Taskpool, bigSize: int, smallSize: int) =

  echo "----------------------------------------"
  echo "target size: " & ($bigSize) & "/" & ($smallSize)

  let updateRounds = (smallSize - 4) div 9
  let staticRounds = (bigSize - updateRounds*9 - 4) div 9

  let origDir    = string(paths.getCurrentDir())
  let buildDir   = origDir  & "/build"
  let timingsDir = origDir  & "/timings"
  let circuitDir = origDir  & "/circuit"   
  let circomFile = buildDir & "/main.circom"
  let graphFile  = buildDir & "/main.graph"

  createDir(buildDir)
  setCurrentDir(buildDir)

  # --- create the `main.circom` file ---

  echoPhase("creating `main.circom`...")
  let argText  = "(" & ($staticRounds) & "," & ($updateRounds) & ")"
  let mainText = "pragma circom 2.1.1;\n" &
                 "include \"parametric_bench.circom\";\n" &
                 "component main { public [staticInp,updateInp] } = Benchmark" & argText & ";\n"
  echo "writing " & circomFile
  writeFile(circomFile, mainText)

  # --- compile the circuit ---

  echoPhase("compiling the circuit...")
  nimExecProcess( circomExe, ["--O2", "--r1cs", "--wasm", "-l"&circuitDir, circomFile] )

  # --- extract the witness computation graph ---

  echoPhase("generating the witness...")
  nimExecProcess( graphExe, ["--O2", circomFile, graphFile, "-l"&circuitDir] )

  # --- generate partial and full input files ---

  echoPhase("writing input and partial input json...")
  let x = randFr()
  let y = randFr()
  let z = randFr()
  let fullInputsSeq    = @[ ("staticInp",x), ("secretInp", y), ("updateInp",z) ]
  let partialInputsSeq = @[ ("staticInp",x), ("secretInp", y)                  ]
  exportJsonFr( "input.json"   , fullInputsSeq    )
  exportJsonFr( "partial.json" , partialInputsSeq )

  # --- generate the witness ---
  
  echoPhase("generating the witness...")
  let witnessJsDir = buildDir & "/main_js" 
  setCurrentDir(witnessJsDir)
  nimExecProcess( "node", ["generate_witness.js", "main.wasm", "../input.json", "../main.wtns"] )
  setCurrentDir(buildDir)

  # --- extracting the partial witness ---

  echoPhase("extracting partial witness...")
  let graph = loadGraph(graphFile)
  let fullInputsTable    : Table[string,seq[F]] = inputSeqToTable(    fullInputsSeq )
  let partialInputsTable : Table[string,seq[F]] = inputSeqToTable( partialInputsSeq )

  var ourFullWitness : wgentypes.Witness
  let witnessgen_time = measuringTime(true, "witness generation (using the witness graph)"):
    ourFullWitness = generateWitness(graph, fullInputsTable)

  let partialWitness1 : seq[Option[F]] = generatePartialWitness(graph, partialInputsTable)
  let partialWitness2 : ptypes.PartialWitness = ptypes.PartialWitness(values: partialWitness1)
  exportWitness( buildDir & "/our.wtns", ourFullWitness )

  # --- permuting the R1CS rows ---

  echoPhase("permuting the R1CS rows...")
  let witnessDeltaMask = notBoolSeq(partialWitnessMask(partialWitness2))
  let origR1CSFile     = buildDir & "/main.r1cs"
  let permutedR1CSFile = buildDir & "/permuted.r1cs"
  let origR1CS = parseR1CS(origR1CSFile)
  let subgroupSize = exportPermutedR1CS( permutedR1CSFile, origR1CS, witnessDeltaMask ) 
  let permutedR1CS = parseR1CS(permutedR1CSFile)
  printR1CSMetaData(permutedR1CS)

  # --- snarkjs trusted setup ---

  echoPhase("creating Groth16 trusted setup via snarkjs...")

  let ptauFile = guessPtauFile( permutedR1CS.nConstr )
  echo "ptau file = `" & ptauFile & "`"

  let cmd1 = "NODE_OPTIONS=\"--max-old-space-size=8192\" " & malleableExe & " groth16 setup permuted.r1cs " & ptauFile & " tmp_0000.zkey"
  let cmd2 = "echo \"some_entropy_17528v3a7rcawcsyiur\" | NODE_OPTIONS=\"--max-old-space-size=8192\" " & malleableExe & " zkey contribute tmp_0000.zkey tmp_0001.zkey --name=\"1st Contributor Name\""
  let cmd3 = snarkjsExe & " zkey export verificationkey permuted.zkey permuted_verification_key.json"

  nimExecCmd( cmd1 )
  nimExecCmd( cmd2 )
  removeFile("tmp_0000.zkey")
  moveFile("tmp_0001.zkey", "permuted.zkey")
  nimExecCmd( cmd3 )

  # --- loading zkey (and witness) --- 

  echoPhase("loading zkey file...")
  let zkeyFile = buildDir & "/permuted.zkey"
  let permutedZKey = parseZKey( zkeyFile )
  let permutedVKey = extractVKey( permutedZKey )

  let circomFullWitness = parseWitness("main.wtns")

  # --- computing the delta image ---

  let (partialAB, deltaImages) = computeDeltaImages(permutedZKey, partialWitness2, false)
  let imageAB   = deltaImages.imageAB
  let imageSize = countTrues(imageAB)
  echo "number of constraints = " & ($permutedR1CS.nConstr)
  echo "size of imageAB       = " & ($imageSize)

  # --- standard Groth16 proof --- 

  echoPhase("normal Groth16 proving...")
  var theProof: Proof
  let standardProofTime = measuringTime(true, "standard Groth16 proof"):
    theProof = generateProof( permutedZKey, circomFullWitness, pool, false )

  let verifyTime = measuringTime(true, "verifying the proof"):
    let ok = verifyProof( permutedVKey, theProof )
    echo "the proof verified ok = " & $ok
    assert(ok, "VERIFYING THE PROOF FAILED!")

  let groth16_timings = 
        Groth16Timings( proofTime  : standardProofTime ,
                        verifyTime : verifyTime        )

  # echo $standardProoftime
  # echo $verifyTime

  # --- partial proof --- 

  echoPhase("partial Groth16 proofs...")

  var partialProof: PartialProof
  let partialProofTime = measuringTime(true, "partial Groth16 proof"):
    partialProof = generatePartialProof( permutedZKey, partialWitness2, pool, false )

  var finishedPartialProof: Proof
  let partialFinishTime = measuringTime(true, "finishing partial Groth16 proof"):
    finishedPartialProof = finishPartialProof( permutedZKey, circomFullWitness, partialProof, pool, false )
  
  block:
    let ok = verifyProof( permutedVKey, finishedPartialProof )
    assert(ok, "PARTIAL PROOF FAILED TO VERIFY!")

  let partial_timings =
        PartialTimings( partialTime : partialProofTime  , 
                        finishTime  : partialFinishTime )

  # --- dyna proof V1 --- 

  echoPhase("dnyamic Groth16 proofs (V1, our version)...")

  var v1_setup: DynaSetupV1
  let v1_SetupTime = measuringTime(true, "dynamic setup V1"):
    v1_setup = dynaSetupV1FromZKey( permutedZkey, pool )

  var v1_preproof: DynaPreProofV1
  let v1_PreProofTime = measuringTime(true, "dynamic preprocess V1"):
    v1_preproof = dynaPreProofV1( permutedZKey, v1_setup, partialWitness2, pool, false )

  var v1_proof: Proof
  let v1_FinishTime = measuringTime(true, "dynamic finish V1"):
    v1_proof = finishDynaProofV1( permutedZKey, circomFullWitness, v1_preproof, pool, false )

  block:
    let ok = verifyProof( permutedVKey, v1_proof )
    assert(ok, "V1 DYNAMIC PROOF FAILED TO VERIFY!")

  let dynaV1_timings = 
        DynaV1Timings( setupTime      : v1_SetupTime    ,
                       preprocessTime : v1_PreProofTime ,
                       finishTime     : v1_FinishTime   )

  # --- dyna proof V3 --- 

  echoPhase("dnyamic Groth16 proofs (V3, our version)...")

  echo "subgroup size = " & ($subgroupSize)

  var v3_setup: DynaSetupV3 
  let v3_SetupTime = measuringTime(true, "dynamic setup V3"):
    v3_setup = dynaSetupV3FromZKey( permutedZkey, subgroupSize, pool )

  var v3_preproof: DynaPreProofV3 
  let v3_PreProofTime = measuringTime(true, "dynamic preprocess V3"):
    v3_preproof = dynaPreProofV3( permutedZKey, v3_setup, partialWitness2, pool, false )

  var v3_proof: Proof
  let v3_FinishTime = measuringTime(true, "dynamic finish V3"):
    v3_proof = finishDynaProofV3( permutedZKey, circomFullWitness, v3_preproof, pool, false )

  block:
    let ok = verifyProof( permutedVKey, v3_proof )
    assert(ok, "V3 DYNAMIC PROOF FAILED TO VERIFY!")

  let dynaV3_timings = 
        DynaV3Timings( setupTime      : v3_SetupTime    ,
                       preprocessTime : v3_PreProofTime ,
                       finishTime     : v3_FinishTime   )

  # --- export timings --- 

  echoPhase("export timings...")

  let tgtParams = Params( fullSize      : bigSize   ,
                          witnessUpdate : smallSize ,
                          imageUpdate   : smallSize )

  let realParams = Params( fullSize      : permutedR1CS.nConstr         ,
                           witnessUpdate : countTrues(witnessDeltaMask) , 
                           imageUpdate   : imageSize                    )

  let timings = Timings( tgtParams  : tgtParams       , 
                         realParams : realParams      ,
                         witnessgen : witnessgen_time ,
                         groth16    : groth16_timings ,
                         partial    : partial_timings ,
                         dynaV1     : dynaV1_timings  ,
                         dynaV3     : dynaV3_timings  )

  createDir(timingsDir)
  let fname = timingsDir & "/timings_" & ($bigSize) & "_" & ($smallSize) & ".json"
  writeFile( fname, $$timings )      # marshal / $$ is generic JSON serialization

  # --- exit ---

  echoPhase("done.")
  setCurrentDir(origDir)

#-------------------------------------------------------------------------------

when isMainModule:
  setStdIoUnbuffered()

  let nthreads = 1
  var pool  = Taskpool.new(nthreads)
 
  # benchmarkTarget(pool, 2048, 128)
  # benchmarkTarget(pool, 8192, 512)

  let smallSizes = @[128,256,512,1024]
  let bigSizes   = @[2048,4096,8192,16384,32768]
  for big in bigSizes:
    for small in smallSizes:
      benchmarkTarget(pool, big, small)

