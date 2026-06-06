
import groth16/test_proof
import groth16/files/export_json

#-------------------------------------------------------------------------------

proc exampleProveAndVerify() =
  let zkey_fname : string = "./build/product.zkey"
  let wtns_fname : string = "./build/product.wtns"
  let (_, proof) = testProveAndVerify( zkey_fname, wtns_fname)

  exportPublicIO( "./build/nim_public.json" , proof )
  exportProof(    "./build/nim_proof.json"  , proof )

#-------------------------------------------------------------------------------

when isMainModule:
  exampleProveAndVerify()