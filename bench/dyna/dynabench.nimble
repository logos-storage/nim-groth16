version       = "0.0.0"
author        = "Balazs Komuves"
description   = "benchmarking dynamic Groth16 proofs"
license       = "MIT or Apache-2.0"
srcDir        = "src"
bin           = @["dyna_bench"]

requires "nim >= 2.0.0"
requires "https://github.com/mratsim/constantine#bc3845aa492b52f7fef047503b1592e830d1a774"
requires "https://github.com/logos-storage/circom-witnessgen#461a7b14d4c2bf76f1f94cc3b91d2beb9d5652fa"
requires "https://github.com/logos-storage/nim-groth16#73b5ae2734050d64157afbc29d364345ff0ec211"