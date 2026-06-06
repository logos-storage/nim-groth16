
version     = "0.1.1"
author      = "Balazs Komuves"
description = "Groth16 proof system"
license     = "MIT OR Apache-2.0"

skipDirs    = @["groth16/example"]
binDir      = "build"
namedBin    = {"cli/cli_main": "nim-groth16"}.toTable()
installExt  = @["nim"]

requires "nim >= 2.2.0"
requires "https://github.com/status-im/nim-taskpools >= 0.0.5"
requires "https://github.com/mratsim/constantine"
# requires "https://github.com/mratsim/constantine#bc3845aa492b52f7fef047503b1592e830d1a774"
