
#
# the container format used by `circom` / `snarkjs`
# see <https://github.com/iden3/binfileutils>
#
# format:
# =======
#
# global header:
# --------------
#   magic              : word32
#   version            : word32
#   number of sections : word32
#
# for each section:
# -----------------
#   section id         : word32 
#   section size       : word64
#   section data       : <section_size> number of bytes
#

#-------------------------------------------------------------------------------

import std/streams

import sugar

import constantine/math/arithmetic
import constantine/math/io/io_bigints

#-------------------------------------------------------------------------------

type 
  SectionCallback*[T] = proc (stream: Stream, sectId: int, sectLen: int, user: var T) {.closure.}

#-------------------------------------------------------------------------------

func magicWord(magic: string): uint32 = 
  assert( magic.len == 4, "magicWord: expecting a string of 4 characters" )
  var w : uint32 = 0 
  for i in 0..3:
    let a = uint32(ord(magic[i])) 
    w += a shl (8*i)
  return w

#-------------------------------------------------------------------------------

proc parsePrimeField*( stream: Stream ) : (int, BigInt[256]) = 
  let n8p = int( stream.readUint32() )
  assert( n8p <= 32 , "at most 256 bit primes are allowed" )
  var p_bytes : array[32, uint8] 
  discard stream.readData( addr(p_bytes), n8p )
  var p : BigInt[256]
  unmarshal(p, p_bytes, littleEndian);
  return (n8p, p)

#-------------------------------------------------------------------------------

proc readSection[T] ( expectedMagic: string
                    , expectedVersion: int
                    , stream: Stream
                    , user: var T
                    , callback: SectionCallback[T] 
                    , filt: proc(_: int): bool {.raises: [IOError, OSError].}
                    ) {.raises: [IOError, OSError].} =

  let sectId  = int( stream.readUint32() )
  let sectLen = int( stream.readUint64() )
  let oldpos = stream.getPosition()
  if filt(sectId):
    callback(stream, sectId, sectLen, user)
  stream.setPosition(oldpos + sectLen)

#-------------------------------------------------------------------------------

proc parseContainer*[T] ( expectedMagic: string
                        , expectedVersion: int
                        , fname: string
                        , user: var T
                        , callback: SectionCallback[T] 
                        , filt: proc(_: int): bool {.raises: [IOError, OSError].}
                        ) {.raises: [IOError, OSError].} =

  let stream = newFileStream(fname, mode = fmRead)
  defer: stream.close()

  let magic = stream.readUint32()
  assert( magic == magicWord(expectedMagic) , "not a `" & expectedMagic & "` file" )
  let version = stream.readUint32()
  assert( version == uint32(expectedVersion) , "not a version " & ($expectedVersion) & " `" & expectedMagic & "` file" )
  let nsections = stream.readUint32()
  # echo("number of sections = ",nsections)

  for i in 1..nsections:
    readSection(expectedMagic, expectedVersion, stream, user, callback, filt)

#-------------------------------------------------------------------------------

