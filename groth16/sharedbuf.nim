#
# SharedBuf[T] — non-owning (ptr, len) view over a contiguous buffer.
#
# Purpose: cross a `Taskpool.spawn` boundary safely under any GC mode
# (refc / arc / orc) without sending a `seq[T]` through the spawn closure.
#
# Contract:
#   - The pointed-to memory is owned by the caller. The caller MUST keep
#     the underlying storage alive until every worker holding the view
#     has finished (i.e. all FlowVars from `pool.spawn` have been `sync`ed).
#   - Element type `T` must have no GC fields (workers read/write payload
#     bytes through a raw pointer; the seq's GC machinery is never touched).
#

{.push raises: [], gcsafe.}

template makeUncheckedArray[T](p: ptr T): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](p)

type SharedBuf*[T] = object
  payload*: ptr UncheckedArray[T]
  len*:     int

func view*[T](_: type SharedBuf, v: openArray[T]): SharedBuf[T] =
  if v.len > 0:
    SharedBuf[T](payload: makeUncheckedArray(addr v[0]), len: v.len)
  else:
    default(SharedBuf[T])

template checkIdx(v: SharedBuf, i: int) =
  doAssert i >= 0 and i < v.len

func `[]`*[T](v: SharedBuf[T], i: int): var T =
  v.checkIdx(i)
  v.payload[i]

template toOpenArray*[T](v: SharedBuf[T]): var openArray[T] =
  v.payload.toOpenArray(0, v.len - 1)

template toOpenArray*[T](v: SharedBuf[T], s, e: int): var openArray[T] =
  v.toOpenArray().toOpenArray(s, e)

{.pop.}
