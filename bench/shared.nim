
import strformat
import times, strutils

#-------------------------------------------------------------------------------

func seconds*(x: float): string = fmt"{x:.4f} seconds"

func quoted*(s: string): string = fmt"`{s:s}`"

template withMeasureTime*(doPrint: bool, text: string, code: untyped) =
  block:
    if doPrint:
      let t0 = epochTime()
      code
      let elapsed = epochTime() - t0
      let elapsedStr = elapsed.formatFloat(format = ffDecimal, precision = 4)
      echo ( text & " took " & elapsedStr & " seconds" )
    else:
      code

#-------------------------------------------------------------------------------

