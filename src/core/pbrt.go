package pbrt

import (
    "math"
)

const (
	GOPBRT_VERSION = "0.0.0"
	PBRT_VERSION = "2.0.0"
)

type Options struct {
    NumCores int
    QuickRender, Quiet, Verbose, OpenWindow bool
    ImageFile string
}

// Global Inline Functions
func Lerp(t, v1, v2 float64) float64 {
    return (1.0 - t) * v1 + t * v2
}

func Clamp(f, fmin, fmax float64) float64 {
    if f < fmin {
        return fmin
    } else if f > fmax {
        return fmax
    }
    return f
}

func Clampi(val, low, high int) int {
    if val < low {
    	return low
   	} else if val > high { 
   		return high
	} else {
		return val
	}	 
}

func Mod(a, b int) int {
    n := int(a/b)
    a -= n*b
    if a < 0 { a += b }
    return a
}

func Radians(angle float64) float64 {
    return (angle * math.Pi) / 180.0
}

func Degrees(rad float64) float64 {
    return (180.0 / math.Pi) * rad
}

func Log2(x float64) float64 {
    return math.Log2(x)
}

func Floor2Int(v float64) int {
	return int(math.Floor(v))
}

func Log2Int(v float64) int {
    return Floor2Int(Log2(v))
}

func IsPowerOf2(v int) bool {
    return (v != 0) && ((v & (v - 1)) == 0)
}

func RoundUpPow2(v uint32) uint32 {
    v--
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v+1
}

func Round2Int(val float64) int {
    return Floor2Int(val + 0.5)
}

func Float2Int(val float64) int {
    return int(val)
}

func Ceil2Int(val float64) int {
    return int(math.Ceil(val))
}

func NumSystemCores() int {
    return 1
}

func ParseFile(filename string) bool {
    return false
}