package core

import (
	"math"
)

const (
	GOPBRT_VERSION = "0.0.0"
	PBRT_VERSION   = "2.0.0"
)

const (
	INV_TWOPI = 1.0 / (2.0 * math.Pi)
	INV_PI = 1.0 / math.Pi
)

type Options struct {
	NumCores                                int
	QuickRender, Quiet, Verbose, OpenWindow, Debug bool
	ImageFile                               string
}

// Global Inline Functions
func Lerp(t, v1, v2 float64) float64 {
	return (1.0-t)*v1 + t*v2
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

func Maxi(v1, v2 int) int {
    if v1 > v2 { return v1 }
    return v2
}


func Mini(v1, v2 int) int {
    if v1 < v2 { return v1 }
    return v2
}

func Mod(a, b int) int {
	n := int(a / b)
	a -= n * b
	if a < 0 {
		a += b
	}
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
	return v + 1
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

func Quadratic(A, B, C float64) (ok bool, t0, t1 float64) {
    // Find quadratic discriminant
    discrim := B * B - 4.0 * A * C
    if (discrim < 0.0) {
         return false, 0.0, 0.0
    }
    rootDiscrim := math.Sqrt(discrim)

    // Compute quadratic _t_ values
    var q float64
    if B < 0.0 { 
        q = -0.5 * (B - rootDiscrim)
    } else {
        q = -0.5 * (B + rootDiscrim)
    }
    t0 = q / A
    t1 = C / q
    if t0 > t1 {
        t0, t1 = t1, t0
    }
    return true, t0, t1
}

func Xor(a, b bool) bool {
    if (a && !b) || (!a && b) { return true }
    return false
}
func NumSystemCores() int {
	return 1
}

func AtomicAdd(dest *float64, delta float64) {
    // TODO: port this function
    *dest = *dest + delta
}
func AtomicAddf(dest *float32, delta float32) {
    // TODO: port this function
    *dest = *dest + delta
}
