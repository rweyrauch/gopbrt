package pbrt

import (
	"math"
	//"strings"
)

type (
	Filter interface {
		Evaluate(x, y float64) float64
		XWidth() float64
		YWidth() float64
		InvXWidth() float64
		InvYWidth() float64
	}

	FilterData struct {
		xWidth, yWidth       float64
		invXWidth, invYWidth float64
	}

	BoxFilter struct {
		FilterData
	}

	GaussianFilter struct {
		FilterData
		alpha, expX, expY float64
	}

	MitchellFilter struct {
		FilterData
		B, C float64
	}

	LanczosSincFilter struct {
		FilterData
		tau float64
	}

	TriangleFilter struct {
		FilterData
	}
)

func (f *BoxFilter) Evaluate(x, y float64) float64 {
	return 1.0
}

func (f *BoxFilter) XWidth() float64    { return f.xWidth }
func (f *BoxFilter) YWidth() float64    { return f.yWidth }
func (f *BoxFilter) InvXWidth() float64 { return f.invXWidth }
func (f *BoxFilter) InvYWidth() float64 { return f.invYWidth }

func CreateBoxFilter(params *ParamSet) *BoxFilter {
	// "xwidth", 0.5
	// "ywidth", 0.5
	return nil
}

func (f *GaussianFilter) Evaluate(x, y float64) float64 {
	return f.gaussian(x, f.expX) * f.gaussian(y, f.expY)
}

func (f *GaussianFilter) XWidth() float64    { return f.xWidth }
func (f *GaussianFilter) YWidth() float64    { return f.yWidth }
func (f *GaussianFilter) InvXWidth() float64 { return f.invXWidth }
func (f *GaussianFilter) InvYWidth() float64 { return f.invYWidth }

func (f *GaussianFilter) gaussian(d, expv float64) float64 {
	return math.Max(0.0, math.Exp(-f.alpha*d*d)-expv)
}

func CreateGaussianFilter(params *ParamSet) *GaussianFilter {
	// "xwidth", 2.0
	// "ywidth", 2.0
	// "alpha", 2.0
	return nil
}

func (f *MitchellFilter) Evaluate(x, y float64) float64 {
	return f.mitchell1D(x*f.invXWidth) * f.mitchell1D(y*f.invXWidth)
}

func (f *MitchellFilter) XWidth() float64    { return f.xWidth }
func (f *MitchellFilter) YWidth() float64    { return f.yWidth }
func (f *MitchellFilter) InvXWidth() float64 { return f.invXWidth }
func (f *MitchellFilter) InvYWidth() float64 { return f.invYWidth }

func (f *MitchellFilter) mitchell1D(x float64) float64 {
	x = math.Abs(2.0 * x)
	if x > 1.0 {
		return ((-f.B-6.0*f.C)*x*x*x + (6.0*f.B+30.0*f.C)*x*x +
			(-12.0*f.B-48.0*f.C)*x + (8.0*f.B + 24.0*f.C)) * (1.0 / 6.0)
	} else {
		return ((12.0-9.0*f.B-6.0*f.C)*x*x*x +
			(-18.0+12.0*f.B+6.0*f.C)*x*x +
			(6.0 - 2.0*f.B)) * (1.0 / 6.0)
	}
}

func CreateMitchellFilter(params *ParamSet) *MitchellFilter {
	// "xwidth", 2.0
	// "ywidth", 2.0
	// "B", 1/3
	// "C", 1/3
	return nil
}

func (f *LanczosSincFilter) Evaluate(x, y float64) float64 {
	return f.sinc1D(x*f.invXWidth) * f.sinc1D(y*f.invYWidth)
}

func (f *LanczosSincFilter) XWidth() float64    { return f.xWidth }
func (f *LanczosSincFilter) YWidth() float64    { return f.yWidth }
func (f *LanczosSincFilter) InvXWidth() float64 { return f.invXWidth }
func (f *LanczosSincFilter) InvYWidth() float64 { return f.invYWidth }

func (f *LanczosSincFilter) sinc1D(x float64) float64 {
	x = math.Abs(x)
	if x < 1.0e-5 {
		return 1.0
	}
	if x > 1.0 {
		return 0.0
	}
	x *= math.Pi
	sinc := math.Sin(x) / x
	lanczos := math.Sin(x*f.tau) / (x * f.tau)
	return sinc * lanczos
}

func CreateLanczosSincFilter(params *ParamSet) *LanczosSincFilter {
	// "xwidth", 4.0
	// "ywidth", 4.0
	// "tau", 3.0
	return nil
}

func (f *TriangleFilter) Evaluate(x, y float64) float64 {
	return math.Max(0.0, f.xWidth-math.Abs(x)) * math.Max(0.0, f.yWidth-math.Abs(y))
}

func (f *TriangleFilter) XWidth() float64    { return f.xWidth }
func (f *TriangleFilter) YWidth() float64    { return f.yWidth }
func (f *TriangleFilter) InvXWidth() float64 { return f.invXWidth }
func (f *TriangleFilter) InvYWidth() float64 { return f.invYWidth }

func CreateTriangleFilter(params *ParamSet) *TriangleFilter {
	// "xwidth", 2.0
	// "ywidth", 2.0
	return nil
}
