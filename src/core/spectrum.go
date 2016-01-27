package core

import (
	"fmt"
	"io"
	"math"
	"sort"
)

func XYZToRGB(xyz [3]float64) [3]float64 {
	rgb := [3]float64{0.0, 0.0, 0.0}
	rgb[0] = 3.240479*xyz[0] - 1.537150*xyz[1] - 0.498535*xyz[2]
	rgb[1] = -0.969256*xyz[0] + 1.875991*xyz[1] + 0.041556*xyz[2]
	rgb[2] = 0.055648*xyz[0] - 0.204043*xyz[1] + 1.057311*xyz[2]
	return rgb
}

func RGBToXYZ(rgb [3]float64) [3]float64 {
	xyz := [3]float64{0.0, 0.0, 0.0}
	xyz[0] = 0.412453*rgb[0] + 0.357580*rgb[1] + 0.180423*rgb[2]
	xyz[1] = 0.212671*rgb[0] + 0.715160*rgb[1] + 0.072169*rgb[2]
	xyz[2] = 0.019334*rgb[0] + 0.119193*rgb[1] + 0.950227*rgb[2]
	return xyz
}

type Spectrum struct {
	c [3]float64
}

type SampledSpectrum struct {
	c []float64
}

func (rgb *Spectrum) String() string {
	return fmt.Sprintf("rgb[%2.4f,%2.4f,%2.4f]", rgb.c[0], rgb.c[1], rgb.c[2])
}

func NewSpectrum1(v float64) *Spectrum {
	return &Spectrum{[3]float64{v, v, v}}
}
func NewSpectrum(v [3]float64) *Spectrum {
	return &Spectrum{[3]float64{v[0], v[1], v[2]}}
}
func NewSpectrumRGB(r, g, b float64) *Spectrum {
	return &Spectrum{[3]float64{r, g, b}}
}
func (rgb *Spectrum) Add(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{rgb.c[0] + s2.c[0], rgb.c[1] + s2.c[1], rgb.c[2] + s2.c[2]}}
}

func (rgb *Spectrum) Sub(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{rgb.c[0] - s2.c[0], rgb.c[1] - s2.c[1], rgb.c[2] - s2.c[2]}}
}

func (rgb *Spectrum) Divide(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{rgb.c[0] / s2.c[0], rgb.c[1] / s2.c[1], rgb.c[2] / s2.c[2]}}
}

func (rgb *Spectrum) Mult(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{rgb.c[0] * s2.c[0], rgb.c[1] * s2.c[1], rgb.c[2] * s2.c[2]}}
}

func (rgb *Spectrum) Scale(s float64) *Spectrum {
	return &Spectrum{[3]float64{rgb.c[0] * s, rgb.c[1] * s, rgb.c[2] * s}}
}

func (rgb *Spectrum) InvScale(invs float64) *Spectrum {
	s := 1.0 / invs
	return &Spectrum{[3]float64{rgb.c[0] * s, rgb.c[1] * s, rgb.c[2] * s}}
}

func (rgb *Spectrum) Equal(s2 *Spectrum) bool {
	for i, v := range rgb.c {
		if v != s2.c[i] {
			return false
		}
	}
	return true
}

func (rgb *Spectrum) NotEqual(s2 *Spectrum) bool {
	return !rgb.Equal(s2)
}

func (rgb *Spectrum) IsBlack() bool {
	for _, v := range rgb.c {
		if v != 0.0 {
			return false
		}
	}
	return true
}

func SqrtSpectrum(s *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{math.Sqrt(s.c[0]), math.Sqrt(s.c[1]), math.Sqrt(s.c[2])}}
}

func PowSpectrum(s *Spectrum, e float64) *Spectrum {
	return &Spectrum{[3]float64{math.Pow(s.c[0], e), math.Pow(s.c[1], e), math.Pow(s.c[2], e)}}
}

func (rgb *Spectrum) Negate() *Spectrum {
	return &Spectrum{[3]float64{-rgb.c[0], -rgb.c[1], -rgb.c[2]}}
}

func ExpSpectrum(s *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{math.Exp(s.c[0]), math.Exp(s.c[1]), math.Exp(s.c[2])}}
}

func (rgb *Spectrum) Clamp(low, high float64) *Spectrum {
	return &Spectrum{[3]float64{Clamp(rgb.c[0], low, high), Clamp(rgb.c[1], low, high), Clamp(rgb.c[2], low, high)}}
}

func (rgb *Spectrum) HasNaNs() bool {
	for _, v := range rgb.c {
		if math.IsNaN(v) {
			return true
		}
	}
	return false
}
func (rgb *Spectrum) Write(w io.Writer) bool {
	_, err := fmt.Fprintf(w, "%f %f %f", rgb.c[0], rgb.c[1], rgb.c[2])
	return err == nil
}
func (rgb *Spectrum) Read(r io.Reader) bool {
	_, err := fmt.Fscanf(r, "%f %f %f", rgb.c[0], rgb.c[1], rgb.c[2])
	return err == nil
}

func (rgb *Spectrum) ToXYZ() [3]float64 {
	return RGBToXYZ(rgb.c)
}

func SpectrumFromXYZ(xyz [3]float64) *Spectrum {
	return &Spectrum{c: XYZToRGB(xyz)}
}

func (rgb *Spectrum) Y() float64 {
	YWeight := []float64{0.212671, 0.715160, 0.072169}
	return YWeight[0]*rgb.c[0] + YWeight[1]*rgb.c[1] + YWeight[2]*rgb.c[2]
}

func SpectrumFromSampled(lambdas, vals []float64) *Spectrum {
	// Sort samples if unordered, use sorted for returned spectrum
	if !spectrumSamplesSorted(lambdas) {
		slambdas, svals := sortSpectrumSamples(lambdas, vals)
		return SpectrumFromSampled(slambdas, svals)
	}
	xyz := [3]float64{0, 0, 0}
	var yint float64 = 0.0
	for i := 0; i < nCIESamples; i++ {
		yint += CIE_Y[i]
		val := interpolateSpectrumSamples(lambdas, vals, CIE_lambda[i])
		xyz[0] += val * CIE_X[i]
		xyz[1] += val * CIE_Y[i]
		xyz[2] += val * CIE_Z[i]
	}
	xyz[0] /= yint
	xyz[1] /= yint
	xyz[2] /= yint
	return SpectrumFromXYZ(xyz)
}

func LerpSpectrum(t float64, s1, s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float64{Lerp(t, s1.c[0], s2.c[0]), Lerp(t, s1.c[1], s2.c[1]), Lerp(t, s1.c[2], s2.c[2])}}
}

func spectrumSamplesSorted(lambda []float64) bool {
	for i := 0; i < len(lambda)-1; i++ {
		if lambda[i] > lambda[i+1] {
			return false
		}
	}
	return true
}

type sampledSpectrumSorter struct {
	lambdas, vals []float64
}

func (s *sampledSpectrumSorter) Len() int {
	return len(s.lambdas)
}
func (s *sampledSpectrumSorter) Swap(i, j int) {
	s.lambdas[i], s.lambdas[j] = s.lambdas[j], s.lambdas[i]
	s.vals[i], s.vals[j] = s.vals[j], s.vals[i]
}
func (s *sampledSpectrumSorter) Less(i, j int) bool {
	return s.lambdas[i] < s.lambdas[j]
}

func sortSpectrumSamples(lambdas, vals []float64) (slambdas, svals []float64) {
	slambdas = make([]float64, len(lambdas), len(lambdas))
	svals = make([]float64, len(vals), len(vals))
	copy(slambdas, lambdas)
	copy(svals, vals)

	spectrumSorter := &sampledSpectrumSorter{slambdas, svals}
	sort.Sort(spectrumSorter)

	return slambdas, svals
}

func interpolateSpectrumSamples(lambdas, vals []float64, l float64) float64 {
	Assert(spectrumSamplesSorted(lambdas))
	if l <= lambdas[0] {
		return vals[0]
	}
	if l >= lambdas[len(lambdas)-1] {
		return vals[len(vals)-1]
	}
	for i := 0; i < len(lambdas)-1; i++ {
		if l >= lambdas[i] && l <= lambdas[i+1] {
			t := (l - lambdas[i]) / (lambdas[i+1] - lambdas[i])
			return Lerp(t, vals[i], vals[i+1])
		}
	}
	Severe("Fatal logic error in InterpolateSpectrumSamples()")
	return 0.0
}

func Blackbody(wl []float64, temp float64) []float64 {
	vals := make([]float64, len(wl), len(wl))
    if temp <= 0 {
        for i, _ := range vals { vals[i] = 0.0 }
        return vals
    }
    C2 := 1.4388e7
    norm := math.Pow(555.0, 5.0) * (math.Exp(C2/(555.0*temp)) - 1.0)
    for i, w := range wl {
        vals[i] = norm/(math.Pow(w, 5.0)*(math.Exp(C2/(w*temp))-1.0))
	}
	return vals
}
