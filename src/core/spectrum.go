package core

import (
	"fmt"
	"io"
)

func XYZToRGB(xyz [3]float32) [3]float32 {
	rgb := [3]float32{0.0, 0.0, 0.0}
	rgb[0] = 3.240479*xyz[0] - 1.537150*xyz[1] - 0.498535*xyz[2]
	rgb[1] = -0.969256*xyz[0] + 1.875991*xyz[1] + 0.041556*xyz[2]
	rgb[2] = 0.055648*xyz[0] - 0.204043*xyz[1] + 1.057311*xyz[2]
	return rgb
}

func RGBToXYZ(rgb [3]float32) [3]float32 {
	xyz := [3]float32{0.0, 0.0, 0.0}
	xyz[0] = 0.412453*rgb[0] + 0.357580*rgb[1] + 0.180423*rgb[2]
	xyz[1] = 0.212671*rgb[0] + 0.715160*rgb[1] + 0.072169*rgb[2]
	xyz[2] = 0.019334*rgb[0] + 0.119193*rgb[1] + 0.950227*rgb[2]
	return xyz
}

type Spectrum struct {
	c [3]float32
}

type SampledSpectrum struct {
	c []float32
}

func (rgb *Spectrum) String() string {
	return fmt.Sprintf("rgb[%2.4f,%2.4f,%2.4f]", rgb.c[0], rgb.c[1], rgb.c[2])
}

func NewSpectrum1(v float32) *Spectrum {
	return &Spectrum{[3]float32{v, v, v}}
}
func NewSpectrum(v [3]float32) *Spectrum {
	return &Spectrum{[3]float32{v[0], v[1], v[2]}}
}
func NewSpectrumRGB(r, g, b float32) *Spectrum {
	return &Spectrum{[3]float32{r, g, b}}
}
func (rgb *Spectrum) Add(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{rgb.c[0] + s2.c[0], rgb.c[1] + s2.c[1], rgb.c[2] + s2.c[2]}}
}

func (rgb *Spectrum) Sub(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{rgb.c[0] - s2.c[0], rgb.c[1] - s2.c[1], rgb.c[2] - s2.c[2]}}
}

func (rgb *Spectrum) Divide(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{rgb.c[0] / s2.c[0], rgb.c[1] / s2.c[1], rgb.c[2] / s2.c[2]}}
}

func (rgb *Spectrum) Mult(s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{rgb.c[0] * s2.c[0], rgb.c[1] * s2.c[1], rgb.c[2] * s2.c[2]}}
}

func (rgb *Spectrum) Scale(s float32) *Spectrum {
	return &Spectrum{[3]float32{rgb.c[0] * s, rgb.c[1] * s, rgb.c[2] * s}}
}

func (rgb *Spectrum) InvScale(invs float32) *Spectrum {
	s := 1.0 / invs
	return &Spectrum{[3]float32{rgb.c[0] * s, rgb.c[1] * s, rgb.c[2] * s}}
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
	return &Spectrum{[3]float32{Sqrtf(s.c[0]), Sqrtf(s.c[1]), Sqrtf(s.c[2])}}
}

func PowSpectrum(s *Spectrum, e float32) *Spectrum {
	return &Spectrum{[3]float32{Powf(s.c[0], e), Powf(s.c[1], e), Powf(s.c[2], e)}}
}

func (rgb *Spectrum) Negate() *Spectrum {
	return &Spectrum{[3]float32{-rgb.c[0], -rgb.c[1], -rgb.c[2]}}
}

func ExpSpectrum(s *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{Expf(s.c[0]), Expf(s.c[1]), Expf(s.c[2])}}
}

func (rgb *Spectrum) Clamp(low, high float32) *Spectrum {
	return &Spectrum{[3]float32{Clampf(rgb.c[0], low, high), Clampf(rgb.c[1], low, high), Clampf(rgb.c[2], low, high)}}
}

func (rgb *Spectrum) HasNaNs() bool {
	for _, v := range rgb.c {
		if IsNaNf(v) {
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

func (rgb *Spectrum) ToXYZ() [3]float32 {
	return RGBToXYZ(rgb.c)
}

func SpectrumFromXYZ(xyz [3]float32) *Spectrum {
	return &Spectrum{c: XYZToRGB(xyz)}
}

func (rgb *Spectrum) Y() float32 {
	YWeight := []float32{0.212671, 0.715160, 0.072169}
	return YWeight[0]*rgb.c[0] + YWeight[1]*rgb.c[1] + YWeight[2]*rgb.c[2]
}

func LerpSpectrum(t float32, s1, s2 *Spectrum) *Spectrum {
	return &Spectrum{[3]float32{Lerpf(t, s1.c[0], s2.c[0]), Lerpf(t, s1.c[1], s2.c[1]), Lerpf(t, s1.c[2], s2.c[2])}}
}
