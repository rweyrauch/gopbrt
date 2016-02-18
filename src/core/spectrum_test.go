/*
	gopbrt

	Port of pbrt v2.0.0 by Matt Pharr and Greg Humphreys to the go language.
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

	The MIT License (MIT)
	Copyright (c) 2016 Rick Weyrauch

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
package core

import (
	"fmt"
	"math"
	"testing"
)

type rgbOutputData struct {
	rgb [3]float64
	xyz [3]float64
	y   float64
}

var rgbInputVector [][3]float64 = [][3]float64{
	{1.0, 1.0, 1.0},
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{0.5, 0.5, 0.0},
	{0.0, 0.75, 1.0},
}

// Output data captured using pbrt v2 final.
var rgbOutputVector []rgbOutputData = []rgbOutputData{
	{[3]float64{1.0, 1.0, 1.0}, [3]float64{0.9504560232, 1.0000000000, 1.0887540579}, 1.0},
	{[3]float64{0.0, 0.0, 0.0}, [3]float64{0.0000000000, 0.0000000000, 0.0000000000}, 0.0},
	{[3]float64{1.0, 0.0, 0.0}, [3]float64{0.4124529958, 0.2126709968, 0.0193339996}, 0.2126709968},
	{[3]float64{0.5, 0.5, 0.0}, [3]float64{0.3850165009, 0.4639154971, 0.0692635030}, 0.4639154971},
	{[3]float64{0.0, 0.75, 1.0}, [3]float64{0.4486080408, 0.6085390449, 1.0396218300}, 0.6085390449},
}

const (
	tol = 0.0001
)

var blackBodySpectrumInputVector [][2]float64 = [][2]float64{
	{20000.0, 0.8},
	{5500.0, 0.5},
	{12005.0, 1.0},
	{0.0, 1.0},
	{0.0, 0.0},
	{400.0, 0.5},
}

// Output data captured using pbrt v2 final.
var blackBodySpectrumOutputVector [][3]float64 = [][3]float64{
	{0.602281, 0.796446, 1.485405},
	{0.564041, 0.478611, 0.425158},
	{0.823771, 0.991291, 1.562498},
	{0.000000, 0.000000, 0.000000},
	{0.000000, 0.000000, 0.000000},
	{4053.423096, -420.130829, -29.750914},
}

func AlmostEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}
func AlmostEqualRelative(a, b, tol float64) bool {
	diff := math.Abs(a - b)
	largest := math.Max(math.Abs(a), math.Abs(b))
	return diff <= largest*tol
}

func TestSpectrumRGB(t *testing.T) {
	var options *Options
	PbrtInit(options)

	for i, rgbin := range rgbInputVector {
		rgb := NewSpectrum(rgbin)
		if !AlmostEqualRelative(rgb.c[0], rgbOutputVector[i].rgb[0], tol) ||
			!AlmostEqualRelative(rgb.c[1], rgbOutputVector[i].rgb[1], tol) ||
			!AlmostEqualRelative(rgb.c[2], rgbOutputVector[i].rgb[2], tol) {
			t.Errorf("Conversion to RGB failed.")
			t.Logf("RGBin: %2.10f,%2.10f,%2.10f RGBout: %2.10f,%2.10f,%2.10f",
				rgbin[0], rgbin[1], rgbin[2], rgb.c[0], rgb.c[1], rgb.c[2])
		}
		xyz := rgb.ToXYZ()
		y := rgb.Y()
		if !AlmostEqualRelative(xyz[0], rgbOutputVector[i].xyz[0], tol) ||
			!AlmostEqualRelative(xyz[0], rgbOutputVector[i].xyz[0], tol) ||
			!AlmostEqualRelative(xyz[0], rgbOutputVector[i].xyz[0], tol) {
			t.Errorf("XYZ conversion failed.")
			t.Logf("RGBin: %2.10f,%2.10f,%2.10f XYZout: %2.10f,%2.10f,%2.10f", rgbin[0], rgbin[1], rgbin[2], xyz[0], xyz[1], xyz[2])
		}
		if !AlmostEqualRelative(y, rgbOutputVector[i].y, tol) {
			t.Errorf("Lumanance conversion failed.")
			t.Logf("RGBin: %2.10f,%2.10f,%2.10f Y: %2.10f", rgbin[0], rgbin[1], rgbin[2], y)
		}
	}
	PbrtCleanup()
}

func TestSpectrumBlackbody(t *testing.T) {
	var options *Options
	PbrtInit(options)

	for i, bbin := range blackBodySpectrumInputVector {
		v := Blackbody(CIE_lambda, bbin[0])
		rgb := *SpectrumFromSampled(CIE_lambda, v).Scale(bbin[1])
		if !AlmostEqualRelative(rgb.c[0], blackBodySpectrumOutputVector[i][0], tol) ||
			!AlmostEqualRelative(rgb.c[1], blackBodySpectrumOutputVector[i][1], tol) ||
			!AlmostEqualRelative(rgb.c[2], blackBodySpectrumOutputVector[i][2], tol) {
			t.Errorf("Blackbody conversion failed.")
			t.Logf("Blackbody: T: %2.10f Sc: %2.10f  RGBout: %2.10f,%2.10f,%2.10f", bbin[0], bbin[1], rgb.c[0], rgb.c[1], rgb.c[2])
		}
	}
	PbrtCleanup()
}

var Au_eta_spd = []float64{
	298.757050, 1.795000,
	302.400421, 1.812000,
	306.133759, 1.822625,
	309.960449, 1.830000,
	313.884003, 1.837125,
	317.908142, 1.840000,
	322.036835, 1.834250,
	326.274139, 1.824000,
	330.624481, 1.812000,
	335.092377, 1.798000,
	339.682678, 1.782000,
	344.400482, 1.766000,
	349.251221, 1.752500,
	354.240509, 1.740000,
	359.374420, 1.727625,
	364.659332, 1.716000,
	370.102020, 1.705875,
	375.709625, 1.696000,
	381.489777, 1.684750,
	387.450562, 1.674000,
	393.600555, 1.666000,
	399.948975, 1.658000,
	406.505493, 1.647250,
	413.280579, 1.636000,
	420.285339, 1.628000,
	427.531647, 1.616000,
	435.032196, 1.596250,
	442.800629, 1.562000,
	450.851562, 1.502125,
	459.200653, 1.426000,
	467.864838, 1.345875,
	476.862213, 1.242000,
	486.212463, 1.086750,
	495.936707, 0.916000,
	506.057861, 0.754500,
	516.600769, 0.608000,
	527.592224, 0.491750,
	539.061646, 0.402000,
	551.040771, 0.345500,
	563.564453, 0.306000,
	576.670593, 0.267625,
	590.400818, 0.236000,
	604.800842, 0.212375,
	619.920898, 0.194000,
	635.816284, 0.177750,
	652.548279, 0.166000,
	670.184753, 0.161000,
	688.800964, 0.160000,
	708.481018, 0.160875,
	729.318665, 0.164000,
	751.419250, 0.169500,
	774.901123, 0.176000,
	799.897949, 0.181375,
	826.561157, 0.188000,
	855.063293, 0.198125,
	885.601257, 0.210000,
}

var Au_eta_rgb = []float64{0.1431240439, 0.3749396205, 1.4425071478}

func TestSpectrumFromSampled(t *testing.T) {
	var options *Options
	PbrtInit(options)

	wls := make([]float64, 0, len(Au_eta_spd)/2)
	v := make([]float64, 0, len(Au_eta_spd)/2)
	for j := 0; j < len(Au_eta_spd)/2; j++ {
		wls = append(wls, Au_eta_spd[2*j])
		v = append(v, Au_eta_spd[2*j+1])
	}
	value := *SpectrumFromSampled(wls, v)

	fmt.Printf("Au Eta: %v\n", value)

	if !AlmostEqualRelative(value.c[0], Au_eta_rgb[0], tol) ||
		!AlmostEqualRelative(value.c[1], Au_eta_rgb[1], tol) ||
		!AlmostEqualRelative(value.c[2], Au_eta_rgb[2], tol) {
		t.Errorf("Sampled spectrum conversion failed.")
	}

	PbrtCleanup()
}
