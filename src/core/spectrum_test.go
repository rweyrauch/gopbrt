package core

import (
	"testing"
	"math"	
)

type rgbOutputData struct {
	rgb [3]float64
	xyz [3]float64
	y float64
}

var rgbInputVector [][3]float64 = [][3]float64{
	{ 1.0, 1.0, 1.0 },
	{ 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0 },
	{ 0.5, 0.5, 0.0 },
	{ 0.0, 0.75, 1.0 },
}
// Output data captured using pbrt v2 final.
var rgbOutputVector []rgbOutputData = []rgbOutputData{
	{[3]float64{ 1.0, 1.0, 1.0 }, [3]float64{ 0.9504560232,1.0000000000,1.0887540579 }, 1.0},
	{[3]float64{ 0.0, 0.0, 0.0 }, [3]float64{ 0.0000000000,0.0000000000,0.0000000000 }, 0.0},
	{[3]float64{ 1.0, 0.0, 0.0 }, [3]float64{ 0.4124529958,0.2126709968,0.0193339996 }, 0.2126709968},
	{[3]float64{ 0.5, 0.5, 0.0 }, [3]float64{ 0.3850165009,0.4639154971,0.0692635030 }, 0.4639154971},
	{[3]float64{ 0.0, 0.75, 1.0 }, [3]float64{ 0.4486080408,0.6085390449,1.0396218300 }, 0.6085390449},
}

const (
	tol = 0.00001
)

var blackBodySpectrumInputVector [][2]float64 = [][2]float64{
	{ 20000.0, 0.8 },
	{ 5500.0, 0.5 },
	{ 12005.0, 1.0 },
	{ 0.0, 1.0 },
	{ 0.0, 0.0 },
	{ 400.0, 0.5 },
}

// Output data captured using pbrt v2 final.
var blackBodySpectrumOutputVector [][3]float64 = [][3]float64{
	{ 0.602281, 0.796446, 1.485405 },
	{ 0.564041, 0.478611, 0.425158 },
	{ 0.823771, 0.991291, 1.562498 },
	{ 0.000000, 0.000000, 0.000000 },
	{ 0.000000, 0.000000, 0.000000 },
	{ 4053.423096, -420.130829, -29.750914 },
}

func AlmostEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}
func AlmostEqualRelative(a, b, tol float64) bool {
	diff := math.Abs(a-b)
	largest := math.Max(math.Abs(a), math.Abs(b))
	return diff <= largest * tol
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
