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
	"testing"
	//"math"
)

/*
func AlmostEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}
func AlmostEqualRelative(a, b, tol float64) bool {
	diff := math.Abs(a-b)
	largest := math.Max(math.Abs(a), math.Abs(b))
	return diff <= largest * tol
}
*/

type FrDielInputs struct {
	cosi, cost float64
	etai       [3]float64
	etat       [3]float64
}
type FrCondInputs struct {
	cosi float64
	n    [3]float64
	k    [3]float64
}

var frDielInputVector []FrDielInputs = []FrDielInputs{
	{1.0, 1.0, [3]float64{0.5, 0.5, 0.5}, [3]float64{0.25, 0.25, 0.25}},
	{0.5, 0.717, [3]float64{1.0, 0.5, 0.25}, [3]float64{0.1, 0.2, 0.3}},
	{-0.2, 0.0, [3]float64{0.5, 0.85, 0.65}, [3]float64{0.75, 0.05, 0.5}},
}

var frDielOutputVector [][3]float64 = [][3]float64{
	[3]float64{0.111111, 0.111111, 0.111111},
	[3]float64{0.658748, 0.195645, 0.039038},
	[3]float64{1.000000, 1.000000, 1.000000},
}

var frCondInputVector []FrCondInputs = []FrCondInputs{
	{1.0, [3]float64{0.5, 0.5, 0.5}, [3]float64{0.25, 0.25, 0.25}},
	{0.717, [3]float64{1.0, 0.5, 0.25}, [3]float64{0.1, 0.2, 0.3}},
	{-0.2, [3]float64{0.5, 0.85, 0.65}, [3]float64{0.75, 0.05, 0.5}},
}

var frCondOutputVector [][3]float64 = [][3]float64{
	[3]float64{0.135135, 0.135135, 0.135135},
	[3]float64{0.029658, 0.144400, 0.400779},
	[3]float64{1.546754, 2.293469, 1.913613},
}

func TestFresnel(t *testing.T) {
	var options *Options
	PbrtInit(options)

	// Fresnel utility functions.
	for i, fr := range frDielInputVector {
		frd := FrDiel(fr.cosi, fr.cost, NewSpectrum(fr.etai), NewSpectrum(fr.etat))
		if !AlmostEqualRelative(frd.c[0], frDielOutputVector[i][0], tol) ||
			!AlmostEqualRelative(frd.c[1], frDielOutputVector[i][1], tol) ||
			!AlmostEqualRelative(frd.c[2], frDielOutputVector[i][2], tol) {
			t.Errorf("FrDiel function failed.")
			t.Logf("FrDiel In: %2.10f, %2.10f { %2.10f, %2.10f, %2.10f }, { %2.10f, %2.10f, %2.10f },\n",
				fr.cosi, fr.cost, fr.etai[0], fr.etai[1], fr.etai[2], fr.etat[0], fr.etat[1], fr.etat[2])
			t.Logf("Output: %f, %f, %f  Expected: %f, %f, %f", frd.c[0], frd.c[1], frd.c[2], frDielOutputVector[i][0], frDielOutputVector[i][1], frDielOutputVector[i][2])
		}
	}

	for i, fr := range frCondInputVector {
		frc := FrCond(fr.cosi, NewSpectrum(fr.n), NewSpectrum(fr.k))
		if !AlmostEqualRelative(frc.c[0], frCondOutputVector[i][0], tol) ||
			!AlmostEqualRelative(frc.c[1], frCondOutputVector[i][1], tol) ||
			!AlmostEqualRelative(frc.c[2], frCondOutputVector[i][2], tol) {
			t.Errorf("FrCond function failed.")
			t.Logf("FrCond In: %2.10f { %2.10f, %2.10f, %2.10f }, { %2.10f, %2.10f, %2.10f },\n",
				fr.cosi, fr.n[0], fr.n[1], fr.n[2], fr.k[0], fr.k[1], fr.k[2])
			t.Logf("Output: %f, %f, %f  Expected: %f, %f, %f", frc.c[0], frc.c[1], frc.c[2], frCondOutputVector[i][0], frCondOutputVector[i][1], frCondOutputVector[i][2])
		}
	}

	PbrtCleanup()
}

func TestSpecularReflection(t *testing.T) {
	var options *Options
	PbrtInit(options)

	u1 := 0.5
	u2 := 0.25
	wo := CreateVector(0.000000, -0.866025, 0.500000)
	wi := CreateVector(0.000000, 0.866025, 0.500000)

	r := NewSpectrum1(0.333)
	ei := 0.5
	et := 0.25
	f := &FresnelDielectric{ei, et}
	spec := &SpecularReflection{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)}, r, f}

	outF := [3]float64{0.000000, 0.000000, 0.000000}
	outPdf := 0.000000
	outSampF := [3]float64{0.666000, 0.666000, 0.666000}
	outSampPdf := 1.000000

	specF := spec.F(wo, wi)
	if !AlmostEqualRelative(specF.c[0], outF[0], tol) || !AlmostEqualRelative(specF.c[1], outF[1], tol) || !AlmostEqualRelative(specF.c[2], outF[2], tol) {
		t.Errorf("SpecularReflection.F() failed.")
	}
	specPdf := spec.Pdf(wo, wi)
	if !AlmostEqualRelative(specPdf, outPdf, tol) {
		t.Errorf("SpecularReflection.Pdf() failed.")
	}
	_, sampF, pdf := spec.Sample_f(wo, u1, u2)
	if !AlmostEqualRelative(sampF.c[0], outSampF[0], tol) || !AlmostEqualRelative(sampF.c[1], outSampF[1], tol) || !AlmostEqualRelative(sampF.c[2], outSampF[2], tol) ||
		!AlmostEqualRelative(pdf, outSampPdf, tol) {
		t.Errorf("SpecularReflection.Sample_f() failed.")
		t.Logf("Output: F{%f, %f, %f}, pdf: %f  Expected: F{%f, %f, %f}, pdf: %f", sampF.c[0], sampF.c[1], sampF.c[2], pdf, outSampF[0], outSampF[1], outSampF[2], outSampPdf)
	}

	ei = 1.0
	et = 1.0

	tx := NewSpectrum1(0.22)
	trans := &SpecularTransmission{BxDFData{BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)}, tx, ei, et, &FresnelDielectric{ei, et}}

	outTransF := [3]float64{0.000000, 0.000000, 0.000000}
	outTransPdf := 0.0
	outTransSampF := [3]float64{0.440000, 0.440000, 0.440000}
	outTransSampPdf := 1.0

	transF := trans.F(wo, wi)
	if !AlmostEqualRelative(transF.c[0], outTransF[0], tol) || !AlmostEqualRelative(transF.c[1], outTransF[1], tol) || !AlmostEqualRelative(transF.c[2], outTransF[2], tol) {
		t.Errorf("SpecularTransmission.F() failed.")
	}

	transPdf := trans.Pdf(wo, wi)
	if !AlmostEqualRelative(transPdf, outTransPdf, tol) {
		t.Errorf("SpecularTransmission.Pdf() failed.")
	}

	_, transSampF, pdf := trans.Sample_f(wo, u1, u2)
	if !AlmostEqualRelative(transSampF.c[0], outTransSampF[0], tol) || !AlmostEqualRelative(transSampF.c[1], outTransSampF[1], tol) || !AlmostEqualRelative(transSampF.c[2], outTransSampF[2], tol) ||
		!AlmostEqualRelative(pdf, outTransSampPdf, tol) {
		t.Errorf("SpecularTransmission.Sample_f() failed.")
		t.Logf("Output: F{%f, %f, %f}, pdf: %f  Expected: F{%f, %f, %f}, pdf: %f", transSampF.c[0], transSampF.c[1], transSampF.c[2], pdf, outTransSampF[0], outTransSampF[1], outTransSampF[2], outTransSampPdf)
	}

	PbrtCleanup()
}
