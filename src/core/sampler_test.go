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
	"testing"
)

func TestSampler(t *testing.T) {

	randsampler := NewRandomSampler(0, 16, 0, 16, 1, 0.0, 1.0)
	surf := &WhittedIntegrator{4}
	vol := NewEmissionIntegrator(1.0)

	sample := NewSample(randsampler, surf, vol, nil)
	fmt.Printf("Sample: %v\n", sample)

	for i, n := range sample.n1D {
		fmt.Printf("\tIndex: %d  Num: %d\n", i, n)
	}

	fmt.Printf("VolIntr: %v\n", vol)

	// Fetch a sample for the volume
	val := sample.oneD[vol.scatterSampleOffset][0]
	fmt.Printf("VolSample: %f\n", val)

	rng := NewRNG(13)
	for i := 0; i < 1000; i++ {
		v := rng.RandomFloat()
		if v < 0.0 || v > 1.0 {
			t.Fail()
		}
	}
}

func TestLDSampler(t *testing.T) {
	ldsampler := NewLDSampler(0, 16, 0, 16, 8, 0.0, 1.0)

	surf := &WhittedIntegrator{4}
	vol := NewEmissionIntegrator(1.0)

	sample := NewSample(ldsampler, surf, vol, nil)
	fmt.Printf("Sample: %v\n", sample)
	samples := make([]Sample, ldsampler.MaximumSampleCount(), ldsampler.MaximumSampleCount())
	for i := 0; i < len(samples); i++ {
		samples[i] = *sample
	}

	rng := NewRNG(13)
	ns := ldsampler.GetMoreSamples(&samples, rng)
	fmt.Printf("Ns: %d\n", ns)
	for i := 0; i < ns; i++ {
		fmt.Printf("Sample[%d]: %v\n", i, samples[i])
	}
}
