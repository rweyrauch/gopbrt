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
