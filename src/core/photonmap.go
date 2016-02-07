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

type (

	PhotonIntegrator struct {
		nCausticPhotonsWanted, nIndirectPhotonsWanted, nLookup int
		maxDistSquared                                         float64
		maxSpecularDepth, maxPhotonDepth                       int
		finalGather                                            bool
		gatherSamples                                          int
		cosGatherAngle                                         float64

		// Declare sample parameters for light source sampling
		lightSampleOffsets                                []LightSampleOffsets
		bsdfSampleOffsets                                 []BSDFSampleOffsets
		bsdfGatherSampleOffsets, indirGatherSampleOffsets BSDFSampleOffsets
		nCausticPaths, nIndirectPaths                     int
		causticMap, indirectMap                           *KdTree // Photons
		radianceMap                                       *KdTree // RadiancePhoton
	}
	
	photon struct {
		p Point
		alpha Spectrum
		wi Vector
	}
	
	radiancePhoton struct {
		p Point
		n Normal
		Lo Spectrum
	}
	
	closePhoton struct {
		photon *photon
		distanceSquared float64
	}
)


func NewPhotonIntegrator(nCaustic, nIndirect, nUsed, maxSpecularDepth, maxPhotonDepth int, maxDist float64, finalGather bool, gatherSamples int, gatherAngle float64) *PhotonIntegrator {
	Unimplemented()
	return nil
}
func (integrator *PhotonIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (integrator *PhotonIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *PhotonIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	Unimplemented()
	return nil
}

