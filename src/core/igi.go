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
	"math"
)

type (
	IGIIntegrator struct {
		// Declare sample parameters for light source sampling
		lightSampleOffsets      []LightSampleOffsets
		bsdfSampleOffsets       []BSDFSampleOffsets
		nLightPaths, nLightSets int
		gLimit                  float64
		nGatherSamples          int
		rrThreshold             float64
		maxSpecularDepth        int
		vlSetOffset             int
		gatherSampleOffset      BSDFSampleOffsets
		virtualLights           [][]VirtualLight
	}

	VirtualLight struct {
		p           Point
		n           Normal
		pathContrib Spectrum
		rayEpsilon  float64
	}
)

func NewIGIIntegrator(nLightPaths, nLightSets int, rrThresh float64, maxDepth int, glimit float64, gatherSamples int) *IGIIntegrator {
	integrator := new(IGIIntegrator)
	integrator.nLightPaths = int(RoundUpPow2(uint32(nLightPaths)))
	integrator.nLightSets = int(RoundUpPow2(uint32(nLightSets)))
	integrator.rrThreshold = rrThresh
	integrator.maxSpecularDepth = maxDepth
	integrator.virtualLights = make([][]VirtualLight, integrator.nLightSets, integrator.nLightSets)
	integrator.gLimit = glimit
	integrator.nGatherSamples = gatherSamples
	integrator.lightSampleOffsets = nil
	integrator.bsdfSampleOffsets = nil

	return integrator
}

func (integrator *IGIIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer) {
	if len(scene.lights) == 0 {
		return
	}
	var arena *MemoryArena
	rng := NewRNG(33)
	// Compute samples for emitted rays from lights
	lightNum := make([]float64, integrator.nLightPaths*integrator.nLightSets, integrator.nLightPaths*integrator.nLightSets)
	lightSampPos := make([]float64, 2*integrator.nLightPaths*integrator.nLightSets, 2*integrator.nLightPaths*integrator.nLightSets)
	lightSampComp := make([]float64, integrator.nLightPaths*integrator.nLightSets, integrator.nLightPaths*integrator.nLightSets)
	lightSampDir := make([]float64, 2*integrator.nLightPaths*integrator.nLightSets, 2*integrator.nLightPaths*integrator.nLightSets)
	LDShuffleScrambled1D(integrator.nLightPaths, integrator.nLightSets, &lightNum, rng)
	LDShuffleScrambled2D(integrator.nLightPaths, integrator.nLightSets, &lightSampPos, rng)
	LDShuffleScrambled1D(integrator.nLightPaths, integrator.nLightSets, &lightSampComp, rng)
	LDShuffleScrambled2D(integrator.nLightPaths, integrator.nLightSets, &lightSampDir, rng)

	// Precompute information for light sampling densities
	lightDistribution := ComputeLightSamplingCDF(scene)
	for s := 0; s < integrator.nLightSets; s++ {
		for i := 0; i < integrator.nLightPaths; i++ {
			// Follow path _i_ from light to create virtual lights
			sampOffset := s*integrator.nLightPaths + i

			// Choose light source to trace virtual light path from
			ln, lightPdf := lightDistribution.SampleDiscrete(lightNum[sampOffset])
			light := scene.lights[ln]

			// Sample ray leaving light source for virtual light path
			var ls LightSample
			ls.uPos[0] = lightSampPos[2*sampOffset]
			ls.uPos[1] = lightSampPos[2*sampOffset+1]
			ls.uComponent = lightSampComp[sampOffset]

			alpha, rr, Nl, pdf := light.Sample_L2(scene, &ls, lightSampDir[2*sampOffset],
				lightSampDir[2*sampOffset+1], camera.ShutterOpen())
			ray := CreateRayDifferentialFromRay(rr)
			if pdf == 0.0 || alpha.IsBlack() {
				continue
			}
			alpha = alpha.Scale(AbsDotNormalVector(Nl, &ray.Dir) / (pdf * lightPdf))
			hit, isect := scene.Intersect(ray)
			for hit && !alpha.IsBlack() {
				// Create virtual light and sample new ray for path
				alpha = alpha.Mult(renderer.Transmittance(scene, ray, nil, rng, arena))
				wo := ray.Dir.Negate()
				bsdf := isect.GetBSDF(ray, arena)

				// Create virtual light at ray intersection point
				contrib := alpha.Mult(bsdf.rho2(wo, rng, BSDF_ALL, 6).InvScale(math.Pi))
				integrator.virtualLights[s] = append(integrator.virtualLights[s], VirtualLight{*isect.dg.p, *isect.dg.nn, *contrib, isect.rayEpsilon})

				// Sample new ray direction and update weight for virtual light path
				bsdfSample := CreateRandomBSDFSample(rng)
				fr, wi, pdf, _ := bsdf.Sample_f(wo, bsdfSample, BSDF_ALL)
				if fr.IsBlack() || pdf == 0.0 {
					break
				}
				contribScale := fr.Scale(AbsDotVectorNormal(wi, bsdf.dgShading.nn) / pdf)

				// Possibly terminate virtual light path with Russian roulette
				rrProb := math.Min(1.0, contribScale.Y())
				if rng.RandomFloat() > rrProb {
					break
				}
				alpha = alpha.Mult(contribScale.InvScale(rrProb))
				ray = CreateChildRayDifferential(isect.dg.p, wi, CreateRayFromRayDifferential(ray), isect.rayEpsilon, INFINITY)
				hit, isect = scene.Intersect(ray)
			}
			//arena.FreeAll();
		}
	}
}

func (integrator *IGIIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
	// Allocate and request samples for sampling all lights
	nLights := len(scene.lights)
	integrator.lightSampleOffsets = make([]LightSampleOffsets, nLights, nLights)
	integrator.bsdfSampleOffsets = make([]BSDFSampleOffsets, nLights, nLights)
	for i := 0; i < nLights; i++ {
		light := scene.lights[i]
		nSamples := light.NumSamples()
		if sampler != nil {
			nSamples = sampler.RoundSize(nSamples)
		}
		integrator.lightSampleOffsets[i] = *CreateLightSampleOffsets(nSamples, sample)
		integrator.bsdfSampleOffsets[i] = *CreateBSDFSampleOffsets(nSamples, sample)
	}
	integrator.vlSetOffset = sample.Add1D(1)

	if sampler != nil {
		integrator.nGatherSamples = sampler.RoundSize(integrator.nGatherSamples)
	}
	integrator.gatherSampleOffset = *CreateBSDFSampleOffsets(integrator.nGatherSamples, sample)
}

func (integrator *IGIIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	wo := ray.Dir.Negate()
	// Compute emitted light if ray hit an area light source
	L = L.Add(isect.Le(wo))

	// Evaluate BSDF at hit point
	bsdf := isect.GetBSDF(ray, arena)
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n,
		wo, isect.rayEpsilon, ray.Time, bsdf, sample, rng,
		integrator.lightSampleOffsets, integrator.bsdfSampleOffsets))
	// Compute indirect illumination with virtual lights
	lSet := Mini(int(sample.oneD[integrator.vlSetOffset][0]*float64(integrator.nLightSets)), integrator.nLightSets-1)
	for i := 0; i < len(integrator.virtualLights[lSet]); i++ {
		vl := &integrator.virtualLights[lSet][i]
		// Compute virtual light's tentative contribution _Llight_
		d2 := DistanceSquaredPoint(p, &vl.p)
		wi := NormalizeVector(vl.p.Sub(p))
		G := AbsDotVectorNormal(wi, n) * AbsDotVectorNormal(wi, &vl.n) / d2
		G = math.Min(G, integrator.gLimit)
		f := bsdf.f(wo, wi, BSDF_ALL)
		if G == 0.0 || f.IsBlack() {
			continue
		}
		Llight := f.Scale(G).Mult(vl.pathContrib.InvScale(float64(integrator.nLightPaths)))
		connectRay := CreateChildRayDifferential(p, wi, CreateRayFromRayDifferential(ray), isect.rayEpsilon,
			math.Sqrt(d2)*(1.0-vl.rayEpsilon))
		Llight = Llight.Mult(renderer.Transmittance(scene, connectRay, nil, rng, arena))

		// Possibly skip virtual light shadow ray with Russian roulette
		if Llight.Y() < integrator.rrThreshold {
			continueProbability := 0.1
			if rng.RandomFloat() > continueProbability {
				continue
			}
			Llight = Llight.InvScale(continueProbability)
		}

		// Add contribution from _VirtualLight_ _vl_
		if !scene.IntersectP(CreateRayFromRayDifferential(connectRay)) {
			L = L.Add(Llight)
		}
	}
	if ray.Depth < integrator.maxSpecularDepth {
		// Do bias compensation for bounding geometry term
		nSamples := 1
		if ray.Depth == 0 {
			nSamples = integrator.nGatherSamples
		}
		for i := 0; i < nSamples; i++ {
			var bsdfSample *BSDFSample
			if ray.Depth == 0 {
				bsdfSample = CreateBSDFSample(sample, &integrator.gatherSampleOffset, i)
			} else {
				bsdfSample = CreateRandomBSDFSample(rng)
			}
			f, wi, pdf, _ := bsdf.Sample_f(wo, bsdfSample, BxDFType(BSDF_ALL^BSDF_SPECULAR))
			if !f.IsBlack() && pdf > 0.0 {
				// Trace ray for bias compensation gather sample
				maxDist := math.Sqrt(AbsDotVectorNormal(wi, n) / integrator.gLimit)
				gatherRay := CreateChildRayDifferential(p, wi, CreateRayFromRayDifferential(ray), isect.rayEpsilon, maxDist)
				Li, gatherIsect, _ := renderer.Li(scene, gatherRay, sample, rng, arena)
				if Li.IsBlack() {
					continue
				}

				// Add bias compensation ray contribution to radiance sum
				Ggather := AbsDotVectorNormal(wi, n) * AbsDotVectorNormal(wi.Negate(), gatherIsect.dg.nn) /
					DistanceSquaredPoint(p, gatherIsect.dg.p)
				if Ggather-integrator.gLimit > 0.0 && !math.IsInf(Ggather, 0) {
					gs := (Ggather - integrator.gLimit) / Ggather
					L = L.Add(f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) * gs / (float64(nSamples) * pdf))))
				}
			}
		}
	}
	if ray.Depth+1 < integrator.maxSpecularDepth {
		// Trace rays for specular reflection and refraction
		L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
		L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
	}
	return L
}

func CreateIGISurfaceIntegrator(params *ParamSet) *IGIIntegrator {
	nLightPaths := params.FindIntParam("nlights", 64)
	if options.FastRender {
		nLightPaths = Maxi(1, nLightPaths/2)
	} else if options.QuickRender {
		nLightPaths = Maxi(1, nLightPaths/4)
	}
	nLightSets := params.FindIntParam("nsets", 4)
	rrThresh := params.FindFloatParam("rrthreshold", 0.0001)
	maxDepth := params.FindIntParam("maxdepth", 5)
	glimit := params.FindFloatParam("glimit", 10.0)
	gatherSamples := params.FindIntParam("gathersamples", 16)
	return NewIGIIntegrator(nLightPaths, nLightSets, rrThresh, maxDepth, glimit, gatherSamples)
}
