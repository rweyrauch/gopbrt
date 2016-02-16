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
	IrradianceCacheIntegrator struct {
		minSamplePixelSpacing, maxSamplePixelSpacing float64
		minWeight, cosMaxSampleAngleDifference       float64
		nSamples, maxSpecularDepth, maxIndirectDepth int

		// Declare sample parameters for light source sampling
		lightSampleOffsets []LightSampleOffsets
		bsdfSampleOffsets  []BSDFSampleOffsets
		octree             *Octree // IrradianceSample
	}

	irradianceSample struct {
		E       Spectrum
		n       Normal
		p       Point
		wAvg    Vector
		maxDist float64
	}

	irradiancePrimeTask struct {
		scene             *Scene
		camera            Camera
		renderer          Renderer
		sampler           Sampler
		origSample        *Sample
		irradianceCache   *IrradianceCacheIntegrator
		progress          *ProgressReporter
		taskNum, numTasks int
	}

	irradProcess struct {
		p                                             Point
		n                                             Normal
		minWeight, cosMaxSampleAngleDifference, sumWt float64
		nFound                                        int
		E                                             Spectrum
		wAvg                                          Vector
	}
)

func NewIrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing, maxAngle float64,
	maxSpecularDepth, maxIndirectDepth, nSamples int) *IrradianceCacheIntegrator {
	Unimplemented()
	return nil
}

func (integrator *IrradianceCacheIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer) {
	wb := scene.WorldBound()
	delta := (wb.PMax.Sub(&wb.PMin)).Scale(0.01)
	wb.PMin = *wb.PMin.SubVector(delta)
	wb.PMax = *wb.PMax.Add(delta)
	integrator.octree = NewOctree(wb, DEFAULT_OCTREE_MAX_DEPTH) // IrradianceSample
	// Prime irradiance cache
	integrator.minWeight *= 1.50
	xstart, xend, ystart, yend := camera.Film().GetSampleExtent()
	sampler := NewHaltonSampler(xstart, xend, ystart, yend, 1,
		camera.ShutterOpen(), camera.ShutterClose())
	sample := NewSample(sampler, integrator, nil, scene)
	nTasks := 64
	progress := NewProgressReporter(nTasks, "Priming irradiance cache", TerminalWidth())
	for i := 0; i < nTasks; i++ {
		task := newIrradiancePrimeTask(scene, renderer, camera,
			sampler, sample, integrator,
			progress, i, nTasks)
		task.run()
	}
	progress.Done()
	integrator.minWeight /= 1.5
}

func (integrator *IrradianceCacheIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
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
}

func (integrator *IrradianceCacheIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	// Evaluate BSDF at hit point
	bsdf := isect.GetBSDF(ray, arena)
	wo := ray.Dir.Negate()
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	L = L.Add(isect.Le(wo))
	// Compute direct lighting for irradiance cache
	L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n, wo,
		isect.rayEpsilon, ray.Time, bsdf, sample, rng,
		integrator.lightSampleOffsets, integrator.bsdfSampleOffsets))

	// Compute indirect lighting for irradiance cache
	if ray.Depth+1 < integrator.maxSpecularDepth {
		// Trace rays for specular reflection and refraction
		L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
		L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
	}

	// Estimate indirect lighting with irradiance cache
	ng := isect.dg.nn
	ng = FaceforwardNormalVector(ng, wo)

	// Compute pixel spacing in world space at intersection point
	pixelSpacing := math.Sqrt(CrossVector(isect.dg.dpdx, isect.dg.dpdy).Length())
	flags := BSDF_REFLECTION | BSDF_DIFFUSE | BSDF_GLOSSY
	L = L.Add(integrator.indirectLo(p, ng, pixelSpacing, wo, isect.rayEpsilon,
		bsdf, flags, rng, scene, renderer, arena))
	flags = BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY
	L = L.Add(integrator.indirectLo(p, ng.Negate(), pixelSpacing, wo, isect.rayEpsilon,
		bsdf, flags, rng, scene, renderer, arena))
	return L
}

func (integrator *IrradianceCacheIntegrator) indirectLo(p *Point, ng *Normal, pixelSpacing float64, wo *Vector,
	rayEpsilon float64, bsdf *BSDF, flags BxDFType, rng *RNG, scene *Scene, renderer Renderer, arena *MemoryArena) *Spectrum {
	if bsdf.NumComponentsMatching(flags) == 0 {
		return NewSpectrum1(0.0)
	}
	var E *Spectrum
	var wi *Vector
	var ok bool
	// Get irradiance _E_ and average incident direction _wi_ at point _p_
	if ok, E, wi = integrator.interpolateE(scene, p, ng); !ok {
		// Compute irradiance at current point
		//PBRT_IRRADIANCE_CACHE_STARTED_COMPUTING_IRRADIANCE(const_cast<Point *>(&p), const_cast<Normal *>(&ng));
		scramble := [2]uint32{rng.RandomUInt(), rng.RandomUInt()}
		minHitDistance := INFINITY
		wAvg := CreateVector(0, 0, 0)
		LiSum := NewSpectrum1(0.0)
		for i := 0; i < integrator.nSamples; i++ {
			// Sample direction for irradiance estimate ray
			u0, u1 := Sample02(uint32(i), scramble)
			w := CosineSampleHemisphere(u0, u1)
			r := CreateRayDifferential(p, bsdf.LocalToWorld(w), rayEpsilon, INFINITY, 0.0, 0)
			r.Dir = *FaceforwardVectorNormal(&r.Dir, ng)

			// Trace ray to sample radiance for irradiance estimate
			//PBRT_IRRADIANCE_CACHE_STARTED_RAY(&r);
			L := integrator.pathL(r, scene, renderer, rng, arena)
			LiSum = LiSum.Add(L)
			wAvg = wAvg.Add(r.Dir.Scale(L.Y()))
			minHitDistance = math.Min(minHitDistance, r.Maxt)
			//PBRT_IRRADIANCE_CACHE_FINISHED_RAY(&r, r.Maxt, &L)
		}
		E = LiSum.Scale(math.Pi / float64(integrator.nSamples))
		//PBRT_IRRADIANCE_CACHE_FINISHED_COMPUTING_IRRADIANCE(const_cast<Point *>(&p), const_cast<Normal *>(&ng));

		// Add computed irradiance value to cache

		// Compute irradiance sample's contribution extent and bounding box
		maxDist := integrator.maxSamplePixelSpacing * pixelSpacing
		minDist := integrator.minSamplePixelSpacing * pixelSpacing
		contribExtent := Clamp(minHitDistance/2.0, minDist, maxDist)
		sampleExtent := CreateBBoxFromPoint(p)
		sampleExtent.Expand(contribExtent)
		//PBRT_IRRADIANCE_CACHE_ADDED_NEW_SAMPLE(const_cast<Point *>(&p), const_cast<Normal *>(&ng), contribExtent, &E, &wAvg, pixelSpacing);

		// Allocate _IrradianceSample_, get write lock, add to octree
		sample := &irradianceSample{*E, *ng, *p, *wAvg, contribExtent}
		//RWMutexLock lock(*mutex, WRITE);
		integrator.octree.Add(sample, sampleExtent)
		wi = wAvg
	}

	// Compute reflected radiance due to irradiance and BSDF
	if wi.LengthSquared() == 0.0 {
		return NewSpectrum1(0.0)
	}
	return bsdf.f(wo, NormalizeVector(wi), flags).Mult(E)
}

func (integrator *IrradianceCacheIntegrator) interpolateE(scene *Scene,
	p *Point, n *Normal) (ok bool, E *Spectrum, wi *Vector) {
	if integrator.octree == nil {
		return false, nil, nil
	}
	//PBRT_IRRADIANCE_CACHE_STARTED_INTERPOLATION(const_cast<Point *>(&p), const_cast<Normal *>(&n));
	proc := &irradProcess{*p, *n, integrator.minWeight, integrator.cosMaxSampleAngleDifference, 0, 0, *NewSpectrum1(0.0), Vector{0, 0, 0}}
	//RWMutexLock lock(*mutex, READ);

	integrator.octree.Lookup(p, proc.procFunc)
	//PBRT_IRRADIANCE_CACHE_FINISHED_INTERPOLATION(const_cast<Point *>(&p), const_cast<Normal *>(&n), proc.Successful() ? 1 : 0, proc.nFound);
	if !proc.Successful() {
		return false, nil, nil
	}
	E = proc.GetIrradiance()
	wi = proc.GetAverageDirection()
	return true, E, wi
}

func (integrator *IrradianceCacheIntegrator) pathL(r *RayDifferential, scene *Scene,
	renderer Renderer, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	pathThroughput := NewSpectrum1(1.0)
	ray := r
	specularBounce := false
	for pathLength := 0; ; pathLength++ {
		// Find next vertex of path
		var isect *Intersection
		var hit bool
		if hit, isect = scene.Intersect(ray); !hit {
			break
		}
		if pathLength == 0 {
			r.Maxt = ray.Maxt
		}
		pathThroughput = pathThroughput.Mult(renderer.Transmittance(scene, ray, nil, rng, arena))
		// Possibly add emitted light at path vertex
		if specularBounce {
			L = L.Add(pathThroughput.Mult(isect.Le(ray.Dir.Negate())))
		}
		// Evaluate BSDF at hit point
		bsdf := isect.GetBSDF(ray, arena)
		// Sample illumination from lights to find path contribution
		p := bsdf.dgShading.p
		n := bsdf.dgShading.nn
		wo := ray.Dir.Negate()
		L = L.Add(pathThroughput.Mult(
			UniformSampleOneLight(scene, renderer, arena, p, n, wo, isect.rayEpsilon,
				ray.Time, bsdf, nil, rng, 0, nil, nil)))
		if pathLength+1 == integrator.maxIndirectDepth {
			break
		}
		// Sample BSDF to get new path direction
		// Get random numbers for sampling new direction, \mono{bs1}, \mono{bs2}, and \mono{bcs}
		f, wi, pdf, flags := bsdf.Sample_f(wo, CreateRandomBSDFSample(rng), BSDF_ALL)
		if f.IsBlack() || pdf == 0.0 {
			break
		}
		specularBounce = (flags & BSDF_SPECULAR) != 0
		pathThroughput = pathThroughput.Mult(f.Scale(AbsDotVectorNormal(wi, n) / pdf))
		ray = CreateChildRayDifferentialFromRayDifferential(p, wi, ray, isect.rayEpsilon, INFINITY)
		// Possibly terminate the path
		if pathLength > 2 {
			rrProb := math.Min(1.0, pathThroughput.Y())
			if rng.RandomFloat() > rrProb {
				break
			}
			pathThroughput = pathThroughput.InvScale(rrProb)
		}
	}
	return L
}

func newIrradiancePrimeTask(scene *Scene, renderer Renderer, camera Camera, sampler Sampler, sample *Sample,
	integrator *IrradianceCacheIntegrator, prog *ProgressReporter, taskNum, numTasks int) *irradiancePrimeTask {
	task := new(irradiancePrimeTask)
	task.scene = scene
	task.camera = camera
	task.renderer = renderer
	task.sampler = sampler.GetSubSampler(taskNum, numTasks)
	task.origSample = sample
	task.irradianceCache = integrator
	task.progress = prog
	task.taskNum = taskNum
	task.numTasks = numTasks
	return task
}

func (primeTask *irradiancePrimeTask) run() {
	if primeTask.sampler == nil {
		primeTask.progress.Update(1)
		return
	}
	rng := NewRNG(int64(29 * primeTask.taskNum))
	maxSamples := primeTask.sampler.MaximumSampleCount()
	samples := primeTask.origSample.Duplicate(maxSamples)
	for {
		if sampleCount := primeTask.sampler.GetMoreSamples(&samples, rng); sampleCount > 0 {
			for i := 0; i < sampleCount; i++ {
				ray, _ := GenerateRayDifferential(primeTask.camera, &samples[i])
				if ok, isect := primeTask.scene.Intersect(ray); ok {
					primeTask.irradianceCache.Li(primeTask.scene, primeTask.renderer, ray, isect, &samples[i], rng, nil)
				}
			}
		} else {
			break
		}
	}
	primeTask.progress.Update(1)
}

func (ip *irradProcess) Successful() bool {
	return ip.sumWt >= ip.minWeight
}
func (ip *irradProcess) GetIrradiance() *Spectrum {
	return ip.E.InvScale(ip.sumWt)
}
func (ip *irradProcess) GetAverageDirection() *Vector {
	return &ip.wAvg
}
func (ip *irradProcess) procFunc(data Object) bool {
	sample, ok := data.(*irradianceSample)
	if ok {
		// Compute estimate error term and possibly use sample
		perr := DistancePoint(&ip.p, &sample.p) / sample.maxDist
		nerr := math.Sqrt((1.0 - DotNormal(&ip.n, &sample.n)) /
			(1.0 - ip.cosMaxSampleAngleDifference))
		err := math.Max(perr, nerr)
		//PBRT_IRRADIANCE_CACHE_CHECKED_SAMPLE(const_cast<IrradianceSample *>(sample), perr, nerr);
		if err < 1.0 {
			ip.nFound++
			wt := 1.0 - err
			ip.E = *ip.E.Add(sample.E.Scale(wt))
			ip.wAvg = *ip.wAvg.Add(sample.wAvg.Scale(wt))
			ip.sumWt += wt
		}
	}
	return true
}
