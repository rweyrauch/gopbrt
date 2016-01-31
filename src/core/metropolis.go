package core

import (
	"math"
)

type (
	MetropolisRenderer struct {
		camera                                                Camera
		bidirectional                                         bool
		nDirectPixelSamples, nPixelSamples, maxDepth          int
		largeStepsPerPixel, nBootstrap, maxConsecutiveRejects int
		directLighting                                        *DirectLightingIntegrator
		nTasksFinished                                        int
	}

	// Metropolis Local Declarations
	PathSample struct {
		bsdfSample BSDFSample
		rrSample   float64
	}

	LightingSample struct {
		bsdfSample  BSDFSample
		lightNum    float64
		lightSample LightSample
	}

	MLTSample struct {
		cameraSample                        Sample
		lightNumSample                      float64
		lightRaySamples                     [5]float64
		cameraPathSamples, lightPathSamples []PathSample
		lightingSamples                     []LightingSample
	}

	PathVertex struct {
		isect               Intersection
		wPrev, wNext        Vector
		bsdf                *BSDF
		specularBounce      bool
		nSpecularComponents int
		alpha               Spectrum
	}

	MLTTask struct {
		progress                         *ProgressReporter
		progressUpdateFrequency, taskNum int
		dx, dy                           float64
		currentPixelSample               int
		x0, x1, y0, y1                   int
		t0, t1                           float64
		b                                float64
		initialSample                    *MLTSample
		scene                            *Scene
		camera                           Camera
		renderer                         *MetropolisRenderer
		//Mutex *filmMutex;
		lightDistribution *Distribution1D
	}
)

func NewMetropolisRenderer(perPixelSamples, nBootstrap, nDirectPixelSamples int, largeStepProbability float64, doDirectSeparately bool,
	maxrejects, maxdepth int, camera Camera, doBidirectional bool) *MetropolisRenderer {
	renderer := new(MetropolisRenderer)
	renderer.camera = camera
	renderer.bidirectional = doBidirectional
	renderer.nDirectPixelSamples = nDirectPixelSamples
	renderer.nPixelSamples = perPixelSamples
	renderer.maxDepth = maxdepth
	renderer.nBootstrap = nBootstrap
	renderer.maxConsecutiveRejects = maxrejects

	renderer.largeStepsPerPixel = Maxi(1, int(RoundUpPow2(uint32(largeStepProbability*float64(renderer.nPixelSamples)))))
	if renderer.largeStepsPerPixel >= renderer.nPixelSamples {
		renderer.largeStepsPerPixel /= 2
	}
	Assert(renderer.largeStepsPerPixel >= 1 && renderer.largeStepsPerPixel < renderer.nPixelSamples)
	if (renderer.nPixelSamples % renderer.largeStepsPerPixel) != 0 {
		origPixelSamples := renderer.nPixelSamples
		renderer.nPixelSamples += renderer.largeStepsPerPixel - (renderer.nPixelSamples % renderer.largeStepsPerPixel)
		Warning("Rounding up to %d Metropolis samples per pixel (from %d)",
			renderer.nPixelSamples, origPixelSamples)
	}

	if doDirectSeparately {
		renderer.directLighting = &DirectLightingIntegrator{SAMPLE_ALL_UNIFORM, renderer.maxDepth, nil, nil, 0}
	}

	return renderer
}

func (r *MetropolisRenderer) Render(scene *Scene) {
	//PBRT_MLT_STARTED_RENDERING();
	if len(scene.lights) > 0 {
		x0, x1, y0, y1 := r.camera.Film().GetPixelExtent()
		t0, t1 := r.camera.ShutterOpen(), r.camera.ShutterClose()
		lightDistribution := ComputeLightSamplingCDF(scene)

		if r.directLighting != nil {
			//PBRT_MLT_STARTED_DIRECTLIGHTING();
			// Compute direct lighting before Metropolis light transport
			if r.nDirectPixelSamples > 0 {
				sampler := NewLDSampler(x0, x1, y0, y1, r.nDirectPixelSamples, t0, t1)
				sample := NewSample(sampler, r.directLighting, nil, scene)

				nDirectTasks := Maxi(32*NumSystemCores(),
					(r.camera.Film().XResolution()*r.camera.Film().YResolution())/(16*16))
				nDirectTasks = int(RoundUpPow2(uint32(nDirectTasks)))
				directProgress := NewProgressReporter(nDirectTasks, "Direct Lighting", TerminalWidth())
				for i := 0; i < nDirectTasks; i++ {
					srt := newSamplerRendererTask(scene, r, r.camera, directProgress, sampler, sample, false, i, nDirectTasks)
					srt.run()
				}
				//std::reverse(directTasks.begin(), directTasks.end())
				//EnqueueTasks(directTasks)
				//WaitForAllTasks()
				directProgress.Done()

			}
			r.camera.Film().WriteImage(1.0)
			//PBRT_MLT_FINISHED_DIRECTLIGHTING();
		}
		// Take initial set of samples to compute $b$
		//PBRT_MLT_STARTED_BOOTSTRAPPING(nBootstrap);
		rng := NewRNG(0)

		var arena *MemoryArena
		cameraPath := make([]PathVertex, r.maxDepth, r.maxDepth)
		lightPath := make([]PathVertex, r.maxDepth, r.maxDepth)
		sumI := 0.0

		bootstrapI := make([]float64, 0, r.nBootstrap)
		sample := NewMLTSample(r.maxDepth)
		for i := 0; i < r.nBootstrap; i++ {
			// Generate random sample and path radiance for MLT bootstrapping
			x := Lerp(rng.RandomFloat(), float64(x0), float64(x1))
			y := Lerp(rng.RandomFloat(), float64(y0), float64(y1))
			LargeStep(rng, sample, r.maxDepth, x, y, t0, t1, r.bidirectional)
			L := r.PathL(sample, scene, arena, r.camera, lightDistribution,
				cameraPath, lightPath, rng)

			// Compute contribution for random sample for MLT bootstrapping
			I := L.Y()

			sumI += I
			bootstrapI = append(bootstrapI, I)
			//arena.FreeAll()
		}

		b := sumI / float64(r.nBootstrap)
		//PBRT_MLT_FINISHED_BOOTSTRAPPING(b);
		Info("MLT computed b = %f", b)

		// Select initial sample from bootstrap samples
		contribOffset := rng.RandomFloat() * sumI
		rng.Seed(0)
		sumI = 0.0

		initialSample := NewMLTSample(r.maxDepth)
		for i := 0; i < r.nBootstrap; i++ {
			x := Lerp(rng.RandomFloat(), float64(x0), float64(x1))
			y := Lerp(rng.RandomFloat(), float64(y0), float64(y1))
			LargeStep(rng, initialSample, r.maxDepth, x, y, t0, t1, r.bidirectional)
			sumI += bootstrapI[i]
			if sumI > contribOffset {
				break
			}
		}

		// Launch tasks to generate Metropolis samples
		nTasks := r.largeStepsPerPixel
		largeStepRate := r.nPixelSamples / r.largeStepsPerPixel
		Info("MLT running %d tasks, large step rate %d", nTasks, largeStepRate)
		/*
		   ProgressReporter progress(nTasks * largeStepRate, "Metropolis");
		   vector<Task *> tasks;
		   Mutex *filmMutex = Mutex::Create();
		   Assert(IsPowerOf2(nTasks));
		   uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
		   uint32_t pfreq = (x1-x0) * (y1-y0);
		   for (uint32_t i = 0; i < nTasks; ++i) {
		       float d[2];
		       Sample02(i, scramble, d);
		       tasks.push_back(new MLTTask(progress, pfreq, i,
		           d[0], d[1], x0, x1, y0, y1, t0, t1, b, initialSample,
		           scene, camera, this, filmMutex, lightDistribution));
		   }
		   EnqueueTasks(tasks);
		   WaitForAllTasks();
		   for (uint32_t i = 0; i < tasks.size(); ++i)
		       delete tasks[i];
		   progress.Done();
		   Mutex::Destroy(filmMutex);
		   delete lightDistribution;
		*/
	}
	r.camera.Film().WriteImage(1.0)
	//PBRT_MLT_FINISHED_RENDERING();
}

func (r *MetropolisRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) {
	Lo := NewSpectrum1(0.0)
	var hit bool
	if hit, isect = scene.Intersect(ray); hit {
		Lo = r.directLighting.Li(scene, r, ray, isect, sample, rng, arena)
	} else {
		// Handle ray that doesn't intersect any geometry
		for i := 0; i < len(scene.lights); i++ {
			Lo = Lo.Add(scene.lights[i].Le(ray))
		}
	}
	return Lo, isect, NewSpectrum1(1.0)
}

func (*MetropolisRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return NewSpectrum1(1.0)
}

func (r *MetropolisRenderer) PathL(sample *MLTSample, scene *Scene, arena *MemoryArena, camera Camera,
	lightDistribution *Distribution1D, cameraPath []PathVertex, lightPath []PathVertex, rng *RNG) *Spectrum {
	// Generate camera path from camera path samples
	//PBRT_STARTED_GENERATING_CAMERA_RAY((CameraSample *)(&sample.cameraSample));
	cameraRay, cameraWt := camera.GenerateRayDifferential(&sample.cameraSample)
	cameraRay.ScaleDifferentials(1.0 / math.Sqrt(float64(r.nPixelSamples)))
	//PBRT_FINISHED_GENERATING_CAMERA_RAY((CameraSample *)(&sample.cameraSample), &cameraRay, cameraWt);
	cameraLength, escapedRay, escapedAlpha := GeneratePath(cameraRay, NewSpectrum1(cameraWt), scene, arena,
		sample.cameraPathSamples, &cameraPath)
	if !r.bidirectional {
		// Compute radiance along path using path tracing
		return r.Lpath(scene, cameraPath, cameraLength, arena,
			sample.lightingSamples, rng, sample.cameraSample.time,
			lightDistribution, escapedRay, escapedAlpha)
	} else {
		// Sample light ray and apply bidirectional path tracing

		// Choose light and sample ray to start light path
		//PBRT_MLT_STARTED_SAMPLE_LIGHT_FOR_BIDIR();
		lightNum, lightPdf := lightDistribution.SampleDiscrete(sample.lightNumSample)
		light := scene.lights[lightNum]
		var lrs LightSample
		lrs.uPos[0] = sample.lightRaySamples[0]
		lrs.uPos[1] = sample.lightRaySamples[1]
		lrs.uComponent = sample.lightRaySamples[2]
		lightWt, lightRay, Nl, lightRayPdf := light.Sample_L2(scene, &lrs, sample.lightRaySamples[3],
			sample.lightRaySamples[4], sample.cameraSample.time)
		//PBRT_MLT_FINISHED_SAMPLE_LIGHT_FOR_BIDIR();
		if lightWt.IsBlack() || lightRayPdf == 0.0 {
			// Compute radiance along path using path tracing
			return r.Lpath(scene, cameraPath, cameraLength, arena,
				sample.lightingSamples, rng, sample.cameraSample.time,
				lightDistribution, escapedRay, escapedAlpha)
		} else {
			// Compute radiance along paths using bidirectional path tracing
			lightWt = lightWt.Scale(AbsDotNormalVector(NormalizeNormal(Nl), &lightRay.dir) / (lightPdf * lightRayPdf))
			lightLength, _, _ := GeneratePath(CreateRayDifferentialFromRay(lightRay), lightWt,
				scene, arena, sample.lightPathSamples, &lightPath)

			return r.Lbidir(scene, cameraPath, cameraLength, lightPath, lightLength,
				arena, sample.lightingSamples, rng, sample.cameraSample.time,
				lightDistribution, escapedRay, escapedAlpha)
		}
	}
}

func (r *MetropolisRenderer) Lpath(scene *Scene, cameraPath []PathVertex, cameraPathLength int,
	arena *MemoryArena, samples []LightingSample, rng *RNG, time float64, lightDistribution *Distribution1D,
	eRay *RayDifferential, eAlpha *Spectrum) *Spectrum {
	//PBRT_MLT_STARTED_LPATH();
	L := NewSpectrum1(0.0)
	previousSpecular, allSpecular := true, true
	for i := 0; i < cameraPathLength; i++ {
		// Initialize basic variables for camera path vertex
		vc := &cameraPath[i]
		pc := vc.bsdf.dgShading.p
		nc := vc.bsdf.dgShading.nn

		// Add emitted light from vertex if appropriate
		if previousSpecular && (r.directLighting == nil || !allSpecular) {
			L = L.Add(vc.alpha.Mult(vc.isect.Le(&vc.wPrev)))
		}
		// Compute direct illumination for Metropolis path vertex
		Ld := NewSpectrum1(0.0)
		if r.directLighting == nil || !allSpecular {
			// Choose light and call _EstimateDirect()_ for Metropolis vertex
			ls := &samples[i]
			lightNum, lightPdf := lightDistribution.SampleDiscrete(ls.lightNum)
			light := scene.lights[lightNum]
			//PBRT_MLT_STARTED_ESTIMATE_DIRECT();

			Ld = vc.alpha.Mult(EstimateDirect(scene, r, arena, light, pc, nc, &vc.wPrev,
				vc.isect.rayEpsilon, time, vc.bsdf, rng,
				&ls.lightSample, &ls.bsdfSample,
				BxDFType(BSDF_ALL^BSDF_SPECULAR)).InvScale(lightPdf))
			//PBRT_MLT_FINISHED_ESTIMATE_DIRECT();
		}
		previousSpecular = vc.specularBounce
		allSpecular = allSpecular && previousSpecular
		L = L.Add(Ld)
	}
	// Add contribution of escaped ray, if any
	if !eAlpha.IsBlack() && previousSpecular && (r.directLighting == nil || !allSpecular) {
		for i := 0; i < len(scene.lights); i++ {
			L = L.Add(eAlpha.Mult(scene.lights[i].Le(eRay)))
		}
	}
	//PBRT_MLT_FINISHED_LPATH();
	return L
}

func (r *MetropolisRenderer) Lbidir(scene *Scene, cameraPath []PathVertex, cameraPathLength int,
	lightPath []PathVertex, lightPathLength int,
	arena *MemoryArena, samples []LightingSample,
	rng *RNG, time float64, lightDistribution *Distribution1D,
	eRay *RayDifferential, eAlpha *Spectrum) *Spectrum {
	//PBRT_MLT_STARTED_LBIDIR();
	L := NewSpectrum1(0.0)
	previousSpecular, allSpecular := true, true
	// Compute number of specular vertices for each path length
	nVerts := cameraPathLength + lightPathLength + 2
	nSpecularVertices := make([]int, nVerts, nVerts)
	for i := 0; i < cameraPathLength; i++ {
		for j := 0; j < lightPathLength; j++ {
			if cameraPath[i].specularBounce || lightPath[j].specularBounce {
				nSpecularVertices[i+j+2]++
			}
		}
	}
	for i := 0; i < cameraPathLength; i++ {
		// Initialize basic variables for camera path vertex
		vc := &cameraPath[i]
		pc := vc.bsdf.dgShading.p
		nc := vc.bsdf.dgShading.nn

		// Compute reflected light at camera path vertex

		// Add emitted light from vertex if appropriate
		if previousSpecular && (r.directLighting == nil || !allSpecular) {
			L = L.Add(vc.alpha.Mult(vc.isect.Le(&vc.wPrev)))
		}
		// Compute direct illumination for Metropolis path vertex
		Ld := NewSpectrum1(0.0)
		if r.directLighting == nil || !allSpecular {
			// Choose light and call _EstimateDirect()_ for Metropolis vertex
			ls := &samples[i]
			lightNum, lightPdf := lightDistribution.SampleDiscrete(ls.lightNum)
			light := scene.lights[lightNum]
			//PBRT_MLT_STARTED_ESTIMATE_DIRECT();

			Ld = vc.alpha.Mult(EstimateDirect(scene, r, arena, light, pc, nc, &vc.wPrev,
				vc.isect.rayEpsilon, time, vc.bsdf, rng,
				&ls.lightSample, &ls.bsdfSample,
				BxDFType(BSDF_ALL^BSDF_SPECULAR)).InvScale(lightPdf))
			//PBRT_MLT_FINISHED_ESTIMATE_DIRECT();
		}
		previousSpecular = vc.specularBounce
		allSpecular = allSpecular && previousSpecular
		L = L.Add(Ld.InvScale(float64(i + 1 - nSpecularVertices[i+1])))
		if !vc.specularBounce {
			// Loop over light path vertices and connect to camera vertex
			for j := 0; j < lightPathLength; j++ {
				vl := &lightPath[j]
				pl := vl.bsdf.dgShading.p
				nl := vl.bsdf.dgShading.nn
				if !vl.specularBounce {
					// Compute contribution between camera and light vertices
					w := NormalizeVector(pl.Sub(pc))
					fc := vc.bsdf.f(&vc.wPrev, w, BSDF_ALL).Scale(float64(1 + vc.nSpecularComponents))
					fl := vl.bsdf.f(w.Negate(), &vl.wPrev, BSDF_ALL).Scale(float64(1 + vl.nSpecularComponents))
					if fc.IsBlack() || fl.IsBlack() {
						continue
					}
					ray := CreateRay(pc, pl.Sub(pc), 1.0e-3, 0.999, time, 0)
					if !scene.IntersectP(ray) {
						// Compute weight for bidirectional path, _pathWt_
						pathWt := 1.0 / float64(i+j+2-nSpecularVertices[i+j+2])
						G := AbsDotNormalVector(nc, w) * AbsDotNormalVector(nl, w) / DistanceSquaredPoint(pl, pc)
						L = L.Add((vc.alpha.Mult(fc).Mult(fl.Scale(G)).Mult(&vl.alpha)).Scale(pathWt))
					}
				}
			}
		}
	}
	// Add contribution of escaped ray, if any
	if !eAlpha.IsBlack() && previousSpecular && (r.directLighting == nil || !allSpecular) {
		for i := 0; i < len(scene.lights); i++ {
			L = L.Add(eAlpha.Mult(scene.lights[i].Le(eRay)))
		}
	}
	//PBRT_MLT_FINISHED_LBIDIR();
	return L
}

func CreateMetropolisRenderer(params *ParamSet, camera Camera) *MetropolisRenderer {
	largeStepProbability := params.FindFloatParam("largestepprobability", 0.25)
	perPixelSamples := params.FindIntParam("samplesperpixel", 100)
	nBootstrap := params.FindIntParam("bootstrapsamples", 100000)
	nDirectPixelSamples := params.FindIntParam("directsamples", 4)
	doDirectSeparately := params.FindBoolParam("dodirectseparately", true)
	mr := params.FindIntParam("maxconsecutiverejects", 512)
	md := params.FindIntParam("maxdepth", 7)
	doBidirectional := params.FindBoolParam("bidirectional", true)

	if options.QuickRender {
		perPixelSamples = Maxi(1, perPixelSamples/4)
		nBootstrap = Maxi(1, nBootstrap/4)
		nDirectPixelSamples = Maxi(1, nDirectPixelSamples/4)
	}

	return NewMetropolisRenderer(perPixelSamples, nBootstrap,
		nDirectPixelSamples, largeStepProbability, doDirectSeparately,
		mr, md, camera, doBidirectional)
}

func NewMLTSample(maxLength int) *MLTSample {
	sample := new(MLTSample)
	sample.cameraPathSamples = make([]PathSample, 0, maxLength)
	sample.lightPathSamples = make([]PathSample, 0, maxLength)
	sample.lightingSamples = make([]LightingSample, 0, maxLength)
	return sample
}

func LargeStep(rng *RNG, sample *MLTSample, maxDepth int, x, y, t0, t1 float64, bidirectional bool) {
	// Do large step mutation of _cameraSample_
	sample.cameraSample.imageX = x
	sample.cameraSample.imageY = y
	sample.cameraSample.time = Lerp(rng.RandomFloat(), t0, t1)
	sample.cameraSample.lensU = rng.RandomFloat()
	sample.cameraSample.lensV = rng.RandomFloat()
	for i := 0; i < maxDepth; i++ {
		// Apply large step to $i$th camera _PathSample_
		cps := &sample.cameraPathSamples[i]
		cps.bsdfSample.uComponent = rng.RandomFloat()
		cps.bsdfSample.uDir[0] = rng.RandomFloat()
		cps.bsdfSample.uDir[1] = rng.RandomFloat()
		cps.rrSample = rng.RandomFloat()

		// Apply large step to $i$th _LightingSample_
		ls := &sample.lightingSamples[i]
		ls.bsdfSample.uComponent = rng.RandomFloat()
		ls.bsdfSample.uDir[0] = rng.RandomFloat()
		ls.bsdfSample.uDir[1] = rng.RandomFloat()
		ls.lightNum = rng.RandomFloat()
		ls.lightSample.uComponent = rng.RandomFloat()
		ls.lightSample.uPos[0] = rng.RandomFloat()
		ls.lightSample.uPos[1] = rng.RandomFloat()
	}
	if bidirectional {
		// Apply large step to bidirectional light samples
		sample.lightNumSample = rng.RandomFloat()
		for i := 0; i < 5; i++ {
			sample.lightRaySamples[i] = rng.RandomFloat()
		}
		for i := 0; i < maxDepth; i++ {
			// Apply large step to $i$th light _PathSample_
			lps := &sample.lightPathSamples[i]
			lps.bsdfSample.uComponent = rng.RandomFloat()
			lps.bsdfSample.uDir[0] = rng.RandomFloat()
			lps.bsdfSample.uDir[1] = rng.RandomFloat()
			lps.rrSample = rng.RandomFloat()
		}
	}
}

const (
	a = 1.0 / 1024.0
	b = 1.0 / 64.0
)
	
var	logRatio float64 = -math.Log(b/a)


func mutate(rng *RNG, v *float64, min, max float64) {
    if min == max { *v = min; return }
    Assert(min < max)
    delta := (max - min) * b * math.Exp(logRatio * rng.RandomFloat())
    if rng.RandomFloat() < 0.5 {
        *v += delta
        if *v >= max { *v = min + (*v - max) }
    } else {
        *v -= delta
        if *v < min { *v = max - (min - *v) }
    }
    if *v < min || *v >= max { *v = min }
}


func SmallStep(rng *RNG, sample *MLTSample, maxDepth, x0, x1, y0, y1 int, t0, t1 float64,
        bidirectional bool) {
    mutate(rng, &sample.cameraSample.imageX, float64(x0), float64(x1))
    mutate(rng, &sample.cameraSample.imageY, float64(y0), float64(y1))
    mutate(rng, &sample.cameraSample.time, t0, t1)
    mutate(rng, &sample.cameraSample.lensU, 0.0, 1.0)
    mutate(rng, &sample.cameraSample.lensV, 0.0, 1.0)
    // Apply small step mutation to camera, lighting, and light samples
    for i := 0; i < maxDepth; i++ {
        // Apply small step to $i$th camera _PathSample_
        eps := &sample.cameraPathSamples[i]
        mutate(rng, &eps.bsdfSample.uComponent, 0.0, 1.0)
        mutate(rng, &eps.bsdfSample.uDir[0], 0.0, 1.0)
        mutate(rng, &eps.bsdfSample.uDir[1], 0.0, 1.0)
        mutate(rng, &eps.rrSample, 0.0, 1.0)

        // Apply small step to $i$th _LightingSample_
        ls := &sample.lightingSamples[i]
        mutate(rng, &ls.bsdfSample.uComponent, 0.0, 1.0)
        mutate(rng, &ls.bsdfSample.uDir[0], 0.0, 1.0)
        mutate(rng, &ls.bsdfSample.uDir[1], 0.0, 1.0)
        mutate(rng, &ls.lightNum, 0.0, 1.0)
        mutate(rng, &ls.lightSample.uComponent, 0.0, 1.0)
        mutate(rng, &ls.lightSample.uPos[0], 0.0, 1.0)
        mutate(rng, &ls.lightSample.uPos[1], 0.0, 1.0)
    }
    
    if bidirectional {
        mutate(rng, &sample.lightNumSample, 0.0, 1.0)
        for i := 0; i < 5; i++ {
            mutate(rng, &sample.lightRaySamples[i], 0.0, 1.0)
		}            
        for i := 0; i < maxDepth; i++ {
            // Apply small step to $i$th light _PathSample_
            lps := &sample.lightPathSamples[i]
            mutate(rng, &lps.bsdfSample.uComponent, 0.0, 1.0)
            mutate(rng, &lps.bsdfSample.uDir[0], 0.0, 1.0)
            mutate(rng, &lps.bsdfSample.uDir[1], 0.0, 1.0)
            mutate(rng, &lps.rrSample, 0.0, 1.0)
        }
    }
}


func GeneratePath(r *RayDifferential, a *Spectrum, scene *Scene, arena *MemoryArena, samples []PathSample,
	path *[]PathVertex) (pathLen int, escapedRay *RayDifferential, escapedAlpha *Spectrum) {
	//PBRT_MLT_STARTED_GENERATE_PATH();
	ray := r
	alpha := a
	escapedAlpha = NewSpectrum1(0.0)
	length := 0
	var hit bool
	var isect *Intersection
	for ; length < len(samples); length++ {
		// Try to generate next vertex of ray path
		v := (*path)[length]
		if hit, isect = scene.Intersect(ray); !hit {
			// Handle ray that leaves the scene during path generation
			escapedAlpha = alpha
			escapedRay = ray
			break
		}

		// Record information for current path vertex
		v.isect = *isect
		v.alpha = *alpha
		bsdf := v.isect.GetBSDF(ray, arena)
		v.bsdf = bsdf
		v.wPrev = *ray.dir.Negate()

		// Sample direction for outgoing Metropolis path direction
		f, wNext, pdf, flags := bsdf.Sample_f(ray.dir.Negate(), &samples[length].bsdfSample, BSDF_ALL)
		v.wNext = *wNext
		v.specularBounce = (flags & BSDF_SPECULAR) != 0
		v.nSpecularComponents = bsdf.NumComponentsMatching(BxDFType(BSDF_SPECULAR | BSDF_REFLECTION | BSDF_TRANSMISSION))
		if f.IsBlack() || pdf == 0.0 {
			//PBRT_MLT_FINISHED_GENERATE_PATH();
			return length + 1, escapedRay, escapedAlpha
		}

		// Terminate path with RR or prepare for finding next vertex
		p := bsdf.dgShading.p
		n := bsdf.dgShading.nn
		pathScale := f.Scale(AbsDotVectorNormal(&v.wNext, n) / pdf)
		rrSurviveProb := math.Min(1.0, pathScale.Y())
		if samples[length].rrSample > rrSurviveProb {
			//PBRT_MLT_FINISHED_GENERATE_PATH();
			return length + 1, escapedRay, escapedAlpha
		}
		alpha = alpha.Mult(pathScale.InvScale(rrSurviveProb))
		//alpha *= renderer->Transmittance(scene, ray, NULL, rng, arena);
		ray = CreateChildRayDifferential(p, &v.wNext, CreateRayFromRayDifferential(ray), v.isect.rayEpsilon, INFINITY)
	}
	//PBRT_MLT_FINISHED_GENERATE_PATH()
	return length, escapedRay, escapedAlpha
}

func NewMLTTask(prog *ProgressReporter, pfreq, taskNum int,
	ddx, ddy float64, xx0, xx1, yy0, yy1 int, tt0, tt1, bb float64, is *MLTSample, sc *Scene, camera Camera,
	renderer *MetropolisRenderer, lightDistribution *Distribution1D) *MLTTask {
		
	mlt := new(MLTTask)
	
    mlt.progress = prog
    mlt.initialSample = is
    	
    mlt.progressUpdateFrequency = pfreq
    mlt.taskNum = taskNum
    mlt.dx = ddx
    mlt.dy = ddy
    mlt.x0 = xx0
    mlt.x1 = xx1
    mlt.y0 = yy0
    mlt.y1 = yy1
    mlt.t0 = tt0
    mlt.t1 = tt1
    mlt.currentPixelSample = 0
    mlt.b = bb
    mlt.scene = sc
    mlt.camera = camera
    mlt.renderer = renderer
    //mlt.filmMutex = fm;
    mlt.lightDistribution = lightDistribution
    
    return mlt
}

func (mlt *MLTTask) run() {
    //PBRT_MLT_STARTED_MLT_TASK(this);
    // Declare basic _MLTTask_ variables and prepare for sampling
    //PBRT_MLT_STARTED_TASK_INIT();
    nPixels := (mlt.x1-mlt.x0) * (mlt.y1-mlt.y0)
    nPixelSamples := mlt.renderer.nPixelSamples
    largeStepRate := nPixelSamples / mlt.renderer.largeStepsPerPixel
    Assert(largeStepRate > 1)
    nTaskSamples := nPixels * largeStepRate
    consecutiveRejects := 0
    progressCounter := mlt.progressUpdateFrequency

    // Declare variables for storing and computing MLT samples
    var arena *MemoryArena
    rng := NewRNG(int64(mlt.taskNum))
    cameraPath := make([]PathVertex, mlt.renderer.maxDepth, mlt.renderer.maxDepth) 
    lightPath := make([]PathVertex, mlt.renderer.maxDepth, mlt.renderer.maxDepth)
    var samples [2]MLTSample
    samples[0] = *NewMLTSample(mlt.renderer.maxDepth)
    samples[1] = *NewMLTSample(mlt.renderer.maxDepth)
 
    var L [2]Spectrum
    var I [2]float64
    var current uint = 0
    var proposed uint = 1

    // Compute _L[current]_ for initial sample
    samples[current] = *mlt.initialSample
    L[current] = *mlt.renderer.PathL(mlt.initialSample, mlt.scene, arena, mlt.camera,
                     mlt.lightDistribution, cameraPath, lightPath, rng)
    I[current] = L[current].Y()
    //arena.FreeAll();

    // Compute randomly permuted table of pixel indices for large steps
    pixelNumOffset := 0
    largeStepPixelNum := make([]int, 0, nPixels)
    for i := 0; i < nPixels; i++ { largeStepPixelNum = append(largeStepPixelNum, i) }
    //Shuffle(&largeStepPixelNum[0], nPixels, 1, rng)
    //PBRT_MLT_FINISHED_TASK_INIT();
    for s := 0; s < nTaskSamples; s++ {
        // Compute proposed mutation to current sample
        //PBRT_MLT_STARTED_MUTATION();
        samples[proposed] = samples[current]
        largeStep := ((s % largeStepRate) == 0)
        if largeStep {
            x := float64(mlt.x0 + largeStepPixelNum[pixelNumOffset] % (mlt.x1 - mlt.x0))
            y := float64(mlt.y0 + largeStepPixelNum[pixelNumOffset] / (mlt.x1 - mlt.x0))
            LargeStep(rng, &samples[proposed], mlt.renderer.maxDepth,
                      x + mlt.dx, y + mlt.dy, mlt.t0, mlt.t1, mlt.renderer.bidirectional)
            pixelNumOffset++
        } else {
            SmallStep(rng, &samples[proposed], mlt.renderer.maxDepth,
                      mlt.x0, mlt.x1, mlt.y0, mlt.y1, mlt.t0, mlt.t1, mlt.renderer.bidirectional)
		}            
        //PBRT_MLT_FINISHED_MUTATION();

        // Compute contribution of proposed sample
        L[proposed] = *mlt.renderer.PathL(&samples[proposed], mlt.scene, arena, mlt.camera,
                         mlt.lightDistribution, cameraPath, lightPath, rng)
        I[proposed] = L[proposed].Y()
        //arena.FreeAll();

        // Compute acceptance probability for proposed sample
        a := math.Min(1.0, I[proposed] / I[current])

        // Splat current and proposed samples to _Film_
        //PBRT_MLT_STARTED_SAMPLE_SPLAT();
        if I[current] > 0.0 {
            if !math.IsInf(1.0 / I[current], 0) {
            	contrib := L[current].InvScale(I[current]).Scale(b / float64(nPixelSamples))
            	mlt.camera.Film().Splat(&samples[current].cameraSample, contrib.Scale(1.0 - a))
        	}
        }
        if I[proposed] > 0.0 {
            if !math.IsInf(1.0 / I[proposed], 0) {
            	contrib :=  L[proposed].InvScale(I[proposed]).Scale(b / float64(nPixelSamples))
            	mlt.camera.Film().Splat(&samples[proposed].cameraSample, contrib.Scale(a))
        	}
        }
        //PBRT_MLT_FINISHED_SAMPLE_SPLAT();

        // Randomly accept proposed path mutation (or not)
        if consecutiveRejects >= mlt.renderer.maxConsecutiveRejects ||
            rng.RandomFloat() < a {
            //PBRT_MLT_ACCEPTED_MUTATION(a, &samples[current], &samples[proposed]);
            current ^= 1
            proposed ^= 1
            consecutiveRejects = 0
        } else {
            //PBRT_MLT_REJECTED_MUTATION(a, &samples[current], &samples[proposed]);
            consecutiveRejects++
        }
        progressCounter--
        if progressCounter == 0 {
            mlt.progress.Update(1)
            progressCounter = mlt.progressUpdateFrequency
        }
    }
    Assert(pixelNumOffset == nPixels)
    // Update display for recently computed Metropolis samples
    //PBRT_MLT_STARTED_DISPLAY_UPDATE();
    //int ntf = AtomicAdd(&renderer->nTasksFinished, 1);
    mlt.renderer.nTasksFinished++
    ntf := mlt.renderer.nTasksFinished
    totalSamples := nPixels * nPixelSamples
    splatScale := float64(totalSamples) / float64(ntf * nTaskSamples)
    mlt.camera.Film().UpdateDisplay(mlt.x0, mlt.y0, mlt.x1, mlt.y1, splatScale)
    if (mlt.taskNum % 8) == 0 {
        //MutexLock lock(*filmMutex);
        mlt.camera.Film().WriteImage(splatScale)
    }
    //PBRT_MLT_FINISHED_DISPLAY_UPDATE();
    //PBRT_MLT_FINISHED_MLT_TASK(this);
}
