package core

import (
	"math"
)

type (
	SamplerRenderer struct {
		visualizeObjectIds bool
		sampler            Sampler
		surfaceIntegrator  SurfaceIntegrator
		volumeIntegrator   VolumeIntegrator
		camera             Camera
	}
)

func NewSamplerRenderer(sampler Sampler, camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, visIds bool) *SamplerRenderer {
	renderer := new(SamplerRenderer)
	renderer.visualizeObjectIds = visIds
	renderer.sampler = sampler
	renderer.surfaceIntegrator = surf
	renderer.volumeIntegrator = vol
	renderer.camera = camera

	return renderer
}

func (r *SamplerRenderer) Render(scene *Scene) {
	//PBRT_FINISHED_PARSING();
	// Allow integrators to do preprocessing for the scene
	//PBRT_STARTED_PREPROCESSING();
	r.surfaceIntegrator.Preprocess(scene, r.camera, r)
	r.volumeIntegrator.Preprocess(scene, r.camera, r)
	//PBRT_FINISHED_PREPROCESSING();
	//PBRT_STARTED_RENDERING();
	// Allocate and initialize _sample_
	sample := NewSample(r.sampler, r.surfaceIntegrator, r.volumeIntegrator, scene)

	// Create and launch _SamplerRendererTask_s for rendering image

	// Compute number of _SamplerRendererTask_s to create for rendering
	nPixels := r.camera.Film().XResolution() * r.camera.Film().YResolution()
	nTasks := Maxi(32*NumSystemCores(), nPixels/(16*16))
	nTasks = int(RoundUpPow2(uint32(nTasks)))

	reporter := &ProgressReporter{}
	for i := 0; i < nTasks; i++ {
		// Create and 'run' the task synchronously for now.  Port to 'go' routine.
		task := newSamplerRendererTask(scene, r, r.camera, reporter, r.sampler, sample, r.visualizeObjectIds, nTasks-1-i, nTasks)
		task.run()
	}
	//EnqueueTasks(renderTasks);
	//WaitForAllTasks();
	//for (uint32_t i = 0; i < renderTasks.size(); ++i)
	//    delete renderTasks[i];
	reporter.Done()
	//PBRT_FINISHED_RENDERING();
	// Clean up after rendering and store final image
	//delete sample;
	r.camera.Film().WriteImage(1.0)
}

func (r *SamplerRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) {
	//Assert(ray.time == sample->time);
	//Assert(!ray.HasNaNs());
	// Allocate local variables for _isect_ and _T_ if needed
	li = NewSpectrum1(0.0)
	var hit bool
	// TODO: Must make RayBase an interface and define two structs (Ray, RayDifferential)
	if hit, isect = scene.Intersect(CreateRayFromRayDifferential(ray)); hit {
		li = r.surfaceIntegrator.Li(scene, r, ray, isect, sample, rng, arena)
	} else {
		// Handle ray that doesn't intersect any geometry
		for i := 0; i < len(scene.lights); i++ {
			li = li.Add(scene.lights[i].Le(ray))
		}
	}
	var Lvi *Spectrum
	Lvi, T = r.volumeIntegrator.Li(scene, r, ray, sample, rng, arena)
	li = li.Mult(T).Add(Lvi)
	return li, isect, T
}

func (r *SamplerRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return r.volumeIntegrator.Transmittance(scene, r, ray, sample, rng, arena)
}

func CreateSamplerRenderer(sampler Sampler, camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, visIds bool) *SamplerRenderer {
	return nil
}

type ProgressReporter struct{}

func (p *ProgressReporter) Update() {}
func (p *ProgressReporter) Done()   {}

type samplerRendererTask struct {
	scene              *Scene
	renderer           Renderer
	camera             Camera
	mainSampler        Sampler
	reporter           *ProgressReporter
	origSample         *Sample
	visualizeObjectIds bool
	taskNum, taskCount int
}

func newSamplerRendererTask(scene *Scene, renderer Renderer, camera Camera,
	progress *ProgressReporter, mainSampler Sampler, sample *Sample, visIds bool, taskNum, taskCount int) *samplerRendererTask {
	task := &samplerRendererTask{scene, renderer, camera, mainSampler, progress, sample, visIds, taskNum, taskCount}
	return task
}

func (t *samplerRendererTask) run() {
	//PBRT_STARTED_RENDERTASK(taskNum);
	// Get sub-_Sampler_ for _SamplerRendererTask_
	sampler := t.mainSampler.GetSubSampler(t.taskNum, t.taskCount)
	if sampler == nil {
		t.reporter.Update()
		//PBRT_FINISHED_RENDERTASK(taskNum)
		return
	}

	// Declare local variables used for rendering loop
	var arena *MemoryArena
	rng := NewRNG(int64(t.taskNum))

	// Allocate space for samples and intersections
	maxSamples := sampler.MaximumSampleCount()
	samples := t.origSample.Duplicate(maxSamples)
	rays := make([]RayDifferential, maxSamples, maxSamples)
	Ls := make([]Spectrum, maxSamples, maxSamples)
	Ts := make([]Spectrum, maxSamples, maxSamples)
	isects := make([]Intersection, maxSamples, maxSamples)

	// Get samples from _Sampler_ and update image
	sampleCount := sampler.GetMoreSamples(samples, rng)
	for sampleCount > 0 {
		Info("SamplerRender Task: %d Remaining Samples: %d", t.taskNum, sampleCount)
		// Generate camera rays and compute radiance along rays
		for i := 0; i < sampleCount; i++ {
			// Find camera ray for _sample[i]_
			//PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
			var rayWeight float64
			var ray *RayDifferential
			ray, rayWeight = t.camera.GenerateRayDifferential(&samples[i])
			rays[i] = *ray
			rays[i].ScaleDifferentials(1.0 / math.Sqrt(float64(sampler.SamplesPerPixel())))
			//PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight)

			// Evaluate radiance along camera ray
			//PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
			if t.visualizeObjectIds {
				if rayWeight > 0.0 {
					/*
										var hit bool
										hit, isects[i] = tscene.Intersect(rays[i])
										if hit {
											// random shading based on shape id...
											ids := [2]int{ isects[i].shapeId, isects[i].primitiveId };
											h := hash((char *)ids, sizeof(ids));
											float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff),
					                                     float((h >> 16) & 0xff) };
											Ls[i] = Spectrum::FromRGB(rgb);
											Ls[i] /= 255.0
										} else {
											Ls[i] = *NewSpectrum1(0.0)
										}
					*/
				} else {
					Ls[i] = *NewSpectrum1(0.0)
				}
			} else {
				if rayWeight > 0.0 {
					ls, isect, ts := t.renderer.Li(t.scene, &rays[i], &samples[i], rng, arena)
					Ls[i] = *ls
					isects[i] = *isect
					Ts[i] = *ts
					Ls[i] = *Ls[i].Scale(float32(rayWeight))
				} else {
					Ls[i] = *NewSpectrum1(0.0)
					Ts[i] = *NewSpectrum1(1.0)
				}

				// Issue warning if unexpected radiance value returned
				if Ls[i].HasNaNs() {
					Error("Not-a-number radiance value returned for image sample.  Setting to black.")
					Ls[i] = *NewSpectrum1(0.0)
				} else if Ls[i].Y() < -1.0e-5 {
					Error("Negative luminance value, %f, returned for image sample.  Setting to black.", Ls[i].Y())
					Ls[i] = *NewSpectrum1(0.0)
					//} else if math.IsInf(Ls[i].Y(),0) {
					//	Error("Infinite luminance value returned for image sample.  Setting to black.")
					//	Ls[i] = *NewSpectrum1(0.0)
				}
			}
			//PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
		}

		// Report sample results to _Sampler_, add contributions to image
		if sampler.ReportResults(samples, rays, Ls, isects, sampleCount) {
			for i := 0; i < sampleCount; i++ {
				//PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
				t.camera.Film().AddSample(&samples[i], &Ls[i])
				//PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
			}
		}

		// Free _MemoryArena_ memory from computing image sample values
		//arena.FreeAll();

		sampleCount = sampler.GetMoreSamples(samples, rng)
	}

    Info("SamplerRender Task: %d completed.", t.taskNum)
    
	// Clean up after _SamplerRendererTask_ is done with its image region
	xstart, xend, ystart, yend := sampler.PixelRegion()
	t.camera.Film().UpdateDisplay(xstart, ystart, xend+1, yend+1, 1.0)
	//delete sampler;
	//delete[] samples;
	//delete[] rays;
	//delete[] Ls;
	//delete[] Ts;
	//delete[] isects;
	t.reporter.Update()
	//PBRT_FINISHED_RENDERTASK(taskNum);
}
