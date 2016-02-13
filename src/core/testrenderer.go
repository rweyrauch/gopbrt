package core

import (
	"math"
)

type TestRenderer struct {
	sampler            Sampler
	surfaceIntegrator  SurfaceIntegrator
	volumeIntegrator   VolumeIntegrator
	camera             Camera
}

type taskInput struct {
	scene			  *Scene
	renderer		  Renderer
	sampler           Sampler
	origSample        *Sample
	camera            Camera
	taskNum, numTasks int
	reporter		  *ProgressReporter
}
//type taskOutput struct {
//	sample Sample
//	Li Spectrum
//}

func (renderer *TestRenderer) Render(scene *Scene) {
	// Allow integrators to do preprocessing for the scene
	renderer.surfaceIntegrator.Preprocess(scene, renderer.camera, renderer)
	renderer.volumeIntegrator.Preprocess(scene, renderer.camera, renderer)
	
	// Allocate and initialize _sample_
	sample := NewSample(renderer.sampler, renderer.surfaceIntegrator, renderer.volumeIntegrator, scene)

	// Create and launch _SamplerRendererTask_s for rendering image

	// Compute number of _SamplerRendererTask_s to create for rendering
	nPixels := renderer.camera.Film().XResolution() * renderer.camera.Film().YResolution()
	nWorkers := Maxi(NumSystemCores(), 1)
	nTasks := Maxi(32*NumSystemCores(), nPixels/(16*16))
	nTasks = int(RoundUpPow2(uint32(nTasks)))
	Info("Num tasks: %d Workers: %d Pixels: %d", nTasks, nWorkers, nPixels)
	
	jobs := make(chan taskInput, nTasks)
	producedSamples := make(chan taskOutput, 128)
	completed := make(chan bool, nTasks)

	for w := 0; w < nWorkers; w++ {
		go rendererWorker(jobs, producedSamples, completed)
	}
	reporter := NewProgressReporter(nTasks, "Rendering", -1)
	for i := 0; i < nTasks; i++ {
		sampler := renderer.sampler.GetSubSampler(i, nTasks)
		task := taskInput{scene, renderer, sampler, sample, renderer.camera, i, nTasks, reporter}
		jobs <- task
	}
	close(jobs)
	numCompleted := 0
	for numCompleted < nWorkers {
		select {
		case output := <-producedSamples:
			renderer.camera.Film().AddSample(&output.sample, &output.Li)
		case done := <-completed:
			if done {
				numCompleted++
			}
		default:
			break
		}
	}
	// close and drain the queue of any unprocessed samples
	close(producedSamples)
	for output := range producedSamples {
		renderer.camera.Film().AddSample(&output.sample, &output.Li)		
	}
	
	reporter.Done()

	// Clean up after rendering and store final image
	renderer.camera.Film().WriteImage(1.0)
}

func (renderer *TestRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) {
	// Allocate local variables for _isect_ and _T_ if needed
	li = NewSpectrum1(0.0)

	var hit bool
	if hit, isect = scene.Intersect(ray); hit {
		li = renderer.surfaceIntegrator.Li(scene, renderer, ray, isect, sample, rng, arena)
	} else {
		// Handle ray that doesn't intersect any geometry
		for i := 0; i < len(scene.lights); i++ {
			li = li.Add(scene.lights[i].Le(ray))
		}
	}

	lvi, ti := renderer.volumeIntegrator.Li(scene, renderer, ray, sample, rng, arena)
	li = li.Mult(ti).Add(lvi)
	return li, isect, ti
}

func (renderer *TestRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return renderer.volumeIntegrator.Transmittance(scene, renderer, ray, sample, rng, arena)
}

func rendererWorker(workQueue <-chan taskInput, results chan<- taskOutput, completed chan<- bool) {

	for work := range workQueue {
		if work.sampler == nil {
			work.reporter.Update(1)
			continue
		}

		rng := NewRNG(int64(work.taskNum))
		
		maxSamples := work.sampler.MaximumSampleCount()
		samples := work.origSample.Duplicate(maxSamples)
		rays := make([]*RayDifferential, maxSamples, maxSamples)
		Ls := make([]*Spectrum, maxSamples, maxSamples)
		Ts := make([]*Spectrum, maxSamples, maxSamples)
		isects := make([]*Intersection, maxSamples, maxSamples)

		// Get samples from _Sampler_ and update image
		sampleCount := work.sampler.GetMoreSamples(samples, rng)
		for sampleCount > 0 {
			// Generate camera rays and compute radiance along rays
			for i := 0; i < sampleCount; i++ {
				// Find camera ray for _sample[i]_
				var rayWeight float64
				rays[i], rayWeight = work.camera.GenerateRayDifferential(&samples[i])
				rays[i].ScaleDifferentials(1.0 / math.Sqrt(float64(work.sampler.SamplesPerPixel())))
	
				// Evaluate radiance along camera ray
				if rayWeight > 0.0 {
					Ls[i], isects[i], Ts[i] = work.renderer.Li(work.scene, rays[i], &samples[i], rng, nil)
					Ls[i] = Ls[i].Scale(rayWeight)
				} else {
					Ls[i] = NewSpectrum1(0.0)
					isects[i] = nil
					Ts[i] = NewSpectrum1(1.0)
				}
/*
				// Issue warning if unexpected radiance value returned
				if Ls[i].HasNaNs() {
					Error("Not-a-number radiance value returned for image sample.  Setting to black.")
					Ls[i] = NewSpectrum1(0.0)
				} else if Ls[i].Y() < -1.0e-5 {
					Error("Negative luminance value, %f, returned for image sample.  Setting to black.", Ls[i].Y())
					Ls[i] = NewSpectrum1(0.0)
				} else if math.IsInf(Ls[i].Y(), 0) {
					Error("Infinite luminance value returned for image sample.  Setting to black.")
					Ls[i] = NewSpectrum1(0.0)
				}
*/				
			}
			
			// Report sample results to _Sampler_, add contributions to image
			if work.sampler.ReportResults(samples, rays, Ls, isects, sampleCount) {
				for i := 0; i < sampleCount; i++ {
					results <- taskOutput{samples[i], *Ls[i]}
				}
			}

			sampleCount = work.sampler.GetMoreSamples(samples, rng)
		}
		work.reporter.Update(1)
	}
	completed <- true
}
