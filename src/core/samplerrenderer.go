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
	"encoding/binary"
	"hash/adler32"
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
	// Allow integrators to do preprocessing for the scene
	r.surfaceIntegrator.Preprocess(scene, r.camera, r)
	r.volumeIntegrator.Preprocess(scene, r.camera, r)
	// Allocate and initialize _sample_
	sample := NewSample(r.sampler, r.surfaceIntegrator, r.volumeIntegrator, scene)

	// Create and launch _SamplerRendererTask_s for rendering image

	// Compute number of _SamplerRendererTask_s to create for rendering
	nPixels := r.camera.Film().XResolution() * r.camera.Film().YResolution()
	nWorkers := Maxi(NumSystemCores(), 1)
	nTasks := Maxi(32*NumSystemCores(), nPixels/(16*16))
	nTasks = int(RoundUpPow2(uint32(nTasks)))

	jobs := make(chan *samplerRendererTask, nTasks)
	producedSamples := make(chan taskOutput, 128)
	completed := make(chan bool, nTasks)

	for w := 0; w < nWorkers; w++ {
		go samplerRendererWorker(jobs, producedSamples, completed)
	}

	reporter := NewProgressReporter(nTasks, "Rendering", -1)
	for i := 0; i < nTasks; i++ {
		sampler := r.sampler.GetSubSampler(i, nTasks)
		task := newSamplerRendererTask(scene, r, r.camera, reporter, sampler, sample, r.visualizeObjectIds, nTasks-1-i, nTasks)
		jobs <- task
	}
	close(jobs)
	numCompleted := 0
	for numCompleted < nWorkers {
		select {
		case output := <-producedSamples:
			r.camera.Film().AddSample(&output.sample, &output.Li)
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
		r.camera.Film().AddSample(&output.sample, &output.Li)
	}

	reporter.Done()

	// Clean up after rendering and store final image
	r.camera.Film().WriteImage(1.0)
}

func (r *SamplerRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (*Spectrum, *Intersection, *Spectrum) {
	Assert(ray.Time() == sample.time)
	Assert(!ray.HasNaNs())
	// Allocate local variables for _isect_ and _T_ if needed
	li := NewSpectrum1(0.0)
	var isect *Intersection

	var hit bool
	if hit, isect = scene.Intersect(ray); hit {
		li = r.surfaceIntegrator.Li(scene, r, ray, isect, sample, rng, arena)
	} else {
		// Handle ray that doesn't intersect any geometry
		for i := 0; i < len(scene.lights); i++ {
			li = li.Add(scene.lights[i].Le(ray))
		}
	}

	Lvi, T := r.volumeIntegrator.Li(scene, r, ray, sample, rng, arena)
	li = li.Mult(T).Add(Lvi)
	return li, isect, T
}

func (r *SamplerRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return r.volumeIntegrator.Transmittance(scene, r, ray, sample, rng, arena)
}

func CreateSamplerRenderer(sampler Sampler, camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, visIds bool) *SamplerRenderer {
	return nil
}

type samplerRendererTask struct {
	scene              *Scene
	renderer           Renderer
	camera             Camera
	sampler            Sampler
	reporter           *ProgressReporter
	origSample         *Sample
	visualizeObjectIds bool
	taskNum, taskCount int
}
type taskOutput struct {
	sample Sample
	Li     Spectrum
}

func newSamplerRendererTask(scene *Scene, renderer Renderer, camera Camera,
	progress *ProgressReporter, sampler Sampler, sample *Sample, visIds bool, taskNum, taskCount int) *samplerRendererTask {
	task := &samplerRendererTask{scene, renderer, camera, sampler, progress, sample, visIds, taskNum, taskCount}
	return task
}

func samplerRendererWorker(workQueue <-chan *samplerRendererTask, results chan<- taskOutput, completed chan<- bool) {

	for t := range workQueue {
		if t.sampler == nil {
			t.reporter.Update(1)
			continue
		}

		// Declare local variables used for rendering loop
		rng := NewRNG(int64(t.taskNum))
		hash := adler32.New()

		// Allocate space for samples and intersections
		maxSamples := t.sampler.MaximumSampleCount()
		samples := t.origSample.Duplicate(maxSamples)
		rays := make([]*RayDifferential, maxSamples, maxSamples)
		Ls := make([]*Spectrum, maxSamples, maxSamples)
		Ts := make([]*Spectrum, maxSamples, maxSamples)
		isects := make([]*Intersection, maxSamples, maxSamples)

		// Get samples from _Sampler_ and update image
		sampleCount := t.sampler.GetMoreSamples(&samples, rng)
		for sampleCount > 0 {
			// Generate camera rays and compute radiance along rays
			for i := 0; i < sampleCount; i++ {
				// Find camera ray for _sample[i]_
				var rayWeight float64
				rays[i], rayWeight = t.camera.GenerateRayDifferential(&samples[i])
				rays[i].ScaleDifferentials(1.0 / math.Sqrt(float64(t.sampler.SamplesPerPixel())))

				// Evaluate radiance along camera ray
				if t.visualizeObjectIds {
					if rayWeight > 0.0 {
						var hit bool
						hit, isects[i] = t.scene.Intersect(rays[i])
						if hit {
							// random shading based on shape id...
							binary.Write(hash, binary.BigEndian, isects[i].shapeId)
							binary.Write(hash, binary.BigEndian, isects[i].primitiveId)
							h := hash.Sum32()
							rgb := [3]float64{float64(h & 0xff), float64((h >> 8) & 0xff), float64((h >> 16) & 0xff)}
							Ls[i] = NewSpectrumRGB(rgb[0], rgb[1], rgb[2])
							Ls[i].Scale(1.0 / 255.0)
						} else {
							Ls[i] = NewSpectrum1(0.0)
						}
					} else {
						Ls[i] = NewSpectrum1(0.0)
					}
					Ts[i] = NewSpectrum1(1.0)
				} else {
					if rayWeight > 0.0 {
						Ls[i], isects[i], Ts[i] = t.renderer.Li(t.scene, rays[i], &samples[i], rng, nil)
						Ls[i] = Ls[i].Scale(rayWeight)
					} else {
						Ls[i] = NewSpectrum1(0.0)
						isects[i] = nil
						Ts[i] = NewSpectrum1(1.0)
					}

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
				}
			}

			// Report sample results to _Sampler_, add contributions to image
			if t.sampler.ReportResults(samples, rays, Ls, isects, sampleCount) {
				for i := 0; i < sampleCount; i++ {
					results <- taskOutput{samples[i], *Ls[i]}
				}
			}

			sampleCount = t.sampler.GetMoreSamples(&samples, rng)
		}
		t.reporter.Update(1)
	}
	completed <- true
}
