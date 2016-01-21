package core

import()

type (
   SamplerRenderer struct {
        visualizeObjectIds bool
        sampler Sampler
        surfaceIntegrator SurfaceIntegrator
        volumeIntegrator VolumeIntegrator
        camera Camera
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
	if sample != nil {} // temp code to force compiler to allow sample var without use
	
    // Create and launch _SamplerRendererTask_s for rendering image

    // Compute number of _SamplerRendererTask_s to create for rendering
    nPixels := r.camera.Film().XResolution() * r.camera.Film().YResolution()
    nTasks := Maxi(32 * NumSystemCores(), nPixels / (16*16))
    nTasks = int(RoundUpPow2(uint32(nTasks)))
    
    //ProgressReporter reporter(nTasks, "Rendering");
    //vector<Task *> renderTasks;
    for i := 0; i < nTasks; i++ {
    //    renderTasks.push_back(new SamplerRendererTask(scene, this, camera,
    //                                                  reporter, sampler, sample, 
    //                                                  visualizeObjectIds, 
    //                                                  nTasks-1-i, nTasks))
    }                                                  
    //EnqueueTasks(renderTasks);
    //WaitForAllTasks();
    //for (uint32_t i = 0; i < renderTasks.size(); ++i)
    //    delete renderTasks[i];
    //reporter.Done();
    //PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    //delete sample;
    r.camera.Film().WriteImage(1.0)
}
	
func (r *SamplerRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) {
    //Assert(ray.time == sample->time);
    //Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    li = CreateSpectrum1(0.0)
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
