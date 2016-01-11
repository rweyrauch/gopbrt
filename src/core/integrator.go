package pbrt

type Integrator interface {
    Preprocess(scene *Scene, camera *Camera, enderer *Renderer)
    RequestSamples(sampler *Sampler, sample *Sample, scene *Scene)
}

type SurfaceIntegrator interface {
	Integrator
    Li(scene *Scene, renderer *Renderer, ray *RayDifferential, isect *Intersection,
        sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}

