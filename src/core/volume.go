package pbrt

type VolumeRegion interface {
	
}

type VolumeIntegrator interface {
	Integrator
    Li(scene *Scene, renderer *Renderer, ray *RayDifferential, sample *Sample, rng *RNG, transmittance *Spectrum, arena *MemoryArena) *Spectrum
    Transmittance(scene *Scene, renderer *Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}