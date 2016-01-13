package pbrt

type Renderer interface {
	Render(scene *Scene)
    Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum
    Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}