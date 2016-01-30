package core

type (
    CreateRadianceProbes struct {
        surfaceIntegrator SurfaceIntegrator
        volumeIntegrator VolumeIntegrator
        camera Camera
        lmax, nIndirSamples int
        bbox BBox
        includeDirectInProbes, includeIndirectInProbes bool
        time, probeSpacing float64
        filename string       
   }
)

func (r *CreateRadianceProbes) Render(scene *Scene) {}
func (r *CreateRadianceProbes) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum ) { return nil, nil, nil }
func (r *CreateRadianceProbes) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func CreateRadianceProbesRenderer(camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, params *ParamSet) *CreateRadianceProbes { return nil }	
