package core

type Renderer interface {
	Render(scene *Scene)
	Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) 
	Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}

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
    
   SurfacePoint struct {
       p Point
       n Normal
       area, rayEpsilon float64
   }
   
   SurfacePointsRenderer struct {
       minDist, time float64
       pCamera Point
       filename string
       points []SurfacePoint
   }
)

func (r *CreateRadianceProbes) Render(scene *Scene) {}
func (r *CreateRadianceProbes) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum ) { return nil, nil, nil }
func (r *CreateRadianceProbes) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func (r *SurfacePointsRenderer) Render(scene *Scene) {}
func (r *SurfacePointsRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) { return nil, nil, nil }
func (r *SurfacePointsRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func CreateRadianceProbesRenderer(camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, params *ParamSet) *CreateRadianceProbes { return nil }	
func CreateSurfacePointsRenderer(params *ParamSet, pCamera *Point, time float64) *SurfacePointsRenderer { return nil }

func FindPoissonPointDistribution(pCamera *Point, time,minDist float64, scene *Scene, points *[]SurfacePoint) {
	Unimplemented()
}