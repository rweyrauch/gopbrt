package core

type Renderer interface {
	Render(scene *Scene)
	Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum
	Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}

type (
    AggregateTest struct {
        nIterations int
        primitives []Primitive
        bboxes []BBox
    }
   
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
   
   MetropolisRenderer struct {
        camera Camera
        bidirectional bool
        nDirectPixelSamples, nPixelSamples, maxDepth int
        largeStepsPerPixel, nBootstrap, maxConsecutiveRejects int
       // directLighting DirectLightingIntegrator
        nTasksFinished int       
   }
   
   SamplerRenderer struct {
        visualizeObjectIds bool
        sampler Sampler
        surfaceIntegrator SurfaceIntegrator
        volumeIntegrator VolumeIntegrator
        camera Camera
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

func (r *AggregateTest) Render(scene *Scene) {}
func (r *AggregateTest) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *AggregateTest) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func (r *CreateRadianceProbes) Render(scene *Scene) {}
func (r *CreateRadianceProbes) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *CreateRadianceProbes) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func (r *MetropolisRenderer) Render(scene *Scene) {}
func (r *MetropolisRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *MetropolisRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func (r *SamplerRenderer) Render(scene *Scene) {}
func (r *SamplerRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *SamplerRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func (r *SurfacePointsRenderer) Render(scene *Scene) {}
func (r *SurfacePointsRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *SurfacePointsRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func CreateAggregateTestRenderer(param *ParamSet, primitives []Primitive) *AggregateTest { return nil }
func CreateRadianceProbesRenderer(camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, params *ParamSet) *CreateRadianceProbes { return nil }
func CreateMetropolisRenderer(params *ParamSet, camera Camera) *MetropolisRenderer { return nil }
func CreateSamplerRenderer(sampler Sampler, camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, visIds bool) *SamplerRenderer { return nil }
func CreateSurfacePointsRenderer(params *ParamSet, pCamera *Point, time float64) *SurfacePointsRenderer { return nil }