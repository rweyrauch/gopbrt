package core

type (
    
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
	
	surfacePointTask struct {
		taskNum int
		scene *Scene
		origin Point
		time float64
		minSampleDist float64
		maxFails int
		reporter *ProgressReporter
		
	}
)

func NewSurfacePointsRenderer(minDist float64, pCamera *Point, time float64, filename string) *SurfacePointsRenderer {
	renderer := new(SurfacePointsRenderer)
	
	return renderer
}

func (renderer *SurfacePointsRenderer) Render(scene *Scene) {}
func (renderer *SurfacePointsRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) { return nil, nil, nil }
func (renderer *SurfacePointsRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func CreateSurfacePointsRenderer(params *ParamSet, pCamera *Point, time float64) *SurfacePointsRenderer { 
     minDist := params.FindFloatParam("minsampledistance", 0.25)
     filename := params.FindFilenameParam("filename", "")
    if options.QuickRender { minDist *= 4.0 }
    return NewSurfacePointsRenderer(minDist, pCamera, time, filename)
}

func FindPoissonPointDistribution(pCamera *Point, time,minDist float64, scene *Scene, points *[]SurfacePoint) {
	Unimplemented()
}