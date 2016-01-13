package pbrt

type Light interface {
    Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64) (s *Spectrum, wi *Vector, pdf float64, vis *VisibilityTester)
    Power(scene *Scene) *Spectrum
    IsDeltaLight() bool
    Le(ray *RayDifferential) *Spectrum
    Pdf(p *Point, wi *Vector) float64
    Sample_L2(scene *Scene, ls *LightSample, u1,  u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64)
    SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum)
}

type LightData struct {
	nSamples int
	LightToWorld, WorldToLight *Transform
}

type VisibilityTester struct {
    r Ray
}

func (v *VisibilityTester) SetSegment(p1 *Point, eps1 float64, p2 *Point, eps2, time float64) {
    dist := DistancePoint(p1, p2)
    v.r = *CreateRay(p1, p2.Sub(p1).InvScale(dist), eps1, dist * (1.0 - eps2), time, 0)
}
func (v *VisibilityTester) SetRay(p *Point, eps float64, w *Vector, time float64) {
    v.r = *CreateRay(p, w, eps, INFINITY, time, 0)
}
func (v *VisibilityTester) Unoccluded(scene *Scene) bool {
	return !scene.IntersectP(&v.r)
}
func (v *VisibilityTester) Transmittance(scene *Scene, renderer Renderer,
    sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return renderer.Transmittance(scene, CreateRayDifferentialFromRay(&v.r), sample, rng, arena)  	
}

type AreaLight interface {
	Light
	L(p *Point, n *Normal, w *Vector) *Spectrum
}

type LightSample struct {
    uPos [2]float64
   	uComponent float64
}

type LightSampleOffsets struct {
    nSamples, componentOffset, posOffset int
};

func CreateLightSampleOffsets(count int, sample *Sample) *LightSampleOffsets {
	return &LightSampleOffsets{count, sample.Add1D(count), sample.Add2D(count)}
}

func CreateLightSample(sample *Sample, offsets *LightSampleOffsets, n int) *LightSample {
	return &LightSample{[2]float64{sample.twoD[offsets.posOffset][2*n], sample.twoD[offsets.posOffset][2*n+1]}, sample.oneD[offsets.componentOffset][n]}
}

func CreateLightSampleRandom(rng *RNG) *LightSample {
	return &LightSample{[2]float64{rng.RandomFloat(), rng.RandomFloat()}, rng.RandomFloat() }
}

type ShapeSet struct {
    shapes []Shape
    sumArea float64
    areas []float64
    areaDistribution *Distribution1D
}

func CreateShapeSet(s *Shape) *ShapeSet {
	return nil
}
func (s *ShapeSet) Area() float64 { 
	return s.sumArea
}

func (s *ShapeSet) SampleAt(p *Point, ls *LightSample) (*Point, *Normal) {
    sn, _ := s.areaDistribution.SampleDiscrete(ls.uComponent)
    return s.shapes[sn].SampleAt(p, ls.uPos[0], ls.uPos[1])   
}

func (s *ShapeSet) Sample(ls *LightSample) (*Point, *Normal) {
   	sn, _ := s.areaDistribution.SampleDiscrete(ls.uComponent)
    return s.shapes[sn].Sample(ls.uPos[0], ls.uPos[1])
}
func (s *ShapeSet) Pdf2(p *Point, wi *Vector) float64 {
    pdf := 0.0
    for i, shape := range s.shapes {
        pdf += s.areas[i] * shape.Pdf2(p, wi)
	}        
    return pdf / s.sumArea
} 
func (s *ShapeSet) Pdf(p *Point) float64 {
    pdf := 0.0
    for i, shape := range s.shapes {
        pdf += s.areas[i] * shape.Pdf(p)
	}        
    return pdf / (float64(len(s.shapes)) * s.sumArea)
}


