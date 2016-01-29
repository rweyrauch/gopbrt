package core

import (
	"math"
)

type Light interface {
	Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64)
	Power(scene *Scene) *Spectrum
	IsDeltaLight() bool
	Le(ray *RayDifferential) *Spectrum
	Pdf(p *Point, wi *Vector) float64
	Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64)
	SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum)
	NumSamples() int
}

type LightData struct {
	nSamples                   int
	LightToWorld, WorldToLight *Transform
}

func LightSHProject(light Light, p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	coeffs = make([]Spectrum, SHTerms(lmax), SHTerms(lmax))
    for i, _ := range coeffs {
        coeffs[i] = *NewSpectrum1(0.0)
    }
    
    ns := RoundUpPow2(uint32(light.NumSamples()))
    scramble1D := rng.RandomUInt()
    scramble2D := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    Ylm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    var i uint32
    for i = 0; i < ns; i++ {
        // Compute incident radiance sample from _light_, update SH _coeffs_
        u := [2]float64{0.0, 0.0}
        Sample02(i, scramble2D, u[:])
        lightSample := &LightSample{u, VanDerCorput(i, scramble1D)}
        var vis VisibilityTester
        Li, wi, pdf := light.Sample_L(p, pEpsilon, lightSample, time, &vis)
        if !Li.IsBlack() && pdf > 0.0 && (!computeLightVisibility || vis.Unoccluded(scene)) {
            // Add light sample contribution to MC estimate of SH coefficients
            SHEvaluate(wi, lmax, Ylm)
            for j := 0; j < SHTerms(lmax); j++ {
                coeffs[j] = *coeffs[j].Add(Li.Scale(Ylm[j] / (pdf * float64(ns))))
            }    
        }
    }
    return coeffs
}

type VisibilityTester struct {
	r Ray
}

func (v *VisibilityTester) SetSegment(p1 *Point, eps1 float64, p2 *Point, eps2, time float64) {
	dist := DistancePoint(p1, p2)
	v.r = *CreateRay(p1, p2.Sub(p1).InvScale(dist), eps1, dist*(1.0-eps2), time, 0)
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
	uPos       [2]float64
	uComponent float64
}

type LightSampleOffsets struct {
	nSamples, componentOffset, posOffset int
}

func CreateLightSampleOffsets(count int, sample *Sample) *LightSampleOffsets {
	return &LightSampleOffsets{count, sample.Add1D(count), sample.Add2D(count)}
}

func CreateLightSample(sample *Sample, offsets *LightSampleOffsets, n int) *LightSample {
	return &LightSample{[2]float64{sample.twoD[offsets.posOffset][2*n], sample.twoD[offsets.posOffset][2*n+1]}, sample.oneD[offsets.componentOffset][n]}
}

func CreateLightSampleRandom(rng *RNG) *LightSample {
	return &LightSample{[2]float64{rng.RandomFloat(), rng.RandomFloat()}, rng.RandomFloat()}
}

type ShapeSet struct {
	shapes           []Shape
	sumArea          float64
	areas            []float64
	areaDistribution *Distribution1D
}

func NewShapeSet(s Shape) *ShapeSet {
	ss := new(ShapeSet)
	
	ss.shapes = make([]Shape, 0, 4)
	ss.areas = make([]float64, 0, 4)
	
    todo := make([]Shape, 1, 1)
    todo[0] = s
    
    for len(todo) != 0 {
		sh := todo[len(todo)-1]
		todo = todo[:len(todo)-1] // pop off of stack
        if sh.CanIntersect() {
            ss.shapes = append(ss.shapes, sh)
        } else {
            todo = sh.Refine(todo)
		}            
    }
    if len(ss.shapes) > 64 {
        Warning("Area light geometry turned into %d shapes; may be very inefficient.", len(ss.shapes))
	}
    
    // Compute total area of shapes in _ShapeSet_ and area CDF
    ss.sumArea = 0.0
    for _, s := range ss.shapes {
        a := s.Area()
        ss.areas = append(ss.areas, a)
        ss.sumArea += a
    }
    ss.areaDistribution = NewDistribution1D(ss.areas)
    
    return ss
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

type (
	DiffuseAreaLight struct {
		LightData
		Lemit    Spectrum
		shapeSet *ShapeSet
		area     float64
	}

	InfiniteAreaLight struct {
		LightData
		radianceMap  *MIPMapSpectrum
		distribution *Distribution2D
	}

	DistantLight struct {
		LightData
		lightDir Vector
		L        Spectrum
	}

	GonioPhotometricLight struct {
		LightData
		lightPos  Point
		Intensity Spectrum
		mipmap    *MIPMapSpectrum
	}

	PointLight struct {
		LightData
		lightPos  Point
		Intensity Spectrum
	}

	ProjectionLight struct {
		LightData
		projectionMap                          *MIPMapSpectrum
		lightPos                               Point
		Intensity                              Spectrum
		lightProjection                        *Transform
		hither, yon                            float64
		screenX0, screenX1, screenY0, screenY1 float64
		cosTotalWidth                          float64
	}

	SpotLight struct {
		LightData
		lightPos                       Point
		Intensity                      Spectrum
		cosTotalWidth, cosFalloffStart float64
	}
)

func NewDiffuseAreaLight(light2world *Transform, Le *Spectrum, ns int, shape Shape) *DiffuseAreaLight {
	light := new(DiffuseAreaLight)
	light.nSamples = ns
	light.LightToWorld = light2world
	light.WorldToLight = InverseTransform(light2world)
	light.Lemit = *Le
	light.shapeSet = NewShapeSet(shape)
    light.area = light.shapeSet.Area()

	return light
} 

func (l *DiffuseAreaLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (Ls *Spectrum, wi *Vector, pdf float64) {
    //PBRT_AREA_LIGHT_STARTED_SAMPLE()
    ps, ns := l.shapeSet.SampleAt(p, ls)
    wi = NormalizeVector(ps.Sub(p))
    pdf = l.shapeSet.Pdf2(p, wi)
    vis.SetSegment(p, pEpsilon, ps, 1.0e-3, time)
    Ls = l.L(ps, ns, wi.Negate())
    //PBRT_AREA_LIGHT_FINISHED_SAMPLE()
    return Ls, wi, pdf
}

func (l *DiffuseAreaLight) Power(scene *Scene) *Spectrum {
	return l.Lemit.Scale(l.area * math.Pi)
}

func (l *DiffuseAreaLight) IsDeltaLight() bool { return false }

func (l *DiffuseAreaLight) Le(ray *RayDifferential) *Spectrum {
	return NewSpectrum1(0.0)
}

func (l *DiffuseAreaLight) Pdf(p *Point, wi *Vector) float64 {
	return l.shapeSet.Pdf2(p, wi)
}

func (l *DiffuseAreaLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (Ls *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
    //PBRT_AREA_LIGHT_STARTED_SAMPLE();
    org, Ns := l.shapeSet.Sample(ls)
    dir := UniformSampleSphere(u1, u2)
    if DotVectorNormal(dir, Ns) < 0.0  { dir = dir.Negate() }
    ray = CreateRay(org, dir, 1.0e-3, INFINITY, time, 0)
    pdf = l.shapeSet.Pdf(org) * INV_TWOPI
    Ls = l.L(org, Ns, dir)
    //PBRT_AREA_LIGHT_FINISHED_SAMPLE()
    return Ls, ray, Ns, pdf
}
func (l *DiffuseAreaLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)		
}

func (l *DiffuseAreaLight) NumSamples() int {
	return l.nSamples
}

func (l *DiffuseAreaLight) L(p *Point, n *Normal, w *Vector) *Spectrum {
	if DotNormalVector(n, w) > 0.0 {
		return &l.Lemit
	} else {
		return NewSpectrum1(0.0)
	}
}

func CreateDiffuseAreaLight(light2world *Transform, paramSet *ParamSet, shape Shape) AreaLight {
     L := paramSet.FindSpectrumParam("L", *NewSpectrum1(1.0))
     sc := paramSet.FindSpectrumParam("scale", *NewSpectrum1(1.0))
    nSamples := paramSet.FindIntParam("nsamples", 1)
    if options.QuickRender { nSamples = Maxi(1, nSamples / 4) }
    return NewDiffuseAreaLight(light2world, L.Mult(&sc), nSamples, shape)
}

func NewDistantLight(light2world *Transform, radiance *Spectrum, dir *Vector) *DistantLight {
	light := new(DistantLight)
	light.nSamples = 1
	light.LightToWorld = light2world
	light.WorldToLight = InverseTransform(light2world)
	light.lightDir = *NormalizeVector(VectorTransform(light2world, dir))
	light.L = *radiance
	return light
}

func (l *DistantLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (L *Spectrum, wi *Vector, pdf float64) {
	L = &l.L
	wi = &l.lightDir
	pdf = 1.0
	vis.SetRay(p, pEpsilon, wi, time)
	return L, wi, 1.0
}

func (l *DistantLight) Power(scene *Scene) *Spectrum {
	_, worldRadius := scene.bound.BoundingSphere()
	return l.L.Scale(math.Pi * worldRadius * worldRadius)
}

func (l *DistantLight) IsDeltaLight() bool                { return true }
func (l *DistantLight) Le(ray *RayDifferential) *Spectrum { return NewSpectrum1(0.0) }
func (l *DistantLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }

func (l *DistantLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (L *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	worldCenter, worldRadius := scene.bound.BoundingSphere()
	v1, v2 := CoordinateSystem(&l.lightDir)
	d1, d2 := ConcentricSampleDisk(ls.uPos[0], ls.uPos[1])
	Pdisk := worldCenter.Add((v1.Scale(d1).Add(v2.Scale(d2))).Scale(worldRadius))

	// Set ray origin and direction for infinite light ray
	ray = CreateRay(Pdisk.Add(l.lightDir.Scale(worldRadius)), l.lightDir.Negate(), 0.0, INFINITY, time, 0)
	Ns = CreateNormalFromVector(&ray.dir)
	pdf = 1.0 / (math.Pi * worldRadius * worldRadius)
	L = &l.L
	return L, ray, Ns, pdf
}

func (l *DistantLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)	
}

func (l *DistantLight) NumSamples() int {
	return l.nSamples
}

func CreateDistantLight(light2world *Transform, paramSet *ParamSet) *DistantLight {
	L := paramSet.FindSpectrumParam("L", *NewSpectrum1(1.0))
	sc := paramSet.FindSpectrumParam("scale", *NewSpectrum1(1.0))
	from := paramSet.FindPointParam("from", *CreatePoint(0, 0, 0))
	to := paramSet.FindPointParam("to", *CreatePoint(0, 0, 1))
	dir := from.Sub(&to)
	return NewDistantLight(light2world, L.Mult(&sc), dir)
}

func (l *GonioPhotometricLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64) {
	return nil, nil, 0.0
}
func (l *GonioPhotometricLight) Power(scene *Scene) *Spectrum      { return nil }
func (l *GonioPhotometricLight) IsDeltaLight() bool                { return false }
func (l *GonioPhotometricLight) Le(ray *RayDifferential) *Spectrum { return NewSpectrum1(0.0) }
func (l *GonioPhotometricLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *GonioPhotometricLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	return nil, nil, nil, 0.0
}
func (l *GonioPhotometricLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)		
}
func (l *GonioPhotometricLight) NumSamples() int {
	return l.nSamples
}

func CreateGoniometricLight(light2world *Transform, paramSet *ParamSet) *GonioPhotometricLight {
	Unimplemented()
	return nil
}

func NewInfiniteAreaLight(light2world *Transform, L *Spectrum, ns int, texmap string) *InfiniteAreaLight {
	light := new(InfiniteAreaLight)
	light.nSamples = ns
	light.LightToWorld = light2world
	light.WorldToLight = InverseTransform(light2world)

	var width, height int
	var texels []Spectrum
	// Read texel data from _texmap_ into _texels_
	if len(texmap) != 0 {
		texels, width, height = ReadImage(texmap)
		if texels != nil {
			for i := 0; i < width*height; i++ {
				texels[i] = *texels[i].Mult(L)
			}
		}
	}
	if texels == nil {
		width, height = 1, 1
		texels = make([]Spectrum, 1, 1)
		texels[0] = *L
	}
	light.radianceMap = NewMIPMapSpectrum(width, height, texels, false, 8.0, TEXTURE_REPEAT)

	// Initialize sampling PDFs for infinite area light

	// Compute scalar-valued image _img_ from environment map
	filter := 1.0 / float64(Maxi(width, height))
	img := make([]float64, width*height, width*height)
	for v := 0; v < height; v++ {
		vp := float64(v) / float64(height)
		sinTheta := math.Sin(math.Pi * (float64(v) + 0.5) / float64(height))
		for u := 0; u < width; u++ {
			up := float64(u) / float64(width)
			img[u+v*width] = float64(light.radianceMap.Lookup(up, vp, filter).Y())
			img[u+v*width] *= sinTheta
		}
	}

	// Compute sampling distributions for rows and columns of image
	light.distribution = NewDistribution2D(img, width, height)

	return light
}

func (l *InfiniteAreaLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (Ls *Spectrum, wi *Vector, pdf float64) {
	//PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
	// Find $(u,v)$ sample coordinates in infinite light texture
	uv, mapPdf := l.distribution.SampleContinuous(ls.uPos[0], ls.uPos[1])
	if mapPdf == 0.0 {
		return NewSpectrum1(0.0), nil, 0.0
	}

	// Convert infinite light sample point to direction
	theta, phi := uv[1]*math.Pi, uv[0]*2.0*math.Pi
	costheta, sintheta := math.Cos(theta), math.Sin(theta)
	sinphi, cosphi := math.Sin(phi), math.Cos(phi)
	wi = VectorTransform(l.LightToWorld, CreateVector(sintheta*cosphi, sintheta*sinphi, costheta))

	// Compute PDF for sampled infinite light direction
	pdf = mapPdf / (2.0 * math.Pi * math.Pi * sintheta)
	if sintheta == 0.0 {
		pdf = 0.0
	}

	// Return radiance value for infinite light direction
	vis.SetRay(p, pEpsilon, wi, time)
	Ls = l.radianceMap.Lookup(uv[0], uv[1], 0.0)
	//PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
	return Ls, wi, pdf
}

func (l *InfiniteAreaLight) Power(scene *Scene) *Spectrum {
	_, worldRadius := scene.bound.BoundingSphere()
	return l.radianceMap.Lookup(0.5, 0.5, 0.5).Scale(math.Pi * worldRadius * worldRadius)
}

func (l *InfiniteAreaLight) IsDeltaLight() bool { return false }

func (l *InfiniteAreaLight) Le(ray *RayDifferential) *Spectrum {
	wh := NormalizeVector(VectorTransform(l.WorldToLight, &ray.dir))
	s := SphericalPhi(wh) / (2.0 * math.Pi)
	t := SphericalTheta(wh) / math.Pi
	return l.radianceMap.Lookup(s, t, 0.0)
}

func (l *InfiniteAreaLight) Pdf(p *Point, w *Vector) float64 {
	//PBRT_INFINITE_LIGHT_STARTED_PDF();
	wi := VectorTransform(l.WorldToLight, w)
	theta, phi := SphericalTheta(wi), SphericalPhi(wi)
	sintheta := math.Sin(theta)
	if sintheta == 0.0 {
		return 0.0
	}
	pdf := l.distribution.Pdf(phi/(2.0*math.Pi), theta/math.Pi) / (2.0 * math.Pi * math.Pi * sintheta)
	//PBRT_INFINITE_LIGHT_FINISHED_PDF();
	return pdf
}

func (l *InfiniteAreaLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (Ls *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	//PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
	// Compute direction for infinite light sample ray

	// Find $(u,v)$ sample coordinates in infinite light texture
	uv, mapPdf := l.distribution.SampleContinuous(ls.uPos[0], ls.uPos[1])
	if mapPdf == 0.0 {
		return NewSpectrum1(0.0), nil, nil, 0.0
	}

	theta, phi := uv[1]*math.Pi, uv[0]*2.0*math.Pi
	costheta, sintheta := math.Cos(theta), math.Sin(theta)
	sinphi, cosphi := math.Sin(phi), math.Cos(phi)
	d := VectorTransform(l.LightToWorld, CreateVector(sintheta*cosphi, sintheta*sinphi, costheta)).Negate()
	Ns = CreateNormalFromVector(d)

	// Compute origin for infinite light sample ray
	worldCenter, worldRadius := scene.bound.BoundingSphere()
	v1, v2 := CoordinateSystem(d.Negate())
	d1, d2 := ConcentricSampleDisk(u1, u2)
	Pdisk := worldCenter.Add((v1.Scale(d1).Add(v2.Scale(d2))).Scale(worldRadius))
	ray = CreateRay(Pdisk.Add(d.Negate().Scale(worldRadius)), d, 0.0, INFINITY, time, 0)

	// Compute _InfiniteAreaLight_ ray PDF
	directionPdf := mapPdf / (2.0 * math.Pi * math.Pi * sintheta)
	areaPdf := 1.0 / (math.Pi * worldRadius * worldRadius)
	pdf = directionPdf * areaPdf
	if sintheta == 0.0 {
		pdf = 0.0
	}
	Ls = l.radianceMap.Lookup(uv[0], uv[1], 0.0)
	//PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
	return Ls, ray, Ns, pdf
}

func (l *InfiniteAreaLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)	
}

func (l *InfiniteAreaLight) NumSamples() int {
	return l.nSamples
}

func CreateInfiniteLight(light2world *Transform, paramSet *ParamSet) *InfiniteAreaLight {
	L := paramSet.FindSpectrumParam("L", *NewSpectrum1(1.0))
	sc := paramSet.FindSpectrumParam("scale", *NewSpectrum1(1.0))
	texmap := paramSet.FindFilenameParam("mapname", "")
	nSamples := paramSet.FindIntParam("nsamples", 1)
	if options.QuickRender {
		nSamples = Maxi(1, nSamples/4)
	}
	return NewInfiniteAreaLight(light2world, L.Mult(&sc), nSamples, texmap)
}

func NewPointLight(light2world *Transform, intensity *Spectrum) *PointLight {
	light := new(PointLight)
	light.nSamples = 1
	light.LightToWorld = light2world
	light.WorldToLight = InverseTransform(light2world)
	light.lightPos = *PointTransform(light2world, CreatePoint(0, 0, 0))
	light.Intensity = *intensity
	return light
}

func (l *PointLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (L *Spectrum, wi *Vector, pdf float64) {
	vis.SetSegment(p, pEpsilon, &l.lightPos, 0.0, time)
	L = l.Intensity.Scale(DistanceSquaredPoint(&l.lightPos, p))
	wi = NormalizeVector(l.lightPos.Sub(p))
	pdf = 1.0
	return L, wi, pdf
}

func (l *PointLight) Power(scene *Scene) *Spectrum      { return l.Intensity.Scale(4.0 * math.Pi) }
func (l *PointLight) IsDeltaLight() bool                { return true }
func (l *PointLight) Le(ray *RayDifferential) *Spectrum { return NewSpectrum1(0.0) }
func (l *PointLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *PointLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (L *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	L = &l.Intensity
	ray = CreateRay(&l.lightPos, UniformSampleSphere(ls.uPos[0], ls.uPos[1]), 0.0, INFINITY, time, 0)
	Ns = CreateNormalFromVector(&ray.dir)
	pdf = UniformSpherePdf()
	return L, ray, Ns, pdf
}

func (l *PointLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)		
}

func (l *PointLight) NumSamples() int {
	return l.nSamples
}

func CreatePointLight(light2world *Transform, paramSet *ParamSet) *PointLight {
	I := paramSet.FindSpectrumParam("I", *NewSpectrum1(1.0))
	sc := paramSet.FindSpectrumParam("scale", *NewSpectrum1(1.0))
	P := paramSet.FindPointParam("from", *CreatePoint(0, 0, 0))
	l2w := TranslateTransform(CreateVector(P.x, P.y, P.z)).MultTransform(light2world)
	return NewPointLight(l2w, I.Mult(&sc))
}

func (l *ProjectionLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64) {
	Unimplemented()
	return nil, nil, 0.0
}
func (l *ProjectionLight) Power(scene *Scene) *Spectrum      { return nil }
func (l *ProjectionLight) IsDeltaLight() bool                { return false }
func (l *ProjectionLight) Le(ray *RayDifferential) *Spectrum { return nil }
func (l *ProjectionLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *ProjectionLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	Unimplemented()
	return nil, nil, nil, 0.0
}
func (l *ProjectionLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)	
}
func (l *ProjectionLight) NumSamples() int {
	return l.nSamples
}

func CreateProjectionLight(light2world *Transform, paramSet *ParamSet) *ProjectionLight {
	Unimplemented()
	return nil
}

func (l *SpotLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64) {
	Unimplemented()
	return nil, nil, 0.0
}
func (l *SpotLight) Power(scene *Scene) *Spectrum      { return nil }
func (l *SpotLight) IsDeltaLight() bool                { return false }
func (l *SpotLight) Le(ray *RayDifferential) *Spectrum { return NewSpectrum1(0.0) }
func (l *SpotLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *SpotLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	Unimplemented()
	return nil, nil, nil, 0.0
}
func (l *SpotLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG) (coeffs []Spectrum) {
	return LightSHProject(l, p, pEpsilon, lmax, scene, computeLightVisibility, time, rng)	
}

func (l *SpotLight) NumSamples() int {
	return l.nSamples
}

func CreateSpotLight(light2world *Transform, paramSet *ParamSet) *SpotLight {
	Unimplemented()
	return nil
}
