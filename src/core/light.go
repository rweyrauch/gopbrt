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
	SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum)
}

type LightData struct {
	nSamples                   int
	LightToWorld, WorldToLight *Transform
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

type (
	DiffuseAreaLight struct {
		LightData
		Lemit    Spectrum
		shapeSet []Spectrum
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

func CreateDiffuseAreaLight(light2world *Transform, paramSet *ParamSet, shape Shape) AreaLight {
	return nil
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
	return l.L.Scale(float32(math.Pi * worldRadius * worldRadius))
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

func (l *DistantLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
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
func (l *GonioPhotometricLight) Le(ray *RayDifferential) *Spectrum { return nil }
func (l *GonioPhotometricLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *GonioPhotometricLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	return nil, nil, nil, 0.0
}
func (l *GonioPhotometricLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
}

func CreateGoniometricLight(light2world *Transform, paramSet *ParamSet) *GonioPhotometricLight {
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
	return l.radianceMap.Lookup(0.5, 0.5, 0.5).Scale(float32(math.Pi * worldRadius * worldRadius))
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

func (l *InfiniteAreaLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
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
	L = l.Intensity.Scale(float32(DistanceSquaredPoint(&l.lightPos, p)))
	wi = NormalizeVector(l.lightPos.Sub(p))
	pdf = 1.0
	return L, wi, pdf
}

func (l *PointLight) Power(scene *Scene) *Spectrum      { return l.Intensity.Scale(float32(4.0 * math.Pi)) }
func (l *PointLight) IsDeltaLight() bool                { return true }
func (l *PointLight) Le(ray *RayDifferential) *Spectrum { return nil }
func (l *PointLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *PointLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (L *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	L = &l.Intensity
	ray = CreateRay(&l.lightPos, UniformSampleSphere(ls.uPos[0], ls.uPos[1]), 0.0, INFINITY, time, 0)
	Ns = CreateNormalFromVector(&ray.dir)
	pdf = UniformSpherePdf()
	return L, ray, Ns, pdf
}

func (l *PointLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
}

func CreatePointLight(light2world *Transform, paramSet *ParamSet) *PointLight {
	I := paramSet.FindSpectrumParam("I", *NewSpectrum1(1.0))
	sc := paramSet.FindSpectrumParam("scale", *NewSpectrum1(1.0))
	P := paramSet.FindPointParam("from", *CreatePoint(0, 0, 0))
	l2w := TranslateTransform(CreateVector(P.x, P.y, P.z)).MultTransform(light2world)
	return NewPointLight(l2w, I.Mult(&sc))
}

func (l *ProjectionLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64) {
	return nil, nil, 0.0
}
func (l *ProjectionLight) Power(scene *Scene) *Spectrum      { return nil }
func (l *ProjectionLight) IsDeltaLight() bool                { return false }
func (l *ProjectionLight) Le(ray *RayDifferential) *Spectrum { return nil }
func (l *ProjectionLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *ProjectionLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	return nil, nil, nil, 0.0
}
func (l *ProjectionLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
}

func CreateProjectionLight(light2world *Transform, paramSet *ParamSet) *ProjectionLight {
	return nil
}

func (l *SpotLight) Sample_L(p *Point, pEpsilon float64, ls *LightSample, time float64, vis *VisibilityTester) (s *Spectrum, wi *Vector, pdf float64) {
	return nil, nil, 0.0
}
func (l *SpotLight) Power(scene *Scene) *Spectrum      { return nil }
func (l *SpotLight) IsDeltaLight() bool                { return false }
func (l *SpotLight) Le(ray *RayDifferential) *Spectrum { return nil }
func (l *SpotLight) Pdf(p *Point, wi *Vector) float64  { return 0.0 }
func (l *SpotLight) Sample_L2(scene *Scene, ls *LightSample, u1, u2, time float64) (s *Spectrum, ray *Ray, Ns *Normal, pdf float64) {
	return nil, nil, nil, 0.0
}
func (l *SpotLight) SHProject(p *Point, pEpsilon float64, lmax int, scene *Scene, computeLightVisibility bool, time float64, rng *RNG, coeffs *Spectrum) {
}

func CreateSpotLight(light2world *Transform, paramSet *ParamSet) *SpotLight {
	return nil
}
