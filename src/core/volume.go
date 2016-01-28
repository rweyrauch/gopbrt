package core

import (
	"fmt"
	"math"
	"strings"
)

// Volume Scattering Local Definitions
type MeasuredSS struct {
	name                   string
	sigma_prime_s, sigma_a [3]float64 // mm^-1
}

var mss []MeasuredSS = []MeasuredSS{
	// From "A Practical Model for Subsurface Light Transport"
	// Jensen, Marschner, Levoy, Hanrahan
	// Proc SIGGRAPH 2001
	{"Apple", [3]float64{2.29, 2.39, 1.97}, [3]float64{0.0030, 0.0034, 0.046}},
	{"Chicken1", [3]float64{0.15, 0.21, 0.38}, [3]float64{0.015, 0.077, 0.19}},
	{"Chicken2", [3]float64{0.19, 0.25, 0.32}, [3]float64{0.018, 0.088, 0.20}},
	{"Cream", [3]float64{7.38, 5.47, 3.15}, [3]float64{0.0002, 0.0028, 0.0163}},
	{"Ketchup", [3]float64{0.18, 0.07, 0.03}, [3]float64{0.061, 0.97, 1.45}},
	{"Marble", [3]float64{2.19, 2.62, 3.00}, [3]float64{0.0021, 0.0041, 0.0071}},
	{"Potato", [3]float64{0.68, 0.70, 0.55}, [3]float64{0.0024, 0.0090, 0.12}},
	{"Skimmilk", [3]float64{0.70, 1.22, 1.90}, [3]float64{0.0014, 0.0025, 0.0142}},
	{"Skin1", [3]float64{0.74, 0.88, 1.01}, [3]float64{0.032, 0.17, 0.48}},
	{"Skin2", [3]float64{1.09, 1.59, 1.79}, [3]float64{0.013, 0.070, 0.145}},
	{"Spectralon", [3]float64{11.6, 20.4, 14.9}, [3]float64{0.00, 0.00, 0.00}},
	{"Wholemilk", [3]float64{2.55, 3.21, 3.77}, [3]float64{0.0011, 0.0024, 0.014}},

	// From "Acquiring Scattering Properties of Participating Media by Dilution",
	// Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
	// Proc SIGGRAPH 2006
	{"Lowfat Milk", [3]float64{0.912600, 1.074800, 1.250000}, [3]float64{0.000200, 0.000400, 0.000800}},
	{"Reduced Milk", [3]float64{1.075000, 1.221300, 1.394100}, [3]float64{0.000200, 0.000400, 0.001000}},
	{"Regular Milk", [3]float64{1.187400, 1.329600, 1.460200}, [3]float64{0.000100, 0.000300, 0.001300}},
	{"Espresso", [3]float64{0.437600, 0.511500, 0.604800}, [3]float64{0.166900, 0.228700, 0.307800}},
	{"Mint Mocha Coffee", [3]float64{0.190000, 0.260000, 0.350000}, [3]float64{0.098400, 0.151900, 0.204000}},
	{"Lowfat Soy Milk", [3]float64{0.141900, 0.162500, 0.274000}, [3]float64{0.000100, 0.000500, 0.002500}},
	{"Regular Soy Milk", [3]float64{0.243400, 0.271900, 0.459700}, [3]float64{0.000100, 0.000500, 0.003400}},
	{"Lowfat Chocolate Milk", [3]float64{0.428200, 0.501400, 0.579100}, [3]float64{0.000500, 0.001600, 0.006800}},
	{"Regular Chocolate Milk", [3]float64{0.735900, 0.917200, 1.068800}, [3]float64{0.000700, 0.003000, 0.010000}},
	{"Coke", [3]float64{0.714300, 1.168800, 1.716900}, [3]float64{0.696600, 1.148000, 1.716900}},
	{"Pepsi", [3]float64{0.643300, 0.999000, 1.442000}, [3]float64{0.637500, 0.984900, 1.442000}},
	{"Sprite", [3]float64{0.129900, 0.128300, 0.139500}, [3]float64{0.123000, 0.119400, 0.130600}},
	{"Gatorade", [3]float64{0.400900, 0.418500, 0.432400}, [3]float64{0.161700, 0.125800, 0.057900}},
	{"Chardonnay", [3]float64{0.157700, 0.174800, 0.351200}, [3]float64{0.154700, 0.170100, 0.344300}},
	{"White Zinfandel", [3]float64{0.176300, 0.237000, 0.291300}, [3]float64{0.173200, 0.232200, 0.284700}},
	{"Merlot", [3]float64{0.763900, 1.642900, 1.919600}, [3]float64{0.758600, 1.642900, 1.919600}},
	{"Budweiser Beer", [3]float64{0.148600, 0.321000, 0.736000}, [3]float64{0.144900, 0.314100, 0.728600}},
	{"Coors Light Beer", [3]float64{0.029500, 0.066300, 0.152100}, [3]float64{0.026800, 0.060800, 0.152100}},
	{"Clorox", [3]float64{0.160000, 0.250000, 0.330000}, [3]float64{0.017500, 0.077700, 0.137200}},
	{"Apple Juice", [3]float64{0.121500, 0.210100, 0.440700}, [3]float64{0.101400, 0.185800, 0.408400}},
	{"Cranberry Juice", [3]float64{0.270000, 0.630000, 0.830000}, [3]float64{0.257200, 0.614500, 0.810400}},
	{"Grape Juice", [3]float64{0.550000, 1.250000, 1.530000}, [3]float64{0.542800, 1.250000, 1.530000}},
	{"Ruby Grapefruit Juice", [3]float64{0.251300, 0.351700, 0.430500}, [3]float64{0.089600, 0.191100, 0.263600}},
	{"White Grapefruit Juice", [3]float64{0.360900, 0.380000, 0.563200}, [3]float64{0.009600, 0.013100, 0.039500}},
	{"Shampoo", [3]float64{0.028800, 0.071000, 0.095200}, [3]float64{0.018400, 0.059600, 0.080500}},
	{"Strawberry Shampoo", [3]float64{0.021700, 0.078800, 0.102200}, [3]float64{0.018900, 0.075600, 0.098900}},
	{"Head & Shoulders Shampoo", [3]float64{0.367400, 0.452700, 0.521100}, [3]float64{0.088300, 0.163700, 0.212500}},
	{"Lemon Tea", [3]float64{0.340000, 0.580000, 0.880000}, [3]float64{0.260200, 0.490200, 0.772700}},
	{"Orange Juice Powder", [3]float64{0.337700, 0.557300, 1.012200}, [3]float64{0.144900, 0.344100, 0.786300}},
	{"Pink Lemonade", [3]float64{0.240000, 0.370000, 0.450000}, [3]float64{0.116500, 0.236600, 0.319500}},
	{"Cappuccino Powder", [3]float64{0.257400, 0.353600, 0.484000}, [3]float64{0.192000, 0.265400, 0.327200}},
	{"Salt Powder", [3]float64{0.760000, 0.868500, 0.936300}, [3]float64{0.511500, 0.586300, 0.614700}},
	{"Sugar Powder", [3]float64{0.079500, 0.175900, 0.278000}, [3]float64{0.065000, 0.159700, 0.257800}},
	{"Suisse Mocha", [3]float64{0.509800, 0.647600, 0.794400}, [3]float64{0.187500, 0.289300, 0.379600}},
	{"Pacific Ocean Surface Water", [3]float64{3.364500, 3.315800, 3.242800}, [3]float64{3.184500, 3.132400, 3.014700}},
}

func RdIntegral(alphap, A float64) float64 {
	sqrtTerm := math.Sqrt(3.0 * (1.0 - alphap))
	return alphap / 2.0 * (1.0 + math.Exp(-4.0/3.0*A*sqrtTerm)) * math.Exp(-sqrtTerm)
}

func RdToAlphap(reflectance, A float64) float64 {
	alphaLow := 0.0
	alphaHigh := 1.0
	kd0 := RdIntegral(alphaLow, A)
	kd1 := RdIntegral(alphaHigh, A)
	for i := 0; i < 16; i++ {
		Assert(kd0 <= reflectance && kd1 >= reflectance)
		alphaMid := (alphaLow + alphaHigh) * 0.5
		kd := RdIntegral(alphaMid, A)
		if kd < reflectance {
			alphaLow = alphaMid
			kd0 = kd
		} else {
			alphaHigh = alphaMid
			kd1 = kd
		}
	}
	return (alphaLow + alphaHigh) * 0.5
}

type VolumeRegion interface {
	WorldBound() *BBox
	IntersectP(ray *Ray) (hit bool, t0, t1 float64)
	Sigma_a(p *Point, w *Vector, time float64) *Spectrum
	Sigma_s(p *Point, w *Vector, time float64) *Spectrum
	Lve(p *Point, w *Vector, time float64) *Spectrum
	P(p *Point, w, wp *Vector, time float64) float64
	Sigma_t(p *Point, wo *Vector, time float64) *Spectrum
	Tau(ray *Ray, step, offset float64) *Spectrum
}

type DensityRegion interface {
	VolumeRegion
	Density(Pobj *Point) float64
}

type DensityRegionData struct {
	sig_a, sig_s, le *Spectrum
	g                float64
	worldToVolume    *Transform
}

type ExponentialDensity struct {
	DensityRegionData
	extent *BBox
	a, b   float64
	upDir  *Vector
}

func (r *ExponentialDensity) WorldBound() *BBox { return nil }
func (r *ExponentialDensity) IntersectP(ray *Ray) (hit bool, t0, t1 float64) {
	Unimplemented()
	return false, 0.0, 0.0
}
func (r *ExponentialDensity) Sigma_a(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *ExponentialDensity) Sigma_s(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *ExponentialDensity) Lve(p *Point, w *Vector, time float64) *Spectrum      { return nil }
func (r *ExponentialDensity) P(p *Point, w, wp *Vector, time float64) float64      { return 0.0 }
func (r *ExponentialDensity) Sigma_t(p *Point, wo *Vector, time float64) *Spectrum { return nil }
func (r *ExponentialDensity) Tau(ray *Ray, step, offset float64) *Spectrum         { return nil }

func (r *ExponentialDensity) Density(Pobj *Point) float64 {
	if !r.extent.Inside(Pobj) {
		return 0.0
	}
	height := DotVector(Pobj.Sub(&r.extent.pMin), r.upDir)
	return r.a * math.Exp(-r.b*height)
}

func CreateExponentialVolumeRegion(volume2world *Transform, params *ParamSet) *ExponentialDensity {
	Unimplemented()
	return nil
}

type HomogeneousVolumeDensity struct {
	DensityRegionData
	extent *BBox
}

func (r *HomogeneousVolumeDensity) WorldBound() *BBox { return nil }
func (r *HomogeneousVolumeDensity) IntersectP(ray *Ray) (hit bool, t0, t1 float64) {
	Unimplemented()
	return false, 0.0, 0.0
}
func (r *HomogeneousVolumeDensity) Sigma_a(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *HomogeneousVolumeDensity) Sigma_s(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *HomogeneousVolumeDensity) Lve(p *Point, w *Vector, time float64) *Spectrum      { return nil }
func (r *HomogeneousVolumeDensity) P(p *Point, w, wp *Vector, time float64) float64      { return 0.0 }
func (r *HomogeneousVolumeDensity) Sigma_t(p *Point, wo *Vector, time float64) *Spectrum { return nil }
func (r *HomogeneousVolumeDensity) Tau(ray *Ray, step, offset float64) *Spectrum         { return nil }

func CreateHomogeneousVolumeDensityRegion(volume2world *Transform, params *ParamSet) *HomogeneousVolumeDensity {
	Unimplemented()
	return nil
}

type VolumeGridDensity struct {
	DensityRegionData
	density    []float64
	nx, ny, nz int
	extent     *BBox
}

func (r *VolumeGridDensity) WorldBound() *BBox { return nil }
func (r *VolumeGridDensity) IntersectP(ray *Ray) (hit bool, t0, t1 float64) {
	return false, 0.0, 0.0
}
func (r *VolumeGridDensity) Sigma_a(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *VolumeGridDensity) Sigma_s(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *VolumeGridDensity) Lve(p *Point, w *Vector, time float64) *Spectrum      { return nil }
func (r *VolumeGridDensity) P(p *Point, w, wp *Vector, time float64) float64      { return 0.0 }
func (r *VolumeGridDensity) Sigma_t(p *Point, wo *Vector, time float64) *Spectrum { return nil }
func (r *VolumeGridDensity) Tau(ray *Ray, step, offset float64) *Spectrum         { return nil }
func (r *VolumeGridDensity) Density(Pobj *Point) float64                          { return 0.0 }

func CreateGridVolumeRegion(volume2world *Transform, params *ParamSet) *VolumeGridDensity {
	Unimplemented()
	return nil
}

type AggregateVolume struct {
	regions []VolumeRegion
	bound   *BBox
}

func (v *AggregateVolume) WorldBound() *BBox                                    { return nil }
func (v *AggregateVolume) IntersectP(ray *Ray) (hit bool, t0, t1 float64)       { return false, 0.0, 0.0 }
func (v *AggregateVolume) Sigma_a(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (v *AggregateVolume) Sigma_s(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (v *AggregateVolume) Lve(p *Point, w *Vector, time float64) *Spectrum      { return nil }
func (v *AggregateVolume) P(p *Point, w, wp *Vector, time float64) float64      { return 0.0 }
func (v *AggregateVolume) Sigma_t(p *Point, wo *Vector, time float64) *Spectrum { return nil }
func (v *AggregateVolume) Tau(ray *Ray, step, offset float64) *Spectrum         { return nil }

func CreateAggregateVolume(regions []VolumeRegion) *AggregateVolume {
	Unimplemented()
	return nil
}

type VolumeIntegrator interface {
	Integrator
	Li(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li, transmittance *Spectrum)
	Transmittance(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}

type (
	EmissionIntegrator struct {
		stepSize                             float64
		tauSampleOffset, scatterSampleOffset int
	}

	SingleScatteringIntegrator struct {
		stepSize                             float64
		tauSampleOffset, scatterSampleOffset int
	}
)

func NewEmissionIntegrator(stepsize float64) *EmissionIntegrator {
	return &EmissionIntegrator{stepsize, 0, 0}
}

func (integrator *EmissionIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer) {}
func (integrator *EmissionIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
	integrator.tauSampleOffset = sample.Add1D(1)
	integrator.scatterSampleOffset = sample.Add1D(1)
}

func (integrator *EmissionIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li, transmittance *Spectrum) {
	Assert(sample != nil)
	var hit bool
	var t0, t1 float64
	if scene.volumeRegion != nil {
		hit, t0, t1 = scene.volumeRegion.IntersectP(CreateRayFromRayDifferential(ray))
		if !hit || (t1-t0 == 0.0) {
			transmittance = NewSpectrum1(1.0)
			li = NewSpectrum1(0.0)
			return li, transmittance
		}
	} else {
		transmittance = NewSpectrum1(1.0)
		li = NewSpectrum1(0.0)
		return li, transmittance
	}

	// Do emission-only volume integration in _vr_
	Lv := NewSpectrum1(0.0)

	// Prepare for volume integration stepping
	nSamples := Ceil2Int((t1 - t0) / integrator.stepSize)

	step := (t1 - t0) / float64(nSamples)
	Tr := NewSpectrum1(1.0)
	p := ray.PointAt(t0)
	var pPrev *Point
	w := ray.dir.Negate()
	t0 += sample.oneD[integrator.scatterSampleOffset][0] * step
	for i := 0; i < nSamples; i++ {
		// Advance to sample at _t0_ and update _T_
		pPrev = p
		p = ray.PointAt(t0)
		tauRay := CreateRay(pPrev, p.Sub(pPrev), 0.0, 1.0, ray.time, ray.depth)
		stepTau := scene.volumeRegion.Tau(tauRay, 0.5*integrator.stepSize, rng.RandomFloat())
		Tr = Tr.Mult(ExpSpectrum(stepTau.Negate()))

		// Possibly terminate ray marching if transmittance is small
		if Tr.Y() < 1.0e-3 {
			continueProb := 0.5
			if rng.RandomFloat() > continueProb {
				Tr = NewSpectrum1(0.0)
				break
			}
			Tr = Tr.Scale(continueProb)
		}

		// Compute emission-only source term at _p_
		Lv = Lv.Add(Tr.Mult(scene.volumeRegion.Lve(p, w, ray.time)))
		t0 += step
	}
	transmittance = Tr
	li = Lv.Scale(step)
	return li, transmittance
}

func (integrator *EmissionIntegrator) Transmittance(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	if scene.volumeRegion == nil {
		return NewSpectrum1(1.0)
	}
	var step, offset float64
	if sample != nil {
		step = integrator.stepSize
		offset = sample.oneD[integrator.tauSampleOffset][0]
	} else {
		step = 4.0 * integrator.stepSize
		offset = rng.RandomFloat()
	}
	tau := scene.volumeRegion.Tau(CreateRayFromRayDifferential(ray), step, offset)
	return ExpSpectrum(tau.Scale(-1.0))
}

func (integrator *EmissionIntegrator) String() string {
	return fmt.Sprintf("emmission[step: %f tauoff: %d scatteroff: %d]", integrator.stepSize, integrator.tauSampleOffset, integrator.scatterSampleOffset)
}

func (i *SingleScatteringIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (i *SingleScatteringIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (i *SingleScatteringIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li, transmittance *Spectrum) {
	return NewSpectrum1(0.0), NewSpectrum1(1.0)
}
func (i *SingleScatteringIntegrator) Transmittance(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return NewSpectrum1(1.0)
}

func CreateSingleScatteringIntegrator(params *ParamSet) *SingleScatteringIntegrator {
	Unimplemented()
	return nil
}
func CreateEmissionVolumeIntegrator(params *ParamSet) *EmissionIntegrator {
	stepSize := params.FindFloatParam("stepsize", 1.0)
	return NewEmissionIntegrator(stepSize)
}

// Volume Scattering Declarations
func PhaseIsotropic(w, wp *Vector) float64 {
	return 1.0 / (4.0 * math.Pi)
}
func PhaseRayleigh(w, wp *Vector) float64 {
	costheta := DotVector(w, wp)
	return 3.0 / (16.0 * math.Pi) * (1 + costheta*costheta)
}
func PhaseMieHazy(w, wp *Vector) float64 {
	costheta := DotVector(w, wp)
	return (0.5 + 4.5*math.Pow(0.5*(1.0+costheta), 8.0)) / (4.0 * math.Pi)
}
func PhaseMieMurky(w, wp *Vector) float64 {
	costheta := DotVector(w, wp)
	return (0.5 + 16.5*math.Pow(0.5*(1.0+costheta), 32.0)) / (4.0 * math.Pi)
}
func PhaseHG(w, wp *Vector, g float64) float64 {
	costheta := DotVector(w, wp)
	return 1.0 / (4.0 * math.Pi) * (1.0 - g*g) / math.Pow(1.0+g*g-2.0*g*costheta, 1.5)
}
func PhaseSchlick(w, wp *Vector, g float64) float64 {
	// improved g->k mapping derived by Thies Heidecke
	// see http://pbrt.org/bugtracker/view.php?id=102
	alpha := 1.5
	k := alpha*g + (1.0-alpha)*g*g*g
	kcostheta := k * DotVector(w, wp)
	return 1.0 / (4.0 * math.Pi) * (1.0 - k*k) / ((1.0 - kcostheta) * (1.0 - kcostheta))
}

func GetVolumeScatteringProperties(name string) (found bool, sigma_a, sigma_prime_s *Spectrum) {
	for _, m := range mss {
		if strings.Compare(name, m.name) == 0 {
			sigma_a = NewSpectrum(m.sigma_a)
			sigma_prime_s = NewSpectrum(m.sigma_prime_s)
			return true, sigma_a, sigma_prime_s
		}
	}
	return false, nil, nil
}

func SubsurfaceFromDiffuse(Kd *Spectrum, meanPathLength, eta float64) (sigma_a, sigma_prime_s *Spectrum) {
	A := (1.0 + Fdr(eta)) / (1.0 - Fdr(eta))
	var rgb [3]float64
	rgb[0], rgb[1], rgb[2] = Kd.c[0], Kd.c[1], Kd.c[2]
	var sigma_prime_s_rgb, sigma_a_rgb [3]float64
	for i := 0; i < 3; i++ {
		// Compute $\alpha'$ for RGB component, compute scattering properties
		alphap := RdToAlphap(rgb[i], A)
		sigma_tr := 1.0 / meanPathLength
		sigma_prime_t := sigma_tr / math.Sqrt(3.0*(1.0-alphap))
		sigma_prime_s_rgb[i] = alphap * sigma_prime_t
		sigma_a_rgb[i] = sigma_prime_t - sigma_prime_s_rgb[i]
	}
	sigma_a = NewSpectrum(sigma_a_rgb)
	sigma_prime_s = NewSpectrum(sigma_prime_s_rgb)
	return sigma_a, sigma_prime_s
}
