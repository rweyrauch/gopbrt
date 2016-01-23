package core

import (
	"math"
)

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
	return nil
}

type HomogeneousVolumeDensity struct {
	DensityRegionData
	extent *BBox
}

func (r *HomogeneousVolumeDensity) WorldBound() *BBox { return nil }
func (r *HomogeneousVolumeDensity) IntersectP(ray *Ray) (hit bool, t0, t1 float64) {
	return false, 0.0, 0.0
}
func (r *HomogeneousVolumeDensity) Sigma_a(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *HomogeneousVolumeDensity) Sigma_s(p *Point, w *Vector, time float64) *Spectrum  { return nil }
func (r *HomogeneousVolumeDensity) Lve(p *Point, w *Vector, time float64) *Spectrum      { return nil }
func (r *HomogeneousVolumeDensity) P(p *Point, w, wp *Vector, time float64) float64      { return 0.0 }
func (r *HomogeneousVolumeDensity) Sigma_t(p *Point, wo *Vector, time float64) *Spectrum { return nil }
func (r *HomogeneousVolumeDensity) Tau(ray *Ray, step, offset float64) *Spectrum         { return nil }

func CreateHomogeneousVolumeDensityRegion(volume2world *Transform, params *ParamSet) *HomogeneousVolumeDensity {
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

func CreateGridVolumeRegion(volume2world *Transform, params *ParamSet) *VolumeGridDensity { return nil }

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

func CreateAggregateVolume(regions []VolumeRegion) *AggregateVolume { return nil }

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

func (i *EmissionIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (i *EmissionIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (i *EmissionIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li, transmittance *Spectrum) {
	return NewSpectrum1(0.0), NewSpectrum1(1.0)
}
func (i *EmissionIntegrator) Transmittance(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return NewSpectrum1(1.0)
}

func (i *SingleScatteringIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (i *SingleScatteringIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (i *SingleScatteringIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li, transmittance *Spectrum) {
	return NewSpectrum1(0.0), NewSpectrum1(1.0)
}
func (i *SingleScatteringIntegrator) Transmittance(scene *Scene, renderer Renderer, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return NewSpectrum1(1.0)
}

func CreateSingleScatteringIntegrator(params *ParamSet) *SingleScatteringIntegrator { return nil }
func CreateEmissionVolumeIntegrator(params *ParamSet) *EmissionIntegrator           { return nil }

// Volume Scattering Declarations
func PhaseIsotropic(w, wp *Vector) float64          { return 0.0 }
func PhaseRayleigh(w, wp *Vector) float64           { return 0.0 }
func PhaseMieHazy(w, wp *Vector) float64            { return 0.0 }
func PhaseMieMurky(w, wp *Vector) float64           { return 0.0 }
func PhaseHG(w, wp *Vector, g float64) float64      { return 0.0 }
func PhaseSchlick(w, wp *Vector, g float64) float64 { return 0.0 }

func GetVolumeScatteringProperties(name string) (found bool, sigma_a, sigma_prime_s *Spectrum) {
	return false, nil, nil
}
func SubsurfaceFromDiffuse(Kd *Spectrum, meanPathLength, eta float64) (sigma_a, sigma_prime_s *Spectrum) {
	return nil, nil
}
