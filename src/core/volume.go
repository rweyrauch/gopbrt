package core

import (
	"math"
	"fmt"
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
			Tr = Tr.Scale(float32(continueProb))
		}

		// Compute emission-only source term at _p_
		Lv = Lv.Add(Tr.Mult(scene.volumeRegion.Lve(p, w, ray.time)))
		t0 += step
	}
	transmittance = Tr
	li = Lv.Scale(float32(step))
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
func PhaseIsotropic(w, wp *Vector) float64          { return 0.0 }
func PhaseRayleigh(w, wp *Vector) float64           { return 0.0 }
func PhaseMieHazy(w, wp *Vector) float64            { return 0.0 }
func PhaseMieMurky(w, wp *Vector) float64           { return 0.0 }
func PhaseHG(w, wp *Vector, g float64) float64      { return 0.0 }
func PhaseSchlick(w, wp *Vector, g float64) float64 { return 0.0 }

func GetVolumeScatteringProperties(name string) (found bool, sigma_a, sigma_prime_s *Spectrum) {
	Unimplemented()
	return false, nil, nil
}
func SubsurfaceFromDiffuse(Kd *Spectrum, meanPathLength, eta float64) (sigma_a, sigma_prime_s *Spectrum) {
	Unimplemented()
	return nil, nil
}
