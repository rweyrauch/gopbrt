package core

type Sampler interface {
	GetMoreSamples(sample *Sample, rng *RNG) int
	MaximumSampleCount() int
	ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool
	GetSubSampler(num, count int) Sampler
	RoundSize(size int) int
}

type SamplerData struct {
	xPixelStart, xPixelEnd, yPixelStart, yPixelEnd int
	samplesPerPixel                                int
	shutterOpen, shutterClose                      float64
}

func (s *SamplerData) computeSubWindow(num, count int) (xstart, xend, ystart, yend int) {
	// Determine how many tiles to use in each dimension, _nx_ and _ny_
	dx := s.xPixelEnd - s.xPixelStart
	dy := s.yPixelEnd - s.yPixelStart
	nx, ny := count, 1
	for (nx&0x1) == 0 && 2*dx*ny < dy*nx {
		nx >>= 1
		ny <<= 1
	}

	// Compute $x$ and $y$ pixel sample range for sub-window
	xo, yo := num%nx, num/nx
	tx0, tx1 := float64(xo)/float64(nx), float64(xo+1)/float64(nx)
	ty0, ty1 := float64(yo)/float64(ny), float64(yo+1)/float64(ny)
	xstart = Floor2Int(Lerp(tx0, float64(s.xPixelStart), float64(s.xPixelEnd)))
	xend = Floor2Int(Lerp(tx1, float64(s.xPixelStart), float64(s.xPixelEnd)))
	ystart = Floor2Int(Lerp(ty0, float64(s.yPixelStart), float64(s.yPixelEnd)))
	yend = Floor2Int(Lerp(ty1, float64(s.yPixelStart), float64(s.yPixelEnd)))

	return xstart, xend, ystart, yend
}

type CameraSample struct {
	imageX, imageY float64
	lensU, lensV   float64
	time           float64
}

type Sample struct {
	CameraSample
	n1D, n2D []int
	oneD     [][]float64
	twoD     [][]float64
}

func CreateSample(sampler Sampler, surf SurfaceIntegrator, vol VolumeIntegrator, scene Scene) *Sample {
	// TODO: implement this
	return nil
}

func (s *Sample) Add1D(n int) int {
	s.n1D = append(s.n1D, n)
	return len(s.n1D)
}

func (s *Sample) Add2D(n int) int {
	s.n2D = append(s.n2D, n)
	return len(s.n2D)
}

func (s *Sample) Duplicate(count int) *Sample {
	// TODO: implement this
	return nil
}

const (
	ADAPTIVE_COMPARE_SHAPE_ID = iota
	ADAPTIVE_CONTRAST_THRESHOLD
)

type AdaptiveTest int

type (
	AdaptiveSampler struct {
		SamplerData
		xPos, yPos             int
		minSamples, maxSamples int
		sampleBuf              []float64

		method           AdaptiveTest
		supersamplePixel bool
	}

	BestCandidateSampler struct {
		SamplerData
		tableWidth                                 float64
		tableOffset                                int
		xTileStart, xTileEnd, yTileStart, yTileEnd int
		xTile, yTile                               int
		sampleOffsets                              [3]float64
	}

	HaltonSampler struct {
		SamplerData
		wantedSamples, currentSample int
	}

	LDSampler struct {
		SamplerData
		xPos, yPos, nPixelSamples int
		sampleBuf                 []float64
	}

	RandomSampler struct {
		SamplerData
		xPos, yPos, nPixelSamples              int
		imageSamples, lensSamples, timeSamples []float64
		samplePos                              int
	}

	StratifiedSampler struct {
		SamplerData
		xPixelSamples, yPixelSamples int
		jitterSamples                bool
		xPos, yPos                   int
		sampleBuf                    []float64
	}
)

func (s *AdaptiveSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *AdaptiveSampler) MaximumSampleCount() int                     { return 0 }
func (s *AdaptiveSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *AdaptiveSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *AdaptiveSampler) RoundSize(size int) int               { return 0 }

func (s *BestCandidateSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *BestCandidateSampler) MaximumSampleCount() int                     { return 0 }
func (s *BestCandidateSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *BestCandidateSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *BestCandidateSampler) RoundSize(size int) int               { return 0 }

func (s *HaltonSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *HaltonSampler) MaximumSampleCount() int                     { return 0 }
func (s *HaltonSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *HaltonSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *HaltonSampler) RoundSize(size int) int               { return 0 }

func (s *LDSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *LDSampler) MaximumSampleCount() int                     { return 0 }
func (s *LDSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *LDSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *LDSampler) RoundSize(size int) int               { return 0 }

func (s *RandomSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *RandomSampler) MaximumSampleCount() int                     { return 0 }
func (s *RandomSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *RandomSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *RandomSampler) RoundSize(size int) int               { return 0 }

func (s *StratifiedSampler) GetMoreSamples(sample *Sample, rng *RNG) int { return 0 }
func (s *StratifiedSampler) MaximumSampleCount() int                     { return 0 }
func (s *StratifiedSampler) ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool {
	return false
}
func (s *StratifiedSampler) GetSubSampler(num, count int) Sampler { return nil }
func (s *StratifiedSampler) RoundSize(size int) int               { return 0 }

func CreateAdaptiveSampler(params *ParamSet, film Film, camera Camera) *AdaptiveSampler { return nil }
func CreateBestCandidateSampler(params *ParamSet, film Film, camera Camera) *BestCandidateSampler {
	return nil
}
func CreateHaltonSampler(params *ParamSet, film Film, camera Camera) *HaltonSampler     { return nil }
func CreateLowDiscrepancySampler(params *ParamSet, film Film, camera Camera) *LDSampler { return nil }
func CreateRandomSampler(params *ParamSet, film Film, camera Camera) *RandomSampler     { return nil }
func CreateStratifiedSampler(params *ParamSet, film Film, camera Camera) *StratifiedSampler {
	return nil
}
