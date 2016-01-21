package core

import (
	"math"
	"strings"
)

type Sampler interface {
	GetMoreSamples(samples []Sample, rng *RNG) int
	MaximumSampleCount() int
	ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool
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

func NewSample(sampler Sampler, surf SurfaceIntegrator, vol VolumeIntegrator, scene *Scene) *Sample {
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
		xPos, yPos, nSamples                   int
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

func NewAdaptiveSampler(xstart, xend, ystart, yend, minsamp, maxsamp int, method string, sopen, sclose float64) *AdaptiveSampler {
	sampler := new(AdaptiveSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = int(RoundUpPow2(uint32(Maxi(minsamp, maxsamp))))
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose
	sampler.xPos = sampler.xPixelStart
	sampler.yPos = sampler.yPixelStart
	sampler.supersamplePixel = false

	if minsamp > maxsamp {
		minsamp, maxsamp = maxsamp, minsamp
	}

	if !IsPowerOf2(minsamp) {
		Warning("Minimum pixel samples being rounded up to power of 2")
		sampler.minSamples = int(RoundUpPow2(uint32(minsamp)))
	} else {
		sampler.minSamples = minsamp
	}
	if !IsPowerOf2(maxsamp) {
		Warning("Maximum pixel samples being rounded up to power of 2")
		sampler.maxSamples = int(RoundUpPow2(uint32(maxsamp)))
	} else {
		sampler.maxSamples = maxsamp
	}

	if sampler.minSamples < 2 {
		Warning("Adaptive sampler needs at least two initial pixel samples.  Using two.")
		sampler.minSamples = 2
	}
	if sampler.minSamples == sampler.maxSamples {
		sampler.maxSamples *= 2
		Warning("Adaptive sampler must have more maximum samples than minimum.  Using %d - %d",
			sampler.minSamples, sampler.maxSamples)
	}
	if strings.Compare(method, "contrast") == 0 {
		sampler.method = ADAPTIVE_CONTRAST_THRESHOLD
	} else if strings.Compare(method, "shapeid") == 0 {
		sampler.method = ADAPTIVE_COMPARE_SHAPE_ID
	} else {
		Warning("Adaptive sampling metric \"%s\" unknown.  Using \"contrast\".", method)
		sampler.method = ADAPTIVE_CONTRAST_THRESHOLD
	}
	sampler.sampleBuf = nil

	return sampler
}

func (s *AdaptiveSampler) GetMoreSamples(samples []Sample, rng *RNG) int {
	if s.sampleBuf == nil {
		samplesNeeded := LDPixelSampleFloatsNeeded(&samples[0], s.maxSamples)
		s.sampleBuf = make([]float64, samplesNeeded, samplesNeeded)
	}

	if s.supersamplePixel {
		LDPixelSample(s.xPos, s.yPos, s.shutterOpen, s.shutterClose, s.maxSamples, samples, s.sampleBuf, rng)
		return s.maxSamples
	} else {
		if s.yPos == s.yPixelEnd {
			return 0
		}
		LDPixelSample(s.xPos, s.yPos, s.shutterOpen, s.shutterClose, s.minSamples, samples, s.sampleBuf, rng)
		return s.minSamples
	}
}

func (s *AdaptiveSampler) MaximumSampleCount() int { return s.maxSamples }

func (s *AdaptiveSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	if s.supersamplePixel {
		s.supersamplePixel = false
		// Advance to next pixel for sampling for _AdaptiveSampler_
		s.xPos++
		if s.xPos == s.xPixelEnd {
			s.xPos = s.xPixelStart
			s.yPos++
		}
		return true
	} else if s.needsSupersampling(samples, rays, Ls, isects, count) {
		//PBRT_SUPERSAMPLE_PIXEL_YES(xPos, yPos);
		s.supersamplePixel = true
		return false
	} else {
		//PBRT_SUPERSAMPLE_PIXEL_NO(xPos, yPos);
		// Advance to next pixel for sampling for _AdaptiveSampler_
		s.xPos++
		if s.xPos == s.xPixelEnd {
			s.xPos = s.xPixelStart
			s.yPos++
		}
		return true
	}
}

func (s *AdaptiveSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	method := "shapeid"
	if s.method == ADAPTIVE_CONTRAST_THRESHOLD {
		method = "constrast"
	}
	return NewAdaptiveSampler(x0, x1, y0, y1, s.minSamples, s.maxSamples, method, s.shutterOpen, s.shutterClose)
}

func (s *AdaptiveSampler) RoundSize(size int) int {
	return int(RoundUpPow2(uint32(size)))
}

func (s *AdaptiveSampler) needsSupersampling(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	switch s.method {
	case ADAPTIVE_COMPARE_SHAPE_ID:
		// See if any shape ids differ within samples
		for i := 0; i < count-1; i++ {
			if isects[i].shapeId != isects[i+1].shapeId ||
				isects[i].primitiveId != isects[i+1].primitiveId {
				return true
			}
		}
		return false
	case ADAPTIVE_CONTRAST_THRESHOLD:
		// Compare contrast of sample differences to threshold
		Lavg := 0.0
		for i := 0; i < count; i++ {
			Lavg += float64(Ls[i].Y())
		}
		Lavg /= float64(count)
		maxContrast := 0.5
		for i := 0; i < count; i++ {
			if math.Abs(float64(Ls[i].Y())-Lavg)/Lavg > maxContrast {
				return true
			}
		}
		return false
	}
	return false
}

func NewBestCandidateSampler(xstart, xend, ystart, yend, ps int, sopen, sclose float64) *BestCandidateSampler {
	sampler := new(BestCandidateSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = ps
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose

	sampler.tableWidth = SQRT_SAMPLE_TABLE_SIZE / math.Sqrt(float64(ps))
	sampler.xTileStart = Floor2Int(float64(xstart) / sampler.tableWidth)
	sampler.xTileEnd = Floor2Int(float64(xend) / sampler.tableWidth)
	sampler.yTileStart = Floor2Int(float64(ystart) / sampler.tableWidth)
	sampler.yTileEnd = Floor2Int(float64(yend) / sampler.tableWidth)
	sampler.xTile = sampler.xTileStart
	sampler.yTile = sampler.yTileStart
	sampler.tableOffset = 0
	// Update sample shifts
	tileRng := CreateRNG(int64(sampler.xTile + (sampler.yTile << 8)))
	for i := 0; i < 3; i++ {
		sampler.sampleOffsets[i] = tileRng.RandomFloat()
	}

	return sampler
}

func (s *BestCandidateSampler) GetMoreSamples(sample []Sample, rng *RNG) int {
again:
	if s.tableOffset == SAMPLE_TABLE_SIZE {
		// Advance to next best-candidate sample table position
		s.tableOffset = 0
		s.xTile++
		if s.xTile > s.xTileEnd {
			s.xTile = s.xTileStart
			s.yTile++
			if s.yTile > s.yTileEnd {
				return 0
			}
		}

		// Update sample shifts
		tileRng := CreateRNG(int64(s.xTile + (s.yTile << 8)))
		for i := 0; i < 3; i++ {
			s.sampleOffsets[i] = tileRng.RandomFloat()
		}
	}
	// Compute raster sample from table
	WRAP := func(x float64) float64 {
		if x >= 1.0 {
			return x - 1.0
		} else {
			return x
		}
	}
	sample[0].imageX = (float64(s.xTile) + bestCandidateSampleTable[s.tableOffset][0]) * s.tableWidth
	sample[0].imageY = (float64(s.yTile) + bestCandidateSampleTable[s.tableOffset][1]) * s.tableWidth
	sample[0].time = Lerp(WRAP(s.sampleOffsets[0]+bestCandidateSampleTable[s.tableOffset][2]), s.shutterOpen, s.shutterClose)
	sample[0].lensU = WRAP(s.sampleOffsets[1] + bestCandidateSampleTable[s.tableOffset][3])
	sample[0].lensV = WRAP(s.sampleOffsets[2] + bestCandidateSampleTable[s.tableOffset][4])

	// Check sample against crop window, goto _again_ if outside
	if sample[0].imageX < float64(s.xPixelStart) || sample[0].imageX >= float64(s.xPixelEnd) ||
		sample[0].imageY < float64(s.yPixelStart) || sample[0].imageY >= float64(s.yPixelEnd) {
		s.tableOffset++
		goto again
	}

	// Compute integrator samples for best-candidate sample
	for i := 0; i < len(sample[0].n1D); i++ {
		LDShuffleScrambled1D(sample[0].n1D[i], 1, sample[0].oneD[i], rng)
	}
	for i := 0; i < len(sample[0].n2D); i++ {
		LDShuffleScrambled2D(sample[0].n2D[i], 1, sample[0].twoD[i], rng)
	}
	s.tableOffset++
	return 1
}

func (s *BestCandidateSampler) MaximumSampleCount() int { return 1 }

func (s *BestCandidateSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	return false
}

func (s *BestCandidateSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	return NewBestCandidateSampler(x0, x1, y0, y1, s.samplesPerPixel, s.shutterOpen, s.shutterClose)
}

func (s *BestCandidateSampler) RoundSize(size int) int {
	return int(RoundUpPow2(uint32(size)))
}

func NewHaltonSampler(xstart, xend, ystart, yend, ps int, sopen, sclose float64) *HaltonSampler {
	sampler := new(HaltonSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = ps
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose

	delta := Maxi(sampler.xPixelEnd-sampler.xPixelStart, sampler.yPixelEnd-sampler.yPixelStart)
	sampler.wantedSamples = sampler.samplesPerPixel * delta * delta
	sampler.currentSample = 0

	return sampler
}

func (s *HaltonSampler) GetMoreSamples(sample []Sample, rng *RNG) int {
retry:
	if s.currentSample >= s.wantedSamples {
		return 0
	}
	// Generate sample with Halton sequence and reject if outside image extent
	u := float64(RadicalInverse(s.currentSample, 3))
	v := float64(RadicalInverse(s.currentSample, 2))
	lerpDelta := float64(Maxi(s.xPixelEnd-s.xPixelStart, s.yPixelEnd-s.yPixelStart))
	sample[0].imageX = Lerp(u, float64(s.xPixelStart), float64(s.xPixelStart)+lerpDelta)
	sample[0].imageY = Lerp(v, float64(s.yPixelStart), float64(s.yPixelStart)+lerpDelta)
	s.currentSample++
	if sample[0].imageX >= float64(s.xPixelEnd) || sample[0].imageY >= float64(s.yPixelEnd) {
		goto retry
	}

	// Generate lens, time, and integrator samples for _HaltonSampler_
	sample[0].lensU = float64(RadicalInverse(s.currentSample, 5))
	sample[0].lensV = float64(RadicalInverse(s.currentSample, 7))
	sample[0].time = Lerp(float64(RadicalInverse(s.currentSample, 11)), s.shutterOpen, s.shutterClose)
	for i := 0; i < len(sample[0].n1D); i++ {
		LatinHypercube(sample[0].oneD[i], uint32(sample[0].n1D[i]), 1, rng)
	}
	for i := 0; i < len(sample[0].n2D); i++ {
		LatinHypercube(sample[0].twoD[i], uint32(sample[0].n2D[i]), 2, rng)
	}
	return 1
}

func (s *HaltonSampler) MaximumSampleCount() int { return 1 }
func (s *HaltonSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	return false
}
func (s *HaltonSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	return NewHaltonSampler(x0, x1, y0, y1, s.samplesPerPixel, s.shutterOpen, s.shutterClose)
}
func (s *HaltonSampler) RoundSize(size int) int { return size }

func NewLDSampler(xstart, xend, ystart, yend, ps int, sopen, sclose float64) *LDSampler {
	sampler := new(LDSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = int(RoundUpPow2(uint32(ps)))
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose
	if !IsPowerOf2(ps) {
		Warning("Pixel samples being rounded up to power of 2.")
		sampler.nPixelSamples = int(RoundUpPow2(uint32(ps)))
	} else {
		sampler.nPixelSamples = ps
	}
	sampler.xPos = sampler.xPixelStart
	sampler.yPos = sampler.yPixelStart
	sampler.sampleBuf = nil
	return sampler
}

func (s *LDSampler) GetMoreSamples(samples []Sample, rng *RNG) int {
	if s.yPos == s.yPixelEnd {
		return 0
	}
	if s.sampleBuf == nil {
		samplesNeeded := LDPixelSampleFloatsNeeded(&samples[0], s.nPixelSamples)
		s.sampleBuf = make([]float64, samplesNeeded, samplesNeeded)
	}
	LDPixelSample(s.xPos, s.yPos, s.shutterOpen, s.shutterClose, s.nPixelSamples, samples, s.sampleBuf, rng)
	s.xPos++
	if s.xPos == s.xPixelEnd {
		s.xPos = s.xPixelStart
		s.yPos++
	}
	return s.nPixelSamples
}
func (s *LDSampler) MaximumSampleCount() int { return s.nPixelSamples }
func (s *LDSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	return false
}
func (s *LDSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	return NewLDSampler(x0, x1, y0, y1, s.nPixelSamples, s.shutterOpen, s.shutterClose)

}
func (s *LDSampler) RoundSize(size int) int { return int(RoundUpPow2(uint32(size))) }

func NewRandomSampler(xstart, xend, ystart, yend, ns int, sopen, sclose float64) *RandomSampler {
	sampler := new(RandomSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = ns
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose
	sampler.xPos = sampler.xPixelStart
	sampler.yPos = sampler.yPixelStart
	sampler.nSamples = ns

	// Get storage for a pixel's worth of stratified samples
	sampler.imageSamples = make([]float64, 2*sampler.nSamples, 2*sampler.nSamples) // (x,y)
	sampler.lensSamples = make([]float64, 2*sampler.nSamples, 2*sampler.nSamples)  // (x,y)
	sampler.timeSamples = make([]float64, sampler.nSamples, sampler.nSamples)      // t

	rng := CreateRNG(int64(xstart + ystart*(xend-xstart)))
	for i, _ := range sampler.imageSamples {
		sampler.imageSamples[i] = rng.RandomFloat()
	}
	for i, _ := range sampler.lensSamples {
		sampler.lensSamples[i] = rng.RandomFloat()
	}
	for i, _ := range sampler.timeSamples {
		sampler.timeSamples[i] = rng.RandomFloat()
	}

	// Shift image samples to pixel coordinates
	for o := 0; o < 2*sampler.nSamples; o = o + 2 {
		sampler.imageSamples[o] += float64(sampler.xPos)
		sampler.imageSamples[o+1] += float64(sampler.yPos)
	}
	sampler.samplePos = 0
	return sampler
}

func (s *RandomSampler) GetMoreSamples(sample []Sample, rng *RNG) int {
	if s.samplePos == s.nSamples {
		if s.xPixelStart == s.xPixelEnd || s.yPixelStart == s.yPixelEnd {
			return 0
		}
		s.xPos++
		if s.xPos == s.xPixelEnd {
			s.xPos = s.xPixelStart
			s.yPos++
		}
		if s.yPos == s.yPixelEnd {
			return 0
		}

		for i, _ := range s.imageSamples {
			s.imageSamples[i] = rng.RandomFloat()
		}
		for i, _ := range s.lensSamples {
			s.lensSamples[i] = rng.RandomFloat()
		}
		for i, _ := range s.timeSamples {
			s.timeSamples[i] = rng.RandomFloat()
		}

		// Shift image samples to pixel coordinates
		for o := 0; o < 2*s.nSamples; o = o + 2 {
			s.imageSamples[o] += float64(s.xPos)
			s.imageSamples[o+1] += float64(s.yPos)
		}
		s.samplePos = 0
	}
	// Return next \mono{RandomSampler} sample point
	sample[0].imageX = s.imageSamples[2*s.samplePos]
	sample[0].imageY = s.imageSamples[2*s.samplePos+1]
	sample[0].lensU = s.lensSamples[2*s.samplePos]
	sample[0].lensV = s.lensSamples[2*s.samplePos+1]
	sample[0].time = Lerp(s.timeSamples[s.samplePos], s.shutterOpen, s.shutterClose)
	// Generate stratified samples for integrators
	for i := 0; i < len(sample[0].n1D); i++ {
		for j := 0; j < sample[0].n1D[i]; j++ {
			sample[0].oneD[i][j] = rng.RandomFloat()
		}
	}
	for i := 0; i < len(sample[0].n2D); i++ {
		for j := 0; j < 2*sample[0].n2D[i]; j++ {
			sample[0].twoD[i][j] = rng.RandomFloat()
		}
	}
	s.samplePos++
	return 1
}

func (s *RandomSampler) MaximumSampleCount() int { return 1 }
func (s *RandomSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	return false
}
func (s *RandomSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	return NewRandomSampler(x0, x1, y0, y1, s.nSamples, s.shutterOpen, s.shutterClose)
}

func (s *RandomSampler) RoundSize(size int) int { return size }

func NewStratifiedSampler(xstart, xend, ystart, yend, xsamp, ysamp int, jitter bool, sopen, sclose float64) *StratifiedSampler {
	sampler := new(StratifiedSampler)
	sampler.xPixelStart = xstart
	sampler.xPixelEnd = xend
	sampler.yPixelStart = ystart
	sampler.yPixelEnd = yend
	sampler.samplesPerPixel = xsamp * ysamp
	sampler.shutterOpen = sopen
	sampler.shutterClose = sclose
	sampler.xPos = sampler.xPixelStart
	sampler.yPos = sampler.yPixelStart
	sampler.xPixelSamples = xsamp
	sampler.yPixelSamples = ysamp
	sampler.sampleBuf = make([]float64, 5*xsamp*ysamp, 5*xsamp*ysamp)

	return sampler
}
func (s *StratifiedSampler) GetMoreSamples(samples []Sample, rng *RNG) int {
	if s.yPos == s.yPixelEnd {
		return 0
	}
	nSamples := s.xPixelSamples * s.yPixelSamples
	// Generate stratified camera samples for _(xPos, yPos)_

	// Generate initial stratified samples into _sampleBuf_ memory
	imageSamples := s.sampleBuf[:2*nSamples]
	lensSamples := s.sampleBuf[2*nSamples : 2*nSamples]
	timeSamples := s.sampleBuf[4*nSamples:]
	StratifiedSample2D(imageSamples, s.xPixelSamples, s.yPixelSamples, rng, s.jitterSamples)
	StratifiedSample2D(lensSamples, s.xPixelSamples, s.yPixelSamples, rng, s.jitterSamples)
	StratifiedSample1D(timeSamples, s.xPixelSamples*s.yPixelSamples, rng, s.jitterSamples)

	// Shift stratified image samples to pixel coordinates
	for o := 0; o < 2*s.xPixelSamples*s.yPixelSamples; o = o + 2 {
		imageSamples[o] += float64(s.xPos)
		imageSamples[o+1] += float64(s.yPos)
	}

	// Decorrelate sample dimensions
	Shuffle(lensSamples, uint32(nSamples), 2, rng)
	Shuffle(timeSamples, uint32(nSamples), 1, rng)

	// Initialize stratified _samples_ with sample values
	for i := 0; i < nSamples; i++ {
		samples[i].imageX = imageSamples[2*i]
		samples[i].imageY = imageSamples[2*i+1]
		samples[i].lensU = lensSamples[2*i]
		samples[i].lensV = lensSamples[2*i+1]
		samples[i].time = Lerp(timeSamples[i], s.shutterOpen, s.shutterClose)
		// Generate stratified samples for integrators
		for j := 0; j < len(samples[i].n1D); j++ {
			LatinHypercube(samples[i].oneD[j], uint32(samples[i].n1D[j]), 1, rng)
		}
		for j := 0; j < len(samples[i].n2D); j++ {
			LatinHypercube(samples[i].twoD[j], uint32(samples[i].n2D[j]), 2, rng)
		}
	}

	// Advance to next pixel for stratified sampling
	s.xPos++
	if s.xPos == s.xPixelEnd {
		s.xPos = s.xPixelStart
		s.yPos++
	}
	return nSamples
}

func (s *StratifiedSampler) MaximumSampleCount() int { return s.xPixelSamples * s.yPixelSamples }

func (s *StratifiedSampler) ReportResults(samples []Sample, rays []RayDifferential, Ls []Spectrum, isects []Intersection, count int) bool {
	return false
}

func (s *StratifiedSampler) GetSubSampler(num, count int) Sampler {
	x0, x1, y0, y1 := s.computeSubWindow(num, count)
	if x0 == x1 || y0 == y1 {
		return nil
	}
	return NewStratifiedSampler(x0, x1, y0, y1, s.xPixelSamples, s.yPixelSamples, s.jitterSamples, s.shutterOpen, s.shutterClose)
}

func (s *StratifiedSampler) RoundSize(size int) int { return size }

func CreateAdaptiveSampler(params *ParamSet, film Film, camera Camera) *AdaptiveSampler {
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	minsamp := params.FindIntParam("minsamples", 4)
	maxsamp := params.FindIntParam("maxsamples", 32)
	if options.QuickRender {
		minsamp = 2
		maxsamp = 4
	}
	method := params.FindStringParam("method", "contrast")
	return NewAdaptiveSampler(xstart, xend, ystart, yend, minsamp, maxsamp, method, camera.ShutterOpen(), camera.ShutterClose())
}

func CreateBestCandidateSampler(params *ParamSet, film Film, camera Camera) *BestCandidateSampler {
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	nsamp := params.FindIntParam("pixelsamples", 4)
	if options.QuickRender {
		nsamp = 1
	}
	return NewBestCandidateSampler(xstart, xend, ystart, yend, nsamp, camera.ShutterOpen(), camera.ShutterClose())
}

func CreateHaltonSampler(params *ParamSet, film Film, camera Camera) *HaltonSampler {
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	nsamp := params.FindIntParam("pixelsamples", 4)
	if options.QuickRender {
		nsamp = 1
	}
	return NewHaltonSampler(xstart, xend, ystart, yend, nsamp, camera.ShutterOpen(), camera.ShutterClose())
}

func CreateLowDiscrepancySampler(params *ParamSet, film Film, camera Camera) *LDSampler {
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	nsamp := params.FindIntParam("pixelsamples", 4)
	if options.QuickRender {
		nsamp = 1
	}
	return NewLDSampler(xstart, xend, ystart, yend, nsamp, camera.ShutterOpen(), camera.ShutterClose())
}

func CreateRandomSampler(params *ParamSet, film Film, camera Camera) *RandomSampler {
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	nsamp := params.FindIntParam("pixelsamples", 4)
	if options.QuickRender {
		nsamp = 1
	}
	return NewRandomSampler(xstart, xend, ystart, yend, nsamp, camera.ShutterOpen(), camera.ShutterClose())
}

func CreateStratifiedSampler(params *ParamSet, film Film, camera Camera) *StratifiedSampler {
	jitter := params.FindBoolParam("jitter", true)
	// Initialize common sampler parameters
	xstart, xend, ystart, yend := film.GetSampleExtent()
	xsamp := params.FindIntParam("xsamples", 2)
	ysamp := params.FindIntParam("ysamples", 2)
	if options.QuickRender {
		xsamp, ysamp = 1, 1
	}
	return NewStratifiedSampler(xstart, xend, ystart, yend, xsamp, ysamp, jitter, camera.ShutterOpen(), camera.ShutterClose())
}
