package pbrt

type Sampler interface {
	GetMoreSamples(sample *Sample, rng *RNG) int
	MaximumSampleCount() int
    ReportResults(samples []*Sample, rays []*RayDifferential, Ls *Spectrum, isects []*Intersection, count int) bool
    GetSubSampler(num, count int) *Sampler
    RoundSize(size int) int	
}

type SamplerData struct {
    xPixelStart, xPixelEnd, yPixelStart, yPixelEnd int
    samplesPerPixel int
    shutterOpen, shutterClose float64
}

func (s *SamplerData) computeSubWindow(num, count int) (xstart, xend, ystart, yend int) {
     // Determine how many tiles to use in each dimension, _nx_ and _ny_
    dx := s.xPixelEnd - s.xPixelStart
    dy := s.yPixelEnd - s.yPixelStart
    nx, ny := count, 1
    for ((nx & 0x1) == 0 && 2 * dx * ny < dy * nx) {
        nx >>= 1
        ny <<= 1
    }

    // Compute $x$ and $y$ pixel sample range for sub-window
    xo, yo := num % nx, num / nx
    tx0, tx1 := float64(xo) / float64(nx), float64(xo+1) / float64(nx)
    ty0, ty1 := float64(yo) / float64(ny), float64(yo+1) / float64(ny)
    xstart = Floor2Int(Lerp(tx0, float64(s.xPixelStart), float64(s.xPixelEnd)))
    xend   = Floor2Int(Lerp(tx1, float64(s.xPixelStart), float64(s.xPixelEnd)))
    ystart = Floor2Int(Lerp(ty0, float64(s.yPixelStart), float64(s.yPixelEnd)))
    yend   = Floor2Int(Lerp(ty1, float64(s.yPixelStart), float64(s.yPixelEnd)))
    
    return xstart, xend, ystart, yend
}

type CameraSample struct {
	imageX, imageY float64
	lensU, lensV float64
	time float64
}

type Sample struct {
	CameraSample
	n1D, n2D []int
    oneD [][]float64
    twoD [][]float64
}

func CreateSample(sampler *Sampler, surf *SurfaceIntegrator, vol *VolumeIntegrator, scene *Scene) *Sample {
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