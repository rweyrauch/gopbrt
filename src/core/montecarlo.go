package pbrt

const (
	OneMinusEpsilon = 0.9999999403953552 // 0x1.fffffep-1
)

// Monte Carlo Utility Declarations
type Distribution1D struct {
	dfunc, cdf []float64
	dfuncInt float64
}


type Distribution2D struct {
	
}

func CreateDistribution1D(f []float64) {
	nd := new(Distribution1D)
	n := len(f)
	nd.dfunc = make([]float64, n, n)
	num := copy(nd.dfunc, f[0:])
	nd.cdf = make([]float64, n+1, n+1)
	// Compute integral of step function at $x_i$
	nd.cdf[0] = 0.0
	for df, i := range nd.dfunc {
		nd.cdf[i+1] = nd.cdf[i] + df / float64(n)
	}

    // Transform step function integral into CDF
    nd.dfuncInt = nd.cdf[n]
    if nd.dfuncInt == 0.0 {
        for i := 1; i < n+1; i++ {
            nd.cdf[i] = float64(i) / float64(n)
        }
    } else {
        for i := 1; i < n+1; i++ {
            nd.cdf[i] /= nd.dfuncInt
		}            
    }
	
	return nd
}

// See std::upper_bound (return index rather than iterator)
func upper_bound(cdf []float64, u float64) int {
	for v, i := range cdf {
		if v > u {
			return i
		}
	}
	return len(cdf)-1
}

func (d *Distribution1D) SampleContinuous(u float64) (sample float64, pdf float64, offset int) {
    // Find surrounding CDF segments and _offset_
    offset = upper_bound(d.cdf, u)

    // Compute offset along CDF segment
    du := (u - d.cdf[offset]) / (d.cdf[offset+1] - d.cdf[offset])
 
    // Compute PDF for sampled offset
    pdf = d.dfunc[offset] / d.dfuncInt

    // Return $x\in{}[0,1)$ corresponding to sample
    sample = (offset + du) / float64(len(d.dfunc))
    
    return sample, pdf, offset
}

func (d *Distribution1D) SampleDiscrete(u float64) (offset int, pdf float64) {
    // Find surrounding CDF segments and _offset_
    offset = upper_bound(d.cdf, u)

    pdf = d.dfunc[offset] / (d.dfuncInt * float64(len(d.dfunc)))
    
    return offset, pdf
}

