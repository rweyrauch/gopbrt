package core

import (
	"math"
)

const (
	OneMinusEpsilon = 0.9999999403953552 // 0x1.fffffep-1
)

// Monte Carlo Utility Declarations
type Distribution1D struct {
	dfunc, cdf []float64
	dfuncInt   float64
	count      int
}

type Distribution2D struct {
	pConditionalV []*Distribution1D
	pMarginal     *Distribution1D
}

type PermutedHalton struct {
	dims       int
	b, permute []uint32
}

func (ph *PermutedHalton) Sample(n int, out []float64) {
	/*
	   	p := 0
	       for i := 0; i < ph.dims; i++ {
	           out[i] = math.Min(float64(PermutedRadicalInverse(n, ph.b[i], ph.permute[p])), OneMinusEpsilon)
	           p += ph.b[i]
	       }
	*/
}

func NewDistribution1D(f []float64) *Distribution1D {
	nd := new(Distribution1D)
	nd.count = len(f)
	nd.dfunc = make([]float64, nd.count, nd.count)
	_ = copy(nd.dfunc, f[0:])
	nd.cdf = make([]float64, nd.count+1, nd.count+1)
	// Compute integral of step function at $x_i$
	nd.cdf[0] = 0.0
	for i, df := range nd.dfunc {
		nd.cdf[i+1] = nd.cdf[i] + df/float64(nd.count)
	}

	// Transform step function integral into CDF
	nd.dfuncInt = nd.cdf[nd.count]
	if nd.dfuncInt == 0.0 {
		for i := 1; i < nd.count+1; i++ {
			nd.cdf[i] = float64(i) / float64(nd.count)
		}
	} else {
		for i := 1; i < nd.count+1; i++ {
			nd.cdf[i] /= nd.dfuncInt
		}
	}

	return nd
}

// See std::upper_bound (return index rather than iterator)
func upper_bound(cdf []float64, u float64) int {
	for i, v := range cdf {
		if v > u {
			return i
		}
	}
	return len(cdf) - 1
}

func (d *Distribution1D) SampleContinuous(u float64) (sample float64, pdf float64, offset int) {
	// Find surrounding CDF segments and _offset_
	offset = upper_bound(d.cdf, u)

	// Compute offset along CDF segment
	du := (u - d.cdf[offset]) / (d.cdf[offset+1] - d.cdf[offset])

	// Compute PDF for sampled offset
	pdf = d.dfunc[offset] / d.dfuncInt

	// Return $x\in{}[0,1)$ corresponding to sample
	sample = (float64(offset) + du) / float64(len(d.dfunc))

	return sample, pdf, offset
}

func (d *Distribution1D) SampleDiscrete(u float64) (offset int, pdf float64) {
	// Find surrounding CDF segments and _offset_
	offset = upper_bound(d.cdf, u)

	pdf = d.dfunc[offset] / (d.dfuncInt * float64(len(d.dfunc)))

	return offset, pdf
}

var primes = []int{
	// First 1000 prime numbers
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
	233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
	283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
	353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
	419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
	467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
	547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
	607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
	661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
	739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
	811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
	877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
	947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
	1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
	1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
	1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
	1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
	1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
	1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
	1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
	1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
	1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
	1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
	1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
	1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
	1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
	1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
	2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
	2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
	2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
	2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
	2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
	2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
	2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
	2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
	2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
	2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
	2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
	2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
	3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
	3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
	3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
	3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
	3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
	3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
	3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
	3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
	3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
	3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
	3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
	3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
	4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
	4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
	4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
	4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
	4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
	4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
	4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
	4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
	4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
	4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
	4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
	4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
	5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
	5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
	5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
	5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
	5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
	5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
	5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
	5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
	5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
	5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
	5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
	5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
	6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
	6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
	6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
	6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
	6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
	6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
	6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
	6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
	6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
	6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
	6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
	7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
	7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
	7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
	7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
	7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
	7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
	7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
	7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
	7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
	7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
}

func StratifiedSample1D(samp []float64, n int, rng *RNG, jitter bool) {
	invTot := 1.0 / float64(n)
	if jitter {
		for i := 0; i < n; i++ {
			delta := rng.RandomFloat()
			samp[i] = math.Min((float64(i)+delta)*invTot, OneMinusEpsilon)
		}
	} else {
		for i := 0; i < n; i++ {
			samp[i] = math.Min((float64(i)+0.5)*invTot, OneMinusEpsilon)
		}
	}
}

func StratifiedSample2D(samp []float64, nx, ny int, rng *RNG, jitter bool) {
	dx, dy := 1.0/float64(nx), 1.0/float64(ny)
	i := 0
	if jitter {
		for y := 0; y < ny; y++ {
			for x := 0; x < nx; x++ {
				jx := rng.RandomFloat()
				jy := rng.RandomFloat()
				samp[i] = math.Min((float64(x)+jx)*dx, OneMinusEpsilon)
				samp[i+1] = math.Min((float64(y)+jy)*dy, OneMinusEpsilon)
				i += 2
			}
		}
	} else {
		for y := 0; y < ny; y++ {
			for x := 0; x < nx; x++ {
				samp[i] = math.Min((float64(x)+0.5)*dx, OneMinusEpsilon)
				samp[i+1] = math.Min((float64(y)+0.5)*dy, OneMinusEpsilon)
				i += 2
			}
		}
	}
}

func LatinHypercube(samples []float64, nSamples, nDim uint32, rng *RNG) {
	// Generate LHS samples along diagonal
	delta := 1.0 / float64(nSamples)
	var i, j uint32
	for i = 0; i < nSamples; i++ {
		for j = 0; j < nDim; j++ {
			samples[nDim*i+j] = math.Min((float64(i)+(rng.RandomFloat()))*delta, OneMinusEpsilon)
		}
	}
	// Permute LHS samples in each dimension
	for i = 0; i < nDim; i++ {
		for j = 0; j < nSamples; j++ {
			other := j + (rng.RandomUInt() % (nSamples - j))
			samples[nDim*j+i], samples[nDim*other+i] = samples[nDim*other+i], samples[nDim*j+i]
		}
	}
}

func RadicalInverse(n, base int) float64 {
	val := 0.0
	invBase := 1.0 / base
	invBi := invBase

	for n > 0 {
		// Compute next digit of radical inverse
		d_i := (n % base)
		val += float64(d_i * invBi)
		n *= invBase
		invBi *= invBase
	}
	return val
}

func Shufflei(samp []uint32, count, dims uint32, rng *RNG) {
	var i, j uint32
	for i = 0; i < count; i++ {
		other := i + (rng.RandomUInt() % (count - i))
		for j = 0; j < dims; j++ {
			samp[dims*i+j], samp[dims*other+j] = samp[dims*other+j], samp[dims*i+j]
		}
	}
}

func Shuffle(samp []float64, count, dims uint32, rng *RNG) {
	var i, j uint32
	for i = 0; i < count; i++ {
		other := i + (rng.RandomUInt() % (count - i))
		for j = 0; j < dims; j++ {
			samp[dims*i+j], samp[dims*other+j] = samp[dims*other+j], samp[dims*i+j]
		}
	}
}

func GeneratePermutation(buf []uint32, b uint32, rng *RNG) {
	var i uint32
	for i = 0; i < b; i++ {
		buf[i] = i
	}
	Shufflei(buf, b, 1, rng)
}

func PermutedRadicalInverse(n uint32, base uint32, p []uint32) float64 {
	val := 0.0
	invBase := 1.0 / base
	invBi := invBase

	for n > 0 {
		d_i := p[n%base]
		val += float64(d_i * invBi)
		n *= invBase
		invBi *= invBase
	}
	return val
}

func LDPixelSampleFloatsNeeded(sample *Sample, nPixelSamples int) int {
	n := 5 // 2 lens + 2 pixel + time
	for i := 0; i < len(sample.n1D); i++ {
		n += sample.n1D[i]
	}
	for i := 0; i < len(sample.n2D); i++ {
		n += 2 * sample.n2D[i]
	}
	return nPixelSamples * n
}

func LDPixelSample(xPos, yPos int, shutterOpen, shutterClose float64, nPixelSamples int, samples []Sample, buf []float64, rng *RNG) {
	
	Assert(len(buf) >= nPixelSamples*5)
	
	// Prepare temporary array pointers for low-discrepancy camera samples
	imageSamples := buf[:2*nPixelSamples]
	lensSamples := buf[2*nPixelSamples:4*nPixelSamples]
	timeSamples := buf[4*nPixelSamples:5*nPixelSamples]
	Assert(len(imageSamples) == 2*nPixelSamples)
	Assert(len(lensSamples) == 2*nPixelSamples)
	Assert(len(timeSamples) == nPixelSamples)
	
	// Prepare temporary array pointers for low-discrepancy integrator samples
	//count1D := len(samples[0].n1D)
	//count2D := len(samples[0].n2D)
	
	/*
	   const uint32_t *n1D = count1D > 0 ? &samples[0].n1D[0] : NULL;
	   const uint32_t *n2D = count2D > 0 ? &samples[0].n2D[0] : NULL;
	   float **oneDSamples = ALLOCA(float *, count1D);
	   float **twoDSamples = ALLOCA(float *, count2D);
	   for (uint32_t i = 0; i < count1D; ++i) {
	       oneDSamples[i] = buf;
	       buf += n1D[i] * nPixelSamples;
	   }
	   for (uint32_t i = 0; i < count2D; ++i) {
	       twoDSamples[i] = buf;
	       buf += 2 * n2D[i] * nPixelSamples;
	   }
	*/
   // Generate low-discrepancy pixel samples
   LDShuffleScrambled2D(1, nPixelSamples, imageSamples, rng)
   LDShuffleScrambled2D(1, nPixelSamples, lensSamples, rng)
   LDShuffleScrambled1D(1, nPixelSamples, timeSamples, rng)
	 
	/*   
	   for (uint32_t i = 0; i < count1D; ++i)
	       LDShuffleScrambled1D(n1D[i], nPixelSamples, oneDSamples[i], rng);
	   for (uint32_t i = 0; i < count2D; ++i)
	       LDShuffleScrambled2D(n2D[i], nPixelSamples, twoDSamples[i], rng);
	*/
	   // Initialize _samples_ with computed sample values
	   for i := 0; i < nPixelSamples; i++ {
	       samples[i].imageX = float64(xPos) + imageSamples[2*i]
	       samples[i].imageY = float64(yPos) + imageSamples[2*i+1]
	       samples[i].time = Lerp(timeSamples[i], shutterOpen, shutterClose)
	       samples[i].lensU = lensSamples[2*i]
	       samples[i].lensV = lensSamples[2*i+1]
	       /*
	       // Copy integrator samples into _samples[i]_
	       for j := 0; j < count1D; j++ {
	           startSamp := n1D[j] * i;
	           for k := 0; k < n1D[j]; k++ {
	               samples[i].oneD[j][k] = oneDSamples[j][startSamp+k]
	            }   
	       }
	       for j := 0; j < count2D; j++ {
	           startSamp := 2 * n2D[j] * i;
	           for k := 0; k < 2*n2D[j]; k++ {
	               samples[i].twoD[j][k] = twoDSamples[j][startSamp+k]
	           }
	       }
	       */
	   }	
}

// Sampling Inline Functions
func Sample02(n uint32, scramble [2]uint32, sample []float64) {
	sample[0] = VanDerCorput(n, scramble[0])
	sample[1] = Sobol2(n, scramble[1])
}

func VanDerCorput(n, scramble uint32) float64 {
	// Reverse bits of _n_
	n = (n << 16) | (n >> 16)
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8)
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4)
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2)
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1)
	n ^= scramble
	return math.Min(float64((n>>8)&0xffffff)/float64(1<<24), OneMinusEpsilon)
}

func Sobol2(n, scramble uint32) float64 {
	var v uint32
	for v = 1 << 31; n != 0; n >>= 1 {
		if n&0x1 != 0 {
			scramble ^= v
		}
		v ^= v >> 1
	}
	return math.Min((float64((scramble>>8)&0xffffff))/float64(1<<24), OneMinusEpsilon)
}

func LDShuffleScrambled1D(nSamples, nPixel int, samples []float64, rng *RNG) {
	scramble := rng.RandomUInt()
	for i := 0; i < nSamples*nPixel; i++ {
		samples[i] = VanDerCorput(uint32(i), scramble)
	}
	for i := 0; i < nPixel; i++ {
		Shuffle(samples[i*nSamples:], uint32(nSamples), 1, rng)
	}
	Shuffle(samples, uint32(nPixel), uint32(nSamples), rng)
}

func LDShuffleScrambled2D(nSamples, nPixel int, samples []float64, rng *RNG) {
	scramble := [2]uint32{rng.RandomUInt(), rng.RandomUInt()}
	for i := 0; i < nSamples*nPixel; i++ {
		Sample02(uint32(i), scramble, samples[2*i:])
	}

	for i := 0; i < nPixel; i++ {
		Shuffle(samples[2*i*nSamples:], uint32(nSamples), 2, rng)
	}
	Shuffle(samples, uint32(nPixel), uint32(2*nSamples), rng)
}

// Monte Carlo Function Definitions
func RejectionSampleDisk(rng *RNG) (x, y float64) {
	sx, sy := 1.0, 1.0
	for sx*sx+sy*sy > 1.0 {
		sx = 1.0 - 2.0*rng.RandomFloat()
		sy = 1.0 - 2.0*rng.RandomFloat()
	}
	x = sx
	y = sy
	return x, y
}

func UniformSampleHemisphere(u1, u2 float64) *Vector {
	z := u1
	r := math.Sqrt(math.Max(0.0, 1.0-z*z))
	phi := 2 * math.Pi * u2
	x := r * math.Cos(phi)
	y := r * math.Sin(phi)
	return &Vector{x, y, z}
}

func UniformHemispherePdf() float64 {
	return 1.0 / (2.0 * math.Pi)
}

func UniformSampleSphere(u1, u2 float64) *Vector {
	z := 1.0 - 2.0*u1
	r := math.Sqrt(math.Max(0.0, 1.0-z*z))
	phi := 2.0 * math.Pi * u2
	x := r * math.Cos(phi)
	y := r * math.Sin(phi)
	return &Vector{x, y, z}
}

func UniformSpherePdf() float64 {
	return 1.0 / (4.0 * math.Pi)
}

func UniformSampleDisk(u1, u2 float64) (x, y float64) {
	r := math.Sqrt(u1)
	theta := 2.0 * math.Pi * u2
	x = r * math.Cos(theta)
	y = r * math.Sin(theta)
	return x, y
}

func ConcentricSampleDisk(u1, u2 float64) (dx, dy float64) {
	var r, theta float64
	// Map uniform random numbers to $[-1,1]^2$
	sx := 2*u1 - 1
	sy := 2*u2 - 1

	// Map square to $(r,\theta)$

	// Handle degeneracy at the origin
	if sx == 0.0 && sy == 0.0 {
		return 0.0, 0.0
	}
	if sx >= -sy {
		if sx > sy {
			// Handle first region of disk
			r = sx
			if sy > 0.0 {
				theta = sy / r
			} else {
				theta = 8.0 + sy/r
			}
		} else {
			// Handle second region of disk
			r = sy
			theta = 2.0 - sx/r
		}
	} else {
		if sx <= sy {
			// Handle third region of disk
			r = -sx
			theta = 4.0 - sy/r
		} else {
			// Handle fourth region of disk
			r = -sy
			theta = 6.0 + sx/r
		}
	}
	theta *= math.Pi / 4.0
	dx = r * math.Cos(theta)
	dy = r * math.Sin(theta)
	return dx, dy
}

func CosineSampleHemisphere(u1, u2 float64) *Vector {
	x, y := ConcentricSampleDisk(u1, u2)
	z := math.Sqrt(math.Max(0.0, 1.0-x*x-y*y))
	return &Vector{x, y, z}
}

func CosineHemispherePdf(costheta, phi float64) float64 {
	return costheta * 1.0 / math.Pi
}

func UniformSampleTriangle(u1, u2 float64) (u, v float64) {
	su1 := math.Sqrt(u1)
	u = 1.0 - su1
	v = u2 * su1
	return u, v
}

func NewDistribution2D(dfunc []float64, nu, nv int) *Distribution2D {
	d := new(Distribution2D)
	d.pConditionalV = make([]*Distribution1D, 0, nv)

	for v := 0; v < nv; v++ {
		// Compute conditional sampling distribution for $\tilde{v}$
		d.pConditionalV = append(d.pConditionalV, NewDistribution1D(dfunc[v*nu:v*nu+nu]))
	}
	// Compute marginal sampling distribution $p[\tilde{v}]$
	marginalFunc := make([]float64, nv, nv)
	for v := 0; v < nv; v++ {
		marginalFunc[v] = d.pConditionalV[v].dfuncInt
	}
	d.pMarginal = NewDistribution1D(marginalFunc)

	return d
}

func (d *Distribution2D) SampleContinuous(u0, u1 float64) (uv [2]float64, pdf float64) {
	var pdfs [2]float64
	var v int
	uv[1], pdfs[1], v = d.pMarginal.SampleContinuous(u1)
	uv[0], pdfs[0], _ = d.pConditionalV[v].SampleContinuous(u0)
	pdf = pdfs[0] * pdfs[1]
	return uv, pdf
}

func (d *Distribution2D) Pdf(u, v float64) float64 {
	iu := Clampi(Float2Int(u*float64(d.pConditionalV[0].count)), 0, d.pConditionalV[0].count-1)
	iv := Clampi(Float2Int(v*float64(d.pMarginal.count)), 0, d.pMarginal.count-1)
	if d.pConditionalV[iv].dfuncInt*d.pMarginal.dfuncInt == 0.0 {
		return 0.0
	}
	return (d.pConditionalV[iv].dfunc[iu] * d.pMarginal.dfunc[iv]) /
		(d.pConditionalV[iv].dfuncInt * d.pMarginal.dfuncInt)
}

func UniformConePdf(cosThetaMax float64) float64 {
	return 1.0 / (2.0 * math.Pi * (1.0 - cosThetaMax))
}

func UniformSampleCone(u1, u2, costhetamax float64) *Vector {
	costheta := (1.0 - u1) + u1*costhetamax
	sintheta := math.Sqrt(1.0 - costheta*costheta)
	phi := u2 * 2.0 * math.Pi
	return &Vector{math.Cos(phi) * sintheta, math.Sin(phi) * sintheta, costheta}
}

func UniformSampleConeVector(u1, u2, costhetamax float64, x, y, z *Vector) *Vector {
	costheta := Lerp(u1, costhetamax, 1.0)
	sintheta := math.Sqrt(1.0 - costheta*costheta)
	phi := u2 * 2.0 * math.Pi
	return x.Scale(math.Cos(phi) * sintheta).Add(y.Scale(math.Sin(phi) * sintheta)).Add(z.Scale(costheta))
}

func SampleHG(w *Vector, g, u1, u2 float64) *Vector {
	var costheta float64
	if math.Abs(g) < 1.0e-3 {
		costheta = 1.0 - 2.0*u1
	} else {
		sqrTerm := (1.0 - g*g) / (1.0 - g + 2.0*g*u1)
		costheta = (1.0 + g*g - sqrTerm*sqrTerm) / (2.0 * g)
	}
	sintheta := math.Sqrt(math.Max(0.0, 1.0-costheta*costheta))
	phi := 2.0 * math.Pi * u2
	v1, v2 := CoordinateSystem(w)
	return SphericalDirectionVectors(sintheta, costheta, phi, v1, v2, w)
}

func HGPdf(w, wp *Vector, g float64) float64 {
	return PhaseHG(w, wp, g)
}

func PowerHeuristic(nf int, fPdf float64, ng int, gPdf float64) float64 {
    f, g := float64(nf) * fPdf, float64(ng) * gPdf
    return (f*f) / (f*f + g*g)
}
