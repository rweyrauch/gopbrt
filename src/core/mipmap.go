package core

import (
	"math"
)

const (
	WEIGHT_LUT_SIZE = 128
)
const (
	TEXTURE_BLACK WrapMode = iota
	TEXTURE_CLAMP
	TEXTURE_REPEAT
)

type WrapMode int

type ResampleWeight struct {
	firstTexel int
	weight     [4]float32
}

type MIPMapData struct {
	doTrilinear            bool
	maxAnisotropy          float64
	wrapMode               WrapMode
	width, height, nLevels int
}
type MIPMapFloat struct {
	MIPMapData
	pyramid      [][]float32
	pyramidSizes [][2]int
}

type MIPMapSpectrum struct {
	MIPMapData
	pyramid      [][]Spectrum
	pyramidSizes [][2]int
}

func init() {
	// Initialize EWA filter weights if needed
	if weightLut == nil {
		weightLut = make([]float32, WEIGHT_LUT_SIZE, WEIGHT_LUT_SIZE)
		for i := 0; i < WEIGHT_LUT_SIZE; i++ {
			alpha := 2.0
			r2 := float64(i) / float64(WEIGHT_LUT_SIZE-1)
			weightLut[i] = float32(math.Exp(-alpha*r2) - math.Exp(-alpha))
		}
	}
}

var weightLut []float32

func NewMIPMapFloat(sres, tres int, img []float32, doTri bool, maxAniso float64, wrap WrapMode) *MIPMapFloat {
	mip := new(MIPMapFloat)
	mip.doTrilinear = doTri
	mip.maxAnisotropy = maxAniso
	mip.wrapMode = wrap

	var resampledImage []float32

	if !IsPowerOf2(sres) || !IsPowerOf2(tres) {
		// Resample image to power-of-two resolution
		sPow2, tPow2 := int(RoundUpPow2(uint32(sres))), int(RoundUpPow2(uint32(tres)))

		// Resample image in $s$ direction
		sWeights := resampleWeights(sres, sPow2)
		resampledImage = make([]float32, sPow2*tPow2, sPow2*tPow2)

		// Apply _sWeights_ to zoom in $s$ direction
		for t := 0; t < tres; t++ {
			for s := 0; s < sPow2; s++ {
				// Compute texel $(s,t)$ in $s$-zoomed image
				resampledImage[t*sPow2+s] = 0.0
				for j := 0; j < 4; j++ {
					origS := sWeights[s].firstTexel + j
					if mip.wrapMode == TEXTURE_REPEAT {
						origS = Mod(origS, sres)
					} else if mip.wrapMode == TEXTURE_CLAMP {
						origS = Clampi(origS, 0, sres-1)
					}
					if origS >= 0 && origS < sres {
						resampledImage[t*sPow2+s] += sWeights[s].weight[j] * img[t*sres+origS]
					}
				}
			}
		}

		// Resample image in $t$ direction
		tWeights := resampleWeights(tres, tPow2)
		workData := make([]float32, tPow2, tPow2)
		for s := 0; s < sPow2; s++ {
			for t := 0; t < tPow2; t++ {
				workData[t] = 0.0
				for j := 0; j < 4; j++ {
					offset := tWeights[t].firstTexel + j
					if mip.wrapMode == TEXTURE_REPEAT {
						offset = Mod(offset, tres)
					} else if mip.wrapMode == TEXTURE_CLAMP {
						offset = Clampi(offset, 0, tres-1)
					}
					if offset >= 0 && offset < tres {
						workData[t] += tWeights[t].weight[j] * resampledImage[offset*sPow2+s]
					}
				}
			}
			for t := 0; t < tPow2; t++ {
				resampledImage[t*sPow2+s] = mip.clamp(workData[t])
			}
		}

		img = resampledImage
		sres = sPow2
		tres = tPow2
	}

	mip.width = sres
	mip.height = tres
	// Initialize levels of MIPMap from image
	mip.nLevels = 1 + Log2Int(float64(Maxi(sres, tres)))
	mip.pyramid = make([][]float32, mip.nLevels, mip.nLevels)
	mip.pyramidSizes = make([][2]int, mip.nLevels, mip.nLevels)

	// Initialize most detailed level of MIPMap
	mip.pyramid[0] = make([]float32, sres*tres, sres*tres)
	copy(mip.pyramid[0], img)
	mip.pyramidSizes[0] = [2]int{sres, tres}

	sRes, tRes := sres, tres

	for i := 1; i < mip.nLevels; i++ {
		// Initialize $i$th MIPMap level from $i-1$st level
		sRes = Maxi(1, sRes/2)
		tRes = Maxi(1, tRes/2)
		mip.pyramid[i] = make([]float32, sRes*tRes, sRes*tRes)
		mip.pyramidSizes[i] = [2]int{sRes, tRes}

		// Filter four texels from finer level of pyramid
		for t := 0; t < tRes; t++ {
			for s := 0; s < sRes; s++ {
				mip.pyramid[i][t*sRes+s] = 0.25 *
					(mip.Texel(i-1, 2*s, 2*t) + mip.Texel(i-1, 2*s+1, 2*t) +
						mip.Texel(i-1, 2*s, 2*t+1) + mip.Texel(i-1, 2*s+1, 2*t+1))
			}
		}
	}
	return mip
}

func (mip *MIPMapFloat) Texel(level, s, t int) float32 {
	//Assert(level < nLevels);
	l := mip.pyramid[level]
	// Compute texel $(s,t)$ accounting for boundary conditions
	switch mip.wrapMode {
	case TEXTURE_REPEAT:
		s = Mod(s, mip.pyramidSizes[level][0])
		t = Mod(t, mip.pyramidSizes[level][1])
	case TEXTURE_CLAMP:
		s = Clampi(s, 0, mip.pyramidSizes[level][0]-1)
		t = Clampi(t, 0, mip.pyramidSizes[level][1]-1)
	case TEXTURE_BLACK:
		{
			if s < 0 || s >= mip.pyramidSizes[level][0] ||
				t < 0 || t >= mip.pyramidSizes[level][1] {
				return 0.0
			}
		}
	}
	//PBRT_ACCESSED_TEXEL(const_cast<MIPMap<T> *>(this), level, s, t);
	return l[s+t*mip.pyramidSizes[level][0]]
}

func (mip *MIPMapFloat) Lookup(s, t, width float64) float32 {
	// Compute MIPMap level for trilinear filtering
	level := float64(mip.nLevels-1) + math.Log2(math.Max(width, 1.0e-8))

	// Perform trilinear interpolation at appropriate MIPMap level
	//PBRT_MIPMAP_TRILINEAR_FILTER(const_cast<MIPMap<T> *>(this), s, t, width, level, nLevels);
	if level < 0 {
		return mip.triangle(0, s, t)
	} else if level >= float64(mip.nLevels-1) {
		return mip.Texel(mip.nLevels-1, 0, 0)
	} else {
		iLevel := Floor2Int(level)
		delta := float32(level) - float32(iLevel)
		return (1.0-delta)*mip.triangle(iLevel, s, t) + delta*mip.triangle(iLevel+1, s, t)
	}
}
func (mip *MIPMapFloat) LookupEwa(s, t, ds0, dt0, ds1, dt1 float64) float32 {
	if mip.doTrilinear {
		//PBRT_STARTED_TRILINEAR_TEXTURE_LOOKUP(s, t);
		val := mip.Lookup(s, t, 2.0*math.Max(math.Max(math.Abs(ds0), math.Abs(dt0)), math.Max(math.Abs(ds1), math.Abs(dt1))))
		//PBRT_FINISHED_TRILINEAR_TEXTURE_LOOKUP();
		return val
	}
	//PBRT_STARTED_EWA_TEXTURE_LOOKUP(s, t);
	// Compute ellipse minor and major axes
	if ds0*ds0+dt0*dt0 < ds1*ds1+dt1*dt1 {
		ds0, ds1 = ds1, ds0
		dt0, dt1 = dt1, dt0
	}
	majorLength := math.Sqrt(ds0*ds0 + dt0*dt0)
	minorLength := math.Sqrt(ds1*ds1 + dt1*dt1)

	// Clamp ellipse eccentricity if too large
	if minorLength*mip.maxAnisotropy < majorLength && minorLength > 0.0 {
		scale := majorLength / (minorLength * mip.maxAnisotropy)
		ds1 *= scale
		dt1 *= scale
		minorLength *= scale
	}
	if minorLength == 0.0 {
		//PBRT_FINISHED_EWA_TEXTURE_LOOKUP();
		//PBRT_STARTED_TRILINEAR_TEXTURE_LOOKUP(s, t);
		val := mip.triangle(0, s, t)
		//PBRT_FINISHED_TRILINEAR_TEXTURE_LOOKUP();
		return val
	}

	// Choose level of detail for EWA lookup and perform EWA filtering
	lod := math.Max(0.0, float64(mip.nLevels)-1.0+math.Log2(minorLength))
	ilod := Floor2Int(lod)
	//PBRT_MIPMAP_EWA_FILTER(const_cast<MIPMap<T> *>(this), s, t, ds0, ds1, dt0, dt1, minorLength, majorLength, lod, nLevels);
	d := float32(lod) - float32(ilod)
	val := (1.0-d)*mip.ewa(ilod, s, t, ds0, dt0, ds1, dt1) + d*mip.ewa(ilod+1, s, t, ds0, dt0, ds1, dt1)
	//PBRT_FINISHED_EWA_TEXTURE_LOOKUP();
	return val
}

func (mip *MIPMapFloat) triangle(level int, s, t float64) float32 {
	level = Clampi(level, 0, mip.nLevels-1)
	s = s*float64(mip.pyramidSizes[level][0]) - 0.5
	t = t*float64(mip.pyramidSizes[level][1]) - 0.5
	s0, t0 := Floor2Int(s), Floor2Int(t)
	ds, dt := float32(s)-float32(s0), float32(t)-float32(t0)
	return (1.0-ds)*(1.0-dt)*mip.Texel(level, s0, t0) +
		(1.0-ds)*dt*mip.Texel(level, s0, t0+1) +
		ds*(1.0-dt)*mip.Texel(level, s0+1, t0) +
		ds*dt*mip.Texel(level, s0+1, t0+1)
}

func (mip *MIPMapFloat) ewa(level int, s, t, ds0, dt0, ds1, dt1 float64) float32 {
	if level >= mip.nLevels {
		return mip.Texel(mip.nLevels-1, 0, 0)
	}

	// Convert EWA coordinates to appropriate scale for level
	s = s*float64(mip.pyramidSizes[level][0]) - 0.5
	t = t*float64(mip.pyramidSizes[level][1]) - 0.5
	ds0 *= float64(mip.pyramidSizes[level][0])
	dt0 *= float64(mip.pyramidSizes[level][1])
	ds1 *= float64(mip.pyramidSizes[level][0])
	dt1 *= float64(mip.pyramidSizes[level][1])

	// Compute ellipse coefficients to bound EWA filter region
	A := dt0*dt0 + dt1*dt1 + 1
	B := -2.0 * (ds0*dt0 + ds1*dt1)
	C := ds0*ds0 + ds1*ds1 + 1
	invF := 1.0 / (A*C - B*B*0.25)
	A *= invF
	B *= invF
	C *= invF

	// Compute the ellipse's $(s,t)$ bounding box in texture space
	det := -B*B + 4.0*A*C
	invDet := 1.0 / det
	uSqrt, vSqrt := math.Sqrt(det*C), math.Sqrt(A*det)
	s0 := Ceil2Int(s - 2.0*invDet*uSqrt)
	s1 := Floor2Int(s + 2.0*invDet*uSqrt)
	t0 := Ceil2Int(t - 2.0*invDet*vSqrt)
	t1 := Floor2Int(t + 2.0*invDet*vSqrt)

	// Scan over ellipse bound and compute quadratic equation
	sum := 0.0
	sumWts := 0.0
	for it := t0; it <= t1; it++ {
		tt := float64(it) - t
		for is := s0; is <= s1; is++ {
			ss := float64(is) - s
			// Compute squared radius and filter texel if inside ellipse
			r2 := A*ss*ss + B*ss*tt + C*tt*tt
			if r2 < 1.0 {
				weight := weightLut[Mini(Float2Int(r2*WEIGHT_LUT_SIZE), WEIGHT_LUT_SIZE-1)]
				sum += float64(mip.Texel(level, is, it) * weight)
				sumWts += float64(weight)
			}
		}
	}
	return float32(sum / sumWts)
}

func (mip *MIPMapFloat) clamp(v float32) float32 { return float32(Clamp(float64(v), 0.0, INFINITYF)) }

func resampleWeights(oldres, newres int) []ResampleWeight {
	//Assert(newres >= oldres);
	wt := make([]ResampleWeight, newres, newres)
	filterwidth := 2.0
	for i := 0; i < newres; i++ {
		// Compute image resampling weights for _i_th texel
		center := (float64(i) + 0.5) * float64(oldres) / float64(newres)
		wt[i].firstTexel = Floor2Int((center - filterwidth) + 0.5)
		for j := 0; j < 4; j++ {
			pos := float64(wt[i].firstTexel+j) + 0.5
			wt[i].weight[j] = float32(Lanczos((pos-center)/filterwidth, 2.0))
		}

		// Normalize filter weights for texel resampling
		invSumWts := 1.0 / (wt[i].weight[0] + wt[i].weight[1] + wt[i].weight[2] + wt[i].weight[3])
		for j := 0; j < 4; j++ {
			wt[i].weight[j] *= invSumWts
		}
	}
	return wt
}

func NewMIPMapSpectrum(sres, tres int, img []Spectrum, doTri bool, maxAniso float64, wrap WrapMode) *MIPMapSpectrum {
	mip := new(MIPMapSpectrum)
	mip.doTrilinear = doTri
	mip.maxAnisotropy = maxAniso
	mip.wrapMode = wrap

	var resampledImage []Spectrum

	if !IsPowerOf2(sres) || !IsPowerOf2(tres) {
		// Resample image to power-of-two resolution
		sPow2, tPow2 := int(RoundUpPow2(uint32(sres))), int(RoundUpPow2(uint32(tres)))

		// Resample image in $s$ direction
		sWeights := resampleWeights(sres, sPow2)
		resampledImage = make([]Spectrum, sPow2*tPow2, sPow2*tPow2)

		// Apply _sWeights_ to zoom in $s$ direction
		for t := 0; t < tres; t++ {
			for s := 0; s < sPow2; s++ {
				// Compute texel $(s,t)$ in $s$-zoomed image
				resampledImage[t*sPow2+s] = *CreateSpectrum1(0.0)
				for j := 0; j < 4; j++ {
					origS := sWeights[s].firstTexel + j
					if mip.wrapMode == TEXTURE_REPEAT {
						origS = Mod(origS, sres)
					} else if mip.wrapMode == TEXTURE_CLAMP {
						origS = Clampi(origS, 0, sres-1)
					}
					if origS >= 0 && origS < sres {
						resampledImage[t*sPow2+s] = *resampledImage[t*sPow2+s].Add(img[t*sres+origS].Scale(sWeights[s].weight[j]))
					}
				}
			}
		}

		// Resample image in $t$ direction
		tWeights := resampleWeights(tres, tPow2)
		workData := make([]Spectrum, tPow2, tPow2)
		for s := 0; s < sPow2; s++ {
			for t := 0; t < tPow2; t++ {
				workData[t] = *CreateSpectrum1(0.0)
				for j := 0; j < 4; j++ {
					offset := tWeights[t].firstTexel + j
					if mip.wrapMode == TEXTURE_REPEAT {
						offset = Mod(offset, tres)
					} else if mip.wrapMode == TEXTURE_CLAMP {
						offset = Clampi(offset, 0, tres-1)
					}
					if offset >= 0 && offset < tres {
						workData[t] = *workData[t].Add(resampledImage[offset*sPow2+s].Scale(tWeights[t].weight[j]))
					}
				}
			}
			for t := 0; t < tPow2; t++ {
				resampledImage[t*sPow2+s] = mip.clamp(workData[t])
			}
		}

		img = resampledImage
		sres = sPow2
		tres = tPow2
	}

	mip.width = sres
	mip.height = tres
	// Initialize levels of MIPMap from image
	mip.nLevels = 1 + Log2Int(float64(Maxi(sres, tres)))
	mip.pyramid = make([][]Spectrum, mip.nLevels, mip.nLevels)
	mip.pyramidSizes = make([][2]int, mip.nLevels, mip.nLevels)

	// Initialize most detailed level of MIPMap
	mip.pyramid[0] = make([]Spectrum, sres*tres, sres*tres)
	copy(mip.pyramid[0], img)
	mip.pyramidSizes[0] = [2]int{sres, tres}

	sRes, tRes := sres, tres

	for i := 1; i < mip.nLevels; i++ {
		// Initialize $i$th MIPMap level from $i-1$st level
		sRes = Maxi(1, sRes/2)
		tRes = Maxi(1, tRes/2)
		mip.pyramid[i] = make([]Spectrum, sRes*tRes, sRes*tRes)
		mip.pyramidSizes[i] = [2]int{sRes, tRes}

		// Filter four texels from finer level of pyramid
		for t := 0; t < tRes; t++ {
			for s := 0; s < sRes; s++ {
				mip.pyramid[i][t*sRes+s] = *(mip.Texel(i-1, 2*s, 2*t).Add(mip.Texel(i-1, 2*s+1, 2*t)).Add(mip.Texel(i-1, 2*s, 2*t+1).Add(mip.Texel(i-1, 2*s+1, 2*t+1)))).Scale(0.25)
			}
		}
	}
	return mip
}

func (mip *MIPMapSpectrum) Texel(level, s, t int) *Spectrum {
	//Assert(level < nLevels);
	l := mip.pyramid[level]
	// Compute texel $(s,t)$ accounting for boundary conditions
	switch mip.wrapMode {
	case TEXTURE_REPEAT:
		s = Mod(s, mip.pyramidSizes[level][0])
		t = Mod(t, mip.pyramidSizes[level][1])
	case TEXTURE_CLAMP:
		s = Clampi(s, 0, mip.pyramidSizes[level][0]-1)
		t = Clampi(t, 0, mip.pyramidSizes[level][1]-1)
	case TEXTURE_BLACK:
		{
			if s < 0 || s >= mip.pyramidSizes[level][0] ||
				t < 0 || t >= mip.pyramidSizes[level][1] {
				return CreateSpectrum1(0.0)
			}
		}
	}
	//PBRT_ACCESSED_TEXEL(const_cast<MIPMap<T> *>(this), level, s, t);
	return &l[s+t*mip.pyramidSizes[level][0]]
}

func (mip *MIPMapSpectrum) Lookup(s, t, width float64) *Spectrum {
	// Compute MIPMap level for trilinear filtering
	level := float64(mip.nLevels-1) + math.Log2(math.Max(width, 1.0e-8))

	// Perform trilinear interpolation at appropriate MIPMap level
	//PBRT_MIPMAP_TRILINEAR_FILTER(const_cast<MIPMap<T> *>(this), s, t, width, level, nLevels);
	if level < 0 {
		return mip.triangle(0, s, t)
	} else if level >= float64(mip.nLevels-1) {
		return mip.Texel(mip.nLevels-1, 0, 0)
	} else {
		iLevel := Floor2Int(level)
		delta := float32(level) - float32(iLevel)
		return mip.triangle(iLevel, s, t).Scale(1.0 - delta).Add(mip.triangle(iLevel+1, s, t).Scale(delta))
	}
}
func (mip *MIPMapSpectrum) LookupEwa(s, t, ds0, dt0, ds1, dt1 float64) *Spectrum {
	if mip.doTrilinear {
		//PBRT_STARTED_TRILINEAR_TEXTURE_LOOKUP(s, t);
		val := mip.Lookup(s, t, 2.0*math.Max(math.Max(math.Abs(ds0), math.Abs(dt0)), math.Max(math.Abs(ds1), math.Abs(dt1))))
		//PBRT_FINISHED_TRILINEAR_TEXTURE_LOOKUP();
		return val
	}
	//PBRT_STARTED_EWA_TEXTURE_LOOKUP(s, t);
	// Compute ellipse minor and major axes
	if ds0*ds0+dt0*dt0 < ds1*ds1+dt1*dt1 {
		ds0, ds1 = ds1, ds0
		dt0, dt1 = dt1, dt0
	}
	majorLength := math.Sqrt(ds0*ds0 + dt0*dt0)
	minorLength := math.Sqrt(ds1*ds1 + dt1*dt1)

	// Clamp ellipse eccentricity if too large
	if minorLength*mip.maxAnisotropy < majorLength && minorLength > 0.0 {
		scale := majorLength / (minorLength * mip.maxAnisotropy)
		ds1 *= scale
		dt1 *= scale
		minorLength *= scale
	}
	if minorLength == 0.0 {
		//PBRT_FINISHED_EWA_TEXTURE_LOOKUP();
		//PBRT_STARTED_TRILINEAR_TEXTURE_LOOKUP(s, t);
		val := mip.triangle(0, s, t)
		//PBRT_FINISHED_TRILINEAR_TEXTURE_LOOKUP();
		return val
	}

	// Choose level of detail for EWA lookup and perform EWA filtering
	lod := math.Max(0.0, float64(mip.nLevels)-1.0+math.Log2(minorLength))
	ilod := Floor2Int(lod)
	//PBRT_MIPMAP_EWA_FILTER(const_cast<MIPMap<T> *>(this), s, t, ds0, ds1, dt0, dt1, minorLength, majorLength, lod, nLevels);
	d := float32(lod) - float32(ilod)
	val := mip.ewa(ilod, s, t, ds0, dt0, ds1, dt1).Scale(1.0 - d).Add(mip.ewa(ilod+1, s, t, ds0, dt0, ds1, dt1).Scale(d))
	//PBRT_FINISHED_EWA_TEXTURE_LOOKUP();
	return val
}

func (mip *MIPMapSpectrum) clamp(v Spectrum) Spectrum { return *v.Clamp(0.0, INFINITYF) }

func (mip *MIPMapSpectrum) triangle(level int, s, t float64) *Spectrum {
	level = Clampi(level, 0, mip.nLevels-1)
	s = s*float64(mip.pyramidSizes[level][0]) - 0.5
	t = t*float64(mip.pyramidSizes[level][1]) - 0.5
	s0, t0 := Floor2Int(s), Floor2Int(t)
	ds, dt := float32(s)-float32(s0), float32(t)-float32(t0)
	return mip.Texel(level, s0, t0).Scale((1.0 - ds) * (1.0 - dt)).Add(mip.Texel(level, s0, t0+1).Scale((1.0 - ds) * dt)).Add(mip.Texel(level, s0+1, t0).Scale(ds * (1.0 - dt))).Add(mip.Texel(level, s0+1, t0+1).Scale(ds * dt))
}

func (mip *MIPMapSpectrum) ewa(level int, s, t, ds0, dt0, ds1, dt1 float64) *Spectrum {
	if level >= mip.nLevels {
		return mip.Texel(mip.nLevels-1, 0, 0)
	}

	// Convert EWA coordinates to appropriate scale for level
	s = s*float64(mip.pyramidSizes[level][0]) - 0.5
	t = t*float64(mip.pyramidSizes[level][1]) - 0.5
	ds0 *= float64(mip.pyramidSizes[level][0])
	dt0 *= float64(mip.pyramidSizes[level][1])
	ds1 *= float64(mip.pyramidSizes[level][0])
	dt1 *= float64(mip.pyramidSizes[level][1])

	// Compute ellipse coefficients to bound EWA filter region
	A := dt0*dt0 + dt1*dt1 + 1
	B := -2.0 * (ds0*dt0 + ds1*dt1)
	C := ds0*ds0 + ds1*ds1 + 1
	invF := 1.0 / (A*C - B*B*0.25)
	A *= invF
	B *= invF
	C *= invF

	// Compute the ellipse's $(s,t)$ bounding box in texture space
	det := -B*B + 4.0*A*C
	invDet := 1.0 / det
	uSqrt, vSqrt := math.Sqrt(det*C), math.Sqrt(A*det)
	s0 := Ceil2Int(s - 2.0*invDet*uSqrt)
	s1 := Floor2Int(s + 2.0*invDet*uSqrt)
	t0 := Ceil2Int(t - 2.0*invDet*vSqrt)
	t1 := Floor2Int(t + 2.0*invDet*vSqrt)

	// Scan over ellipse bound and compute quadratic equation
	sum := *CreateSpectrum1(0.0)
	sumWts := 0.0
	for it := t0; it <= t1; it++ {
		tt := float64(it) - t
		for is := s0; is <= s1; is++ {
			ss := float64(is) - s
			// Compute squared radius and filter texel if inside ellipse
			r2 := A*ss*ss + B*ss*tt + C*tt*tt
			if r2 < 1.0 {
				weight := weightLut[Mini(Float2Int(r2*WEIGHT_LUT_SIZE), WEIGHT_LUT_SIZE-1)]
				sum = *sum.Add(mip.Texel(level, is, it).Scale(weight))
				sumWts += float64(weight)
			}
		}
	}
	return sum.InvScale(float32(sumWts))
}
