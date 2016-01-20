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
    weight [4]float64
}

type MIPMapData struct {
    doTrilinear bool 
    maxAnisotropy float64
    wrapMode WrapMode
    width, height, nLevels int    
}
type MIPMapFloat struct {
   MIPMapData
    pyramid [][]float64
    }

type MIPMapSpectrum struct {
    MIPMapData
    pyramid [][]Spectrum
}

func init() {
    // Initialize EWA filter weights if needed
    if weightLut == nil {
        weightLut = make([]float64, WEIGHT_LUT_SIZE, WEIGHT_LUT_SIZE)
        for i := 0; i < WEIGHT_LUT_SIZE; i++ {
            alpha := 2.0
            r2 := float64(i) / float64(WEIGHT_LUT_SIZE - 1)
            weightLut[i] = math.Exp(-alpha * r2) - math.Exp(-alpha)
        }
    }
}

var weightLut []float64

func NewMIPMapFloat(sres, tres int, img []float64, doTri bool, maxAniso float64, wrap WrapMode) *MIPMapFloat {
    mip := new(MIPMapFloat)
    mip.doTrilinear = doTri
    mip.maxAnisotropy = maxAniso
    mip.wrapMode = wrap
    
    var resampledImage []float64
    
    if !IsPowerOf2(sres) || !IsPowerOf2(tres) {
        // Resample image to power-of-two resolution
        sPow2, tPow2 := int(RoundUpPow2(uint32(sres))), int(RoundUpPow2(uint32(tres)))

        // Resample image in $s$ direction
        sWeights := resampleWeights(sres, sPow2)
        resampledImage = make([]float64, sPow2 * tPow2, sPow2 * tPow2)

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
                        resampledImage[t*sPow2+s] += sWeights[s].weight[j] *
                                                     img[t*sres + origS]
                    }                                 
                }
            }
        }
        
        // Resample image in $t$ direction
        tWeights := resampleWeights(tres, tPow2)
        workData := make([]float64, tPow2, tPow2)
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
                        workData[t] += tWeights[t].weight[j] * resampledImage[offset*sPow2 + s]
                    }        
                }
            }
            for t := 0; t < tPow2; t++ {
                resampledImage[t*sPow2 + s] = mip.clamp(workData[t])
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
    mip.pyramid = make([][]float64, mip.nLevels, mip.nLevels)

    // Initialize most detailed level of MIPMap
    mip.pyramid[0] = make([]float64, sres * tres, sres * tres)
    copy(mip.pyramid[0], img)
    
    sRes, tRes := sres, tres
    
    for i := 1; i < mip.nLevels; i++ {
        // Initialize $i$th MIPMap level from $i-1$st level
         sRes = Maxi(1, sRes/2)
         tRes = Maxi(1, tRes/2)
        mip.pyramid[i] = make([]float64, sRes * tRes, sRes * tRes)

        // Filter four texels from finer level of pyramid
        for t := 0; t < tRes; t++ {
            for s := 0; s < sRes; s++ {
                mip.pyramid[i][t*sRes+s] = 0.25 *
                   (mip.Texel(i-1, 2*s, 2*t)   + mip.Texel(i-1, 2*s+1, 2*t) +
                    mip.Texel(i-1, 2*s, 2*t+1) + mip.Texel(i-1, 2*s+1, 2*t+1))
            }
        }
    }
    return mip
}

func (mip *MIPMapFloat) Texel(level, s, t int) float64 {
    return 0.0   
}
func (mip *MIPMapFloat) Lookup(s, t, width float64) float64 {
    return 0.0
}
func (mip *MIPMapFloat) LookupEwa(s, t, ds0, dt0, ds1, dt1 float64) float64 {
    return 0.0
}

func (mip *MIPMapFloat) triangle(level int, s, t float64) float64 {
    return 0.0
}
func (mip *MIPMapFloat) ewa(level int, s, t, ds0, dt0, ds1, dt1 float64) float64 {
    return 0.0    
}
func (mip *MIPMapFloat) clamp(v float64) float64 { return Clamp(v, 0.0, INFINITY) }

func resampleWeights(oldres, newres int) []ResampleWeight {
   //Assert(newres >= oldres);
   wt := make([]ResampleWeight, newres, newres)
   filterwidth := 2.0
   for i := 0; i < newres; i++ {
       // Compute image resampling weights for _i_th texel
       center := (float64(i) + 0.5) * float64(oldres) / float64(newres)
       wt[i].firstTexel = Floor2Int((center - filterwidth) + 0.5)
        for j := 0; j < 4; j++ {
            pos := float64(wt[i].firstTexel + j) + 0.5
            wt[i].weight[j] = Lanczos((pos - center) / filterwidth, 2.0)
        }

        // Normalize filter weights for texel resampling
        invSumWts := 1.0 / (wt[i].weight[0] + wt[i].weight[1] + wt[i].weight[2] + wt[i].weight[3])
        for j := 0; j < 4; j++ {
            wt[i].weight[j] *= invSumWts
        }
    }
    return wt
}
    
func NewMIPMapSpectrum(xres, yres int, pixels []Spectrum, doTri bool, maxAniso float64, wrap WrapMode) *MIPMapSpectrum {
    return nil
}

func (mip *MIPMapSpectrum) Texel(level, s, t int) *Spectrum {
    return nil
}
func (mip *MIPMapSpectrum) Lookup(s, t, width float64) Spectrum {
    return *CreateSpectrum1(0.0)
}
func (mip *MIPMapSpectrum) LookupEwa(s, t, ds0, dt0, ds1, dt1 float64) Spectrum {
    return *CreateSpectrum1(0.0)    
}
