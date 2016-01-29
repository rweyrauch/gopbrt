package core

import  (
	"math"
)

// Spherical Harmonics Local Definitions
func legendrep(x float64, lmax int, out []float64) {
	P := func(l,m int) float64 { return out[SHIndex(l,m)] }
	
    // Compute $m=0$ Legendre values using recurrence
    out[SHIndex(0,0)] = 1.0
    out[SHIndex(1,0)] = x
    for l := 2; l <= lmax; l++ {
        out[SHIndex(l,0)] = (float64(2*l-1)*x*P(l-1,0) - float64(l-1)*P(l-2,0)) / float64(l)
        Assert(!math.IsNaN(P(l, 0)))
        Assert(!math.IsInf(P(l, 0), 0))
    }

    // Compute $m=l$ edge using Legendre recurrence
     neg := -1.0
     dfact := 1.0
     xroot := math.Sqrt(math.Max(0.0, 1.0 - x*x))
     xpow := xroot
    for l := 1; l <= lmax; l++ {
        out[SHIndex(l,l)] = neg * dfact * xpow
        Assert(!math.IsNaN(P(l, l)))
        Assert(!math.IsInf(P(l, l), 0))
        neg *= -1.0     // neg = (-1)^l
        dfact *= float64(2*l + 1) // dfact = (2*l-1)!!
        xpow *= xroot    // xpow = powf(1.f - x*x, float(l) * 0.5f);
    }

    // Compute $m=l-1$ edge using Legendre recurrence
    for l := 2; l <= lmax; l++ {
        out[SHIndex(l,l-1)] = x * float64(2*l-1) * P(l-1, l-1)
        Assert(!math.IsNaN(P(l, l-1)))
        Assert(!math.IsInf(P(l, l-1), 0))
    }

    // Compute $m=1, \ldots, l-2$ values using Legendre recurrence
    for l := 3; l <= lmax; l++ {
        for m := 1; m <= l-2; m++ {
            out[SHIndex(l,m)] = (float64(2 * (l-1) + 1) * x * P(l-1,m) - float64(l-1+m) * P(l-2,m)) / float64(l - m)
            Assert(!math.IsNaN(P(l, m)))
            Assert(!math.IsInf(P(l, m), 0))
        }
     }   
}

const (
	INV_FOURPI = 1.0 / (4.0 * math.Pi)
)

func K(l, m int) float64 {
    return math.Sqrt((2.0 * float64(l) + 1.0) * INV_FOURPI * divfact(l, m))
}

func divfact(a, b int) float64 {
    if b == 0 { return 1.0 }
    fa := float64(a)
    fb := math.Abs(float64(b))
    v := 1.0
    for x := fa-fb+1.0; x <= fa+fb; x += 1.0 {
        v *= x
    }    
    return 1.0 / v
}

// n!! = 1 if n==0 or 1, otherwise n * (n-2)!!
func dfact(v float64) float64 {
    if v <= 1.0 { return 1.0 }
    return v * dfact(v - 2.0)
}

func fact(v float64) float64 {
    if v <= 1.0 { return 1.0 }
    return v * fact(v - 1.0)
}

func sinCosIndexed(s, c float64, n int, sout, cout []float64) {
    si, ci := 0.0, 1.0
    for i := 0; i < n; i++ {
        // Compute $\sin{}i\phi$ and $\cos{}i\phi$ using recurrence
        sout[i] = si
        cout[i] = ci
        oldsi := si
        si = si * c + ci * s
        ci = ci * c - oldsi * s
    }
}

func toZYZ(m *Matrix4x4) (alpha, beta, gamma float64) {
	 M := func(a, b int) float64 { return m.m[a][b] }

    sy := math.Sqrt(M(2,1)*M(2,1) + M(2,0)*M(2,0))
    if sy > 16.0*math.Nextafter(0.0, 1.0) {
        gamma = -math.Atan2(M(1,2), -M(0,2))
        beta  = -math.Atan2(sy, M(2,2))
        alpha = -math.Atan2(M(2,1), M(2,0))
    } else {
        gamma =  0.0
        beta  = -math.Atan2(sy, M(2,2))
        alpha = -math.Atan2(-M(1,0), M(1,1))
    }
    return alpha, beta, gamma
}


func lambda(l float64) float64 {
    return math.Sqrt((4.0 * math.Pi) / (2.0 * l + 1.0))
}



// Spherical Harmonics Declarations
func SHTerms(lmax int) int {
    return (lmax + 1) * (lmax + 1)
}

func SHIndex(l, m int) int {
    return l*l + l + m
}

func SHEvaluate(w *Vector, lmax int, out []float64) {
    // Compute Legendre polynomial values for $\cos\theta$
    Assert(w.Length() > 0.995 && w.Length() < 1.005)
    legendrep(w.z, lmax, out)

    // Compute $K_l^m$ coefficients
    Klm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    for l := 0; l <= lmax; l++ {
        for m := -l; m <= l; m++ {
            Klm[SHIndex(l, m)] = K(l, m)
        }    
	}
    
    // Compute $\sin\phi$ and $\cos\phi$ values
    sins := make([]float64, lmax+1, lmax+1) 
    coss := make([]float64, lmax+1, lmax+1)
    xyLen := math.Sqrt(math.Max(0.0, 1.0 - w.z*w.z))
    if xyLen == 0.0 {
        for i := 0; i <= lmax; i++ { sins[i] = 0.0 }
        for i := 0; i <= lmax; i++ { coss[i] = 1.0 }
    } else {
        sinCosIndexed(w.y / xyLen, w.x / xyLen, lmax+1, sins, coss)
	}
    
    // Apply SH definitions to compute final $(l,m)$ values
    sqrt2 := math.Sqrt(2.0)
    for l := 0; l <= lmax; l++ {
        for m := -l; m < 0; m++ {
            out[SHIndex(l, m)] = sqrt2 * Klm[SHIndex(l, m)] * out[SHIndex(l, -m)] * sins[-m]
            Assert(!math.IsNaN(out[SHIndex(l,m)]))
            Assert(!math.IsInf(out[SHIndex(l,m)], 0))
        }
        out[SHIndex(l, 0)] *= Klm[SHIndex(l, 0)]
        for m := 1; m <= l; m++ {
            out[SHIndex(l, m)] *= sqrt2 * Klm[SHIndex(l, m)] * coss[m]
            Assert(!math.IsNaN(out[SHIndex(l,m)]))
            Assert(!math.IsInf(out[SHIndex(l,m)], 0))
        }
    }
}

func SHWriteImage(filename string, c []Spectrum, lmax, yres int) {
    xres := 2 * yres
    rgb := make([]float32, xres * yres * 3, xres * yres * 3)
    Ylm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    for y := 0; y < yres; y++ {
        theta := (float64(y) + 0.5) / float64(yres) * math.Pi
        for x := 0; x < xres; x++ {
            phi := (float64(x) + 0.5) / float64(xres) * 2.0 * math.Pi
            // Compute RGB color for direction $(\theta,\phi)$ from SH coefficients
            w := SphericalDirection(math.Sin(theta), math.Cos(theta), phi)
            SHEvaluate(w, lmax, Ylm)
            val := NewSpectrum1(0.0)
            for i := 0; i < SHTerms(lmax); i++ {
                val = val.Add(c[i].Scale(Ylm[i]))
			}                
            offset := xres * y + x
            rgb[3*offset+0] = float32(val.c[0])
            rgb[3*offset+1] = float32(val.c[1])
            rgb[3*offset+2] = float32(val.c[2])            
        }
    }

    WriteImage(filename, rgb, nil, xres, yres, xres, yres, 0, 0)
}


func SHProjectIncidentDirectRadiance(p *Point, pEpsilon, time float64, arena *MemoryArena, scene *Scene, computeLightVis bool, lmax int, rng *RNG, c_d []Spectrum) {
    // Loop over light sources and sum their SH coefficients
    c := make([]Spectrum, SHTerms(lmax), SHTerms(lmax))
    for i := 0; i < len(scene.lights); i++ {
        light := scene.lights[i];
        c = light.SHProject(p, pEpsilon, lmax, scene, computeLightVis, time, rng)
        for j := 0; j < SHTerms(lmax); j++ {
            c_d[j] = *c_d[j].Add(&c[j])
		}            
    }
    SHReduceRinging(c_d, lmax, 0.005)
}


func SHProjectIncidentIndirectRadiance(p *Point, pEpsilon, time float64, renderer Renderer, origSample *Sample, scene *Scene, lmax int, rng *RNG, ns int, c_i []Spectrum) {
    sample := origSample.Duplicate(1)
    var arena *MemoryArena
    scramble := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    nSamples := RoundUpPow2(uint32(ns))
    Ylm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    var i uint32
    for i = 0; i < nSamples; i++ {
        // Sample incident direction for radiance probe
        u := [2]float64{0.0, 0.0}
        Sample02(i, scramble, u[:])
        wi := UniformSampleSphere(u[0], u[1])
        pdf := UniformSpherePdf()

        // Compute incident radiance along direction for probe
        Li := NewSpectrum1(0.0)
        ray := CreateRayDifferential(p, wi, pEpsilon, INFINITY, time, 0)

        // Fill in values in _sample_ for radiance probe ray
        sample[0].time = time
        for j := 0; j < len(sample[0].n1D); j++ {
            for k := 0; k < sample[0].n1D[j]; k++ {
                sample[0].oneD[j][k] = rng.RandomFloat()
			}
        }    	                
        for j := 0; j < len(sample[0].n2D); j++ {
            for k := 0; k < 2 * sample[0].n2D[j]; k++ {
                sample[0].twoD[j][k] = rng.RandomFloat()
			}
		}            	                
        Li, _, _ = renderer.Li(scene, ray, &sample[0], rng, arena)

        // Update SH coefficients for probe sample point
        SHEvaluate(wi, lmax, Ylm)
        for j := 0; j < SHTerms(lmax); j++ {
            c_i[j] = *c_i[j].Add(Li.Scale(Ylm[j] / pdf * float64(nSamples)))
		}            
    }
}


func SHReduceRinging(c []Spectrum, lmax int, lambda float64) {
    for l := 0; l <= lmax; l++ {
        scale := 1.0 / (1.0 + lambda * float64(l * l * (l + 1) * (l + 1)))
        for m := -l; m <= l; m++ {
            c[SHIndex(l, m)] = *c[SHIndex(l, m)].Scale(scale)
        }    
    }
}


func SHRotate(c_in, c_out []Spectrum, m *Matrix4x4, lmax int, arena *MemoryArena) {
    alpha, beta, gamma := toZYZ(m)
    work := make([]Spectrum, SHTerms(lmax), SHTerms(lmax))
    SHRotateZ(c_in, c_out, gamma, lmax)
    SHRotateXPlus(c_out, work, lmax)
    SHRotateZ(work, c_out, beta, lmax)
    SHRotateXMinus(c_out, work, lmax)
    SHRotateZ(work, c_out, alpha, lmax)
}


func SHRotateZ(c_in, c_out []Spectrum, alpha float64, lmax int) {
    //Assert(c_in != c_out)
    c_out[0] = c_in[0]
    if lmax == 0 { return }
    // Precompute sine and cosine terms for $z$-axis SH rotation
    ct := make([]float64, lmax+1, lmax+1)
    st := make([]float64, lmax+1, lmax+1)
    sinCosIndexed(math.Sin(alpha), math.Cos(alpha), lmax+1, st, ct)
    for l := 1; l <= lmax; l++ {
        // Rotate coefficients for band _l_ about $z$
        for m := -l; m < 0; m++ {
            c_out[SHIndex(l, m)] = *(c_in[SHIndex(l,  m)].Scale(ct[-m]).Add(c_in[SHIndex(l, -m)].Scale(-st[-m])))
		}                
        c_out[SHIndex(l, 0)] = c_in[SHIndex(l, 0)]
        for m := 1; m <= l; m++ {
            c_out[SHIndex(l, m)] = *(c_in[SHIndex(l,  m)].Scale(ct[m]).Add(c_in[SHIndex(l, -m)].Scale(st[m])))
		}                
    }
}


var c_costheta [18]float64 = [18]float64{ 
	0.8862268925, 1.0233267546,
        0.4954159260, 0.0000000000, -0.1107783690, 0.0000000000,
        0.0499271341, 0.0000000000, -0.0285469331, 0.0000000000,
        0.0185080823, 0.0000000000, -0.0129818395, 0.0000000000,
        0.0096125342, 0.0000000000, -0.0074057109, 0.0000000000 }

func SHConvolveCosTheta(lmax int, c_in, c_out []Spectrum) {
    for l := 0; l <= lmax; l++ {
        for m := -l; m <= l; m++ {
            o := SHIndex(l, m)
            if l < 18 {
            	c_out[o] = *c_in[o].Scale(lambda(float64(l)) * c_costheta[l])
            } else  {      
            	c_out[o] = *NewSpectrum1(0.0) 
            	}
        }
	}        
}


func SHConvolvePhong(lmax int, n float64, c_in, c_out []Spectrum) {
    for l := 0; l <= lmax; l++ {
        c_phong := math.Exp(-float64(l*l) / (2.0 * float64(n)))
        for m := -l; m <= l; m++ {
            o := SHIndex(l, m)
            c_out[o] = *c_in[o].Scale(lambda(float64(l)) * c_phong)
        }
    }
}


func SHComputeDiffuseTransfer(p *Point, n *Normal, rayEpsilon float64, scene *Scene, rng *RNG, nSamples, lmax int,  c_transfer []Spectrum) {
    for i := 0; i < SHTerms(lmax); i++ {
        c_transfer[i].c[0], c_transfer[i].c[1], c_transfer[i].c[2] = 0.0, 0.0, 0.0
    }    
    scramble := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    Ylm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    var i uint32
    for i = 0; i < uint32(nSamples); i++ {
        // Sample _i_th direction and compute estimate for transfer coefficients
        u := [2]float64{0.0, 0.0}
        Sample02(i, scramble, u[:])
         w := UniformSampleSphere(u[0], u[1])
        pdf := UniformSpherePdf()
        if DotVectorNormal(w, n) > 0.0 && !scene.IntersectP(CreateRay(p, w, rayEpsilon, INFINITY, 0.0, 0)) {
            // Accumulate contribution of direction $\w{}$ to transfer coefficients
            SHEvaluate(w, lmax, Ylm)
            for j := 0; j < SHTerms(lmax); j++ {
                c_transfer[j] = *c_transfer[j].Add(NewSpectrum1((Ylm[j] * AbsDotVectorNormal(w, n)) / (pdf * float64(nSamples))))
            }    
        }
    }
}


func SHComputeTransferMatrix(p *Point, rayEpsilon float64, scene *Scene, rng *RNG, nSamples, lmax int, T []Spectrum) {
    for i := 0; i < SHTerms(lmax)*SHTerms(lmax); i++ {
        T[i].c[0] = 0.0
        T[i].c[1] = 0.0
        T[i].c[2] = 0.0
	}        
    scramble := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    Ylm := make([]float64, SHTerms(lmax), SHTerms(lmax))
    var i uint32
    for i = 0; i < uint32(nSamples); i++ {
        // Compute Monte Carlo estimate of $i$th sample for transfer matrix
        u := [2]float64{0.0, 0.0}
        Sample02(i, scramble, u[:])
         w := UniformSampleSphere(u[0], u[1])
         pdf := UniformSpherePdf()
        if !scene.IntersectP(CreateRay(p, w, rayEpsilon, INFINITY, 0.0, 0)) {
            // Update transfer matrix for unoccluded direction
            SHEvaluate(w, lmax, Ylm)
            for j := 0; j < SHTerms(lmax); j++ {
                for k := 0; k < SHTerms(lmax); k++ {
                    T[j*SHTerms(lmax)+k] = *T[j*SHTerms(lmax)+k].Add(NewSpectrum1((Ylm[j] * Ylm[k]) / (pdf * float64(nSamples))))
				}
			}                                    	
        }
    }
}


func SHComputeBSDFMatrix(Kd, Ks *Spectrum, roughness float64, rng *RNG, nSamples, lmax int, B []Spectrum) {
    for i := 0; i < SHTerms(lmax)*SHTerms(lmax); i++ {
        B[i].c[0] = 0.0
        B[i].c[1] = 0.0
        B[i].c[2] = 0.0
	}        
    // Create _BSDF_ for computing BSDF transfer matrix
    dg := CreateDiffGeometry(CreatePoint(0,0,0), CreateVector(1,0,0), CreateVector(0,1,0), CreateNormal(0,0,0), CreateNormal(0,0,0), 0, 0, nil)
    bsdf := NewBSDF(dg, CreateNormal(0,0,1), 1.0)
    bsdf.Add(NewLambertian(Kd))
    fresnel := &FresnelDielectric{1.50, 1.0}
    bsdf.Add(NewMicrofacet(Ks, fresnel, &Blinn{1.0 / roughness}))

    // Precompute directions $\w{}$ and SH values for directions
    Ylm := make([]float64, SHTerms(lmax) * nSamples, SHTerms(lmax) * nSamples)
    w := make([]Vector, nSamples, nSamples)
    scramble := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    var i uint32
    for i = 0; i < uint32(nSamples); i++ {
        u := [2]float64{0.0, 0.0}
        Sample02(i, scramble, u[:])
        w[i] = *UniformSampleSphere(u[0], u[1])
        SHEvaluate(&w[int(i)], lmax, Ylm[SHTerms(lmax)*int(i):])
    }

    // Compute double spherical integral for BSDF matrix
    for osamp := 0; osamp < nSamples; osamp++ {
        wo := &w[osamp]
        for isamp := 0; isamp < nSamples; isamp++ {
            wi := &w[isamp]
            // Update BSDF matrix elements for sampled directions
            f := bsdf.f(wo, wi, BSDF_ALL)
            if !f.IsBlack() {
                pdf := UniformSpherePdf() * UniformSpherePdf()
                f = f.Scale(math.Abs(CosTheta(wi)) / (pdf * float64(nSamples * nSamples)))
                for i := 0; i < SHTerms(lmax); i++ {
                    for j := 0; j < SHTerms(lmax); j++ {
                        B[i*SHTerms(lmax)+j] = *B[i*SHTerms(lmax)+j].Add(f.Scale(Ylm[isamp*SHTerms(lmax)+j] * Ylm[osamp*SHTerms(lmax)+i]))
					}
				}                    	                                                				
            }
        }
    }
}


func SHMatrixVectorMultiply(M, v []Spectrum, vout []Spectrum, lmax int) {
    for i := 0; i < SHTerms(lmax); i++ {
        vout[i] = *NewSpectrum1(0.0)
        for j := 0; j < SHTerms(lmax); j++ {
            vout[i] = *vout[i].Add(M[SHTerms(lmax) * i + j].Mult(&v[j]))
        }    
    }
}
