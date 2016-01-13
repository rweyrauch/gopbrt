package pbrt

func Min(x, y int) int {
    if x < y {
        return x
    }
    return y
}

func Max(x, y int) int {
    if x > y {
        return x
    }
    return y
}

type Integrator interface {
    Preprocess(scene *Scene, camera *Camera, enderer *Renderer)
    RequestSamples(sampler *Sampler, sample *Sample, scene *Scene)
}

type SurfaceIntegrator interface {
	Integrator
    Li(scene *Scene, renderer *Renderer, ray *RayDifferential, isect *Intersection,
        sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
}

func UniformSampleAllLights(scene *Scene, renderer *Renderer,
    arena *MemoryArena, p *Point, n *Normal, wo *Vector,
    rayEpsilon, time float64, bsdf *BSDF, sample *Sample, rng *RNG,
    lightSampleOffsets []LightSampleOffsets,
    bsdfSampleOffsets []BSDFSampleOffsets) *Spectrum {
    	
    L := CreateSpectrum1(0.0)
    for i, light := range scene.lights {
        nSamples := 1
        if lightSampleOffsets != nil {
        	nSamples = lightSampleOffsets[i].nSamples
        }
        // Estimate direct lighting from _light_ samples
        Ld := CreateSpectrum1(0.0)
        for j := 0; j < nSamples; j++ {
            // Find light and BSDF sample values for direct lighting estimate
            var lightSample *LightSample
            var bsdfSample *BSDFSample
            if lightSampleOffsets != nil && bsdfSampleOffsets != nil {
                lightSample = CreateLightSample(sample, &lightSampleOffsets[i], j)
                bsdfSample = CreateBSDFSample(sample, &bsdfSampleOffsets[i], j)
            } else {
                lightSample = CreateLightSampleRandom(rng)
                bsdfSample = CreateRandomBSDFSample(rng)
            }
            Ld = Ld.Add(EstimateDirect(scene, renderer, arena, light, p, n, wo,
                rayEpsilon, time, bsdf, rng, lightSample, bsdfSample,
                BxDFType(BSDF_ALL ^ BSDF_SPECULAR)))
        }
        L = L.Add(Ld.InvScale(1.0/ float64(nSamples)))
    }
    return L
}
func UniformSampleOneLight(scene *Scene, renderer *Renderer,
    arena *MemoryArena, p *Point, n *Normal, wo *Vector,
    rayEpsilon, time float64, bsdf *BSDF,
    sample *Sample, rng *RNG, lightNumOffset int,
    lightSampleOffsets *LightSampleOffsets,
    bsdfSampleOffsets *BSDFSampleOffsets) *Spectrum {
        
    return CreateSpectrum1(0.0)	
    /*    
    // Randomly choose a single light to sample, _light_
    nLights := len(scene.lights)
    if nLights == 0 { return CreateSpectrum1(0.0) }
    var lightNum int
    if lightNumOffset != -1 {
        lightNum = Floor2Int(sample.oneD[lightNumOffset][0] * nLights)
    } else {
        lightNum = Floor2Int(rng.RandomFloat() * float64(nLights))
    }
    lightNum = Min(lightNum, nLights-1)
    light := scene.lights[lightNum]

    // Initialize light and bsdf samples for single light sample
    var lightSample *LightSample
    var bsdfSample *BSDFSample
    if lightSampleOffsets != nil && bsdfSampleOffsets != nil {
        lightSample = CreateLightSample(sample, lightSampleOffsets, 0)
        bsdfSample = CreateBSDFSample(sample, bsdfSampleOffsets, 0)
    } else {
        lightSample = CreateLightSampleRandom(rng)
        bsdfSample = CreateRandomBSDFSample(rng)
    }
    return EstimateDirect(scene, renderer, arena, light, p, n, wo,
                       rayEpsilon, time, bsdf, rng, lightSample,
                       bsdfSample, BxDFType(BSDF_ALL ^ BSDF_SPECULAR)).Scale(float64(nLights))
    */                   
}
    
func EstimateDirect(scene *Scene, renderer *Renderer,
    arena *MemoryArena, light Light, p *Point,
    n *Normal, wo *Vector, rayEpsilon, time float64, bsdf *BSDF,
    rng *RNG, lightSample *LightSample, bsdfSample *BSDFSample,
    flags BxDFType) *Spectrum {    	
    Ld := CreateSpectrum1(0.0)
/*    
    // Sample light source with multiple importance sampling
    Li, wi, lightPdf, visibility := light.Sample_L(p, rayEpsilon, lightSample, time)
    if lightPdf > 0.0 && !Li.IsBlack() {
        f := bsdf.f(wo, wi, flags)
        if !f.IsBlack() && visibility.Unoccluded(scene) {
            // Add light's contribution to reflected radiance
            Li = Li.Mult(visibility.Transmittance(scene, renderer, nil, rng, arena))
            if light.IsDeltaLight() {
                Ld += f * Li * (AbsDotVectorNormal(wi, n) / lightPdf)
            } else {
                bsdfPdf := bsdf.Pdf(wo, wi, flags)
                weight := PowerHeuristic(1, lightPdf, 1, bsdfPdf)
                Ld += f * Li * (AbsDotVectorNormal(wi, n) * weight / lightPdf)
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if !light.IsDeltaLight() {
        var sampledType BxDFType
        f := bsdf.Sample_f(wo, &wi, bsdfSample, &bsdfPdf, flags, &sampledType)
        if !f.IsBlack() && bsdfPdf > 0.0 {
            weight := 1.0
            if !(sampledType & BSDF_SPECULAR) {
                lightPdf = light.Pdf(p, wi)
                if lightPdf == 0.0 {
                    return Ld
                }
                weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf)
            }
            // Add light contribution from BSDF sampling
            var lightIsect Intersection
            Li := CreateSpectrum1(0.0)
            ray := CreateRayDifferential(p, wi, rayEpsilon, INFINITY, time)
            if scene.Intersect(ray, &lightIsect) {
                if lightIsect.primitive.GetAreaLight() == light {
                    Li = lightIsect.Le(-wi)
                }
            } else {
                Li = light.Le(ray)
            }    
            if !Li.IsBlack() {
                Li *= renderer.Transmittance(scene, ray, NULL, rng, arena)
                Ld += f * Li * AbsDot(wi, n) * weight / bsdfPdf
            }
        }
    }
*/    
    return Ld
}

func SpecularReflect(ray *RayDifferential, bsdf *BSDF, rng *RNG,
    isect *Intersection, renderer *Renderer, scene *Scene,
    sample *Sample, arena *MemoryArena) *Spectrum {
    L := CreateSpectrum1(0.0)
/*    	
    wo := -ray.d
    var wi Vector
    var pdf float64
    p := bsdf.dgShading.p
    n := bsdf.dgShading.nn
    f := bsdf.Sample_f(wo, &wi, BSDFSample(rng), &pdf, BxDFType(BSDF_REFLECTION | BSDF_SPECULAR))
    if pdf > 0.0 && !f.IsBlack() && AbsDot(wi, n) != 0.0 {
        // Compute ray differential _rd_ for specular reflection
        rd := RayDifferential(p, wi, ray, isect.rayEpsilon)
        if ray.hasDifferentials {
            rd.hasDifferentials = true
            rd.rxOrigin = p + isect.dg.dpdx
            rd.ryOrigin = p + isect.dg.dpdy
            // Compute differential reflected directions
            dndx := bsdf.dgShading.dndu * bsdf.dgShading.dudx +
                          bsdf.dgShading.dndv * bsdf.dgShading.dvdx
            dndy := bsdf.dgShading.dndu * bsdf.dgShading.dudy +
                          bsdf.dgShading.dndv * bsdf.dgShading.dvdy
            dwodx, dwody := -ray.rxDirection - wo, -ray.ryDirection - wo
            dDNdx := Dot(dwodx, n) + Dot(wo, dndx)
            dDNdy := Dot(dwody, n) + Dot(wo, dndy)
            rd.rxDirection = wi - dwodx + 2 * Vector(Dot(wo, n) * dndx + dDNdx * n)
            rd.ryDirection = wi - dwody + 2 * Vector(Dot(wo, n) * dndy + dDNdy * n)
        }
        //PBRT_STARTED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd))
        Li := renderer.Li(scene, rd, sample, rng, arena)
        L = f * Li * AbsDot(wi, n) / pdf
        //PBRT_FINISHED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd))
    }
*/    
    return L
}
    
func SpecularTransmit(ray *RayDifferential, bsdf* BSDF, rng *RNG,
    isect *Intersection, renderer *Renderer, scene *Scene,
    sample *Sample, arena *MemoryArena) *Spectrum {
    	
    L := CreateSpectrum1(0.0)
/*    
    wo := -ray.d
    var Vector wi
    var pdf float64
    p := bsdf.dgShading.p
    n := bsdf.dgShading.nn
    f := bsdf.Sample_f(wo, &wi, BSDFSample(rng), &pdf, BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR))
    if pdf > 0.0 && !f.IsBlack() && AbsDot(wi, n) != 0.0 {
        // Compute ray differential _rd_ for specular transmission
        rd := RayDifferential(p, wi, ray, isect.rayEpsilon)
        if ray.hasDifferentials {
            rd.hasDifferentials = true
            rd.rxOrigin = p + isect.dg.dpdx
            rd.ryOrigin = p + isect.dg.dpdy
        
            eta := bsdf.eta
            w := -wo
            if Dot(wo, n) < 0.0 { eta = 1.0 / eta }
        
            dndx := bsdf.dgShading.dndu * bsdf.dgShading.dudx + bsdf.dgShading.dndv * bsdf.dgShading.dvdx
            dndy := bsdf.dgShading.dndu * bsdf.dgShading.dudy + bsdf.dgShading.dndv * bsdf.dgShading.dvdy
        
            dwodx, dwody := -ray.rxDirection - wo, -ray.ryDirection - wo
            dDNdx := Dot(dwodx, n) + Dot(wo, dndx)
            dDNdy := Dot(dwody, n) + Dot(wo, dndy)
        
            mu := eta * Dot(w, n) - Dot(wi, n)
            dmudx := (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx
            dmudy := (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy
        
            rd.rxDirection = wi + eta * dwodx - Vector(mu * dndx + dmudx * n)
            rd.ryDirection = wi + eta * dwody - Vector(mu * dndy + dmudy * n)
        }
        //PBRT_STARTED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
        Li := renderer.Li(scene, rd, sample, rng, arena)
        L = f * Li * AbsDot(wi, n) / pdf
        //PBRT_FINISHED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
    }
*/    
    return L
}

func ComputeLightSamplingCDF(scene *Scene) *Distribution1D {
    nLights := len(scene.lights)
    lightPower := make([]float64, nLights, nLights)
    for i, light := range scene.lights {
        lightPower[i] = light.Power(scene).Y()
    }
    return CreateDistribution1D(lightPower)
}
