package core

import (
	"math"
	"strings"
	"os"
	"fmt"
)

type LightStrategy int

const (
	SAMPLE_ALL_UNIFORM = iota
	SAMPLE_ONE_UNIFORM
)
const (
	SAMPLE_DEPTH = 3
)

type (
	Integrator interface {
		Preprocess(scene *Scene, camera Camera, renderer Renderer)
		RequestSamples(sampler Sampler, sample *Sample, scene *Scene)
	}

	SurfaceIntegrator interface {
		Integrator
		Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
			sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum
	}

	AmbientOcclusionIntegrator struct {
		nSamples int
		maxDist  float64
	}

	DiffusePRTIntegrator struct {
		lmax, nSamples int
		c_in           []Spectrum
	}

	DirectLightingIntegrator struct {
		strategy LightStrategy
		maxDepth int

		// Declare sample parameters for light source sampling
		lightSampleOffsets []LightSampleOffsets
		bsdfSampleOffsets  []BSDFSampleOffsets
		lightNumOffset     int
	}

	GlossyPRTIntegrator struct {
		Kd, Ks         Spectrum
		roughness      float64
		lmax, nSamples int
		c_in           []Spectrum
		B              []Spectrum
	}

	IrradianceCacheIntegrator struct {
		minSamplePixelSpacing, maxSamplePixelSpacing float64
		minWeight, cosMaxSampleAngleDifference       float64
		nSamples, maxSpecularDepth, maxIndirectDepth int

		// Declare sample parameters for light source sampling
		lightSampleOffsets []LightSampleOffsets
		bsdfSampleOffsets  []BSDFSampleOffsets
		//octree []*OctreeIrradianceSample
	}

	PathIntegrator struct {
		maxDepth           int
		lightSampleOffsets [SAMPLE_DEPTH]LightSampleOffsets
		lightNumOffset     [SAMPLE_DEPTH]int
		bsdfSampleOffsets  [SAMPLE_DEPTH]BSDFSampleOffsets
		pathSampleOffsets  [SAMPLE_DEPTH]BSDFSampleOffsets
	}

	PhotonIntegrator struct {
		nCausticPhotonsWanted, nIndirectPhotonsWanted, nLookup int
		maxDistSquared                                         float64
		maxSpecularDepth, maxPhotonDepth                       int
		finalGather                                            bool
		gatherSamples                                          int
		cosGatherAngle                                         float64

		// Declare sample parameters for light source sampling
		lightSampleOffsets                                []LightSampleOffsets
		bsdfSampleOffsets                                 []BSDFSampleOffsets
		bsdfGatherSampleOffsets, indirGatherSampleOffsets BSDFSampleOffsets
		nCausticPaths, nIndirectPaths                     int
		//causticMap, indirectMap *KdTreePhoton
		//radianceMap *KdTreeRadiancePhoton
	}

	UseRadianceProbes struct {
		bbox                                                 *BBox
		lmax, includeDirectInProbes, includeIndirectInProbes int
		nProbes                                              [3]int
		c_in                                                 []Spectrum

		// Declare sample parameters for light source sampling
		lightSampleOffsets []LightSampleOffsets
		bsdfSampleOffsets  []BSDFSampleOffsets
	}

	WhittedIntegrator struct {
		maxDepth int
	}
)

func (*AmbientOcclusionIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (*AmbientOcclusionIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *AmbientOcclusionIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    bsdf := isect.GetBSDF(ray, arena)
    p := bsdf.dgShading.p
    n := FaceforwardNormalVector(isect.dg.nn, ray.dir.Negate())

    scramble := [2]uint32{ rng.RandomUInt(), rng.RandomUInt() }
    u := [2]float64{0.0, 0.0}
    nClear := 0
    var i uint32
    for i = 0; i < uint32(integrator.nSamples); i++ {
        Sample02(i, scramble, u[:])
        w := UniformSampleSphere(u[0], u[1])
        if DotVectorNormal(w, n) < 0.0 { w = w.Negate() }
        r := CreateRay(p, w, 0.01, integrator.maxDist, 0.0, 0)
        if !scene.IntersectP(r) { nClear++ }
    }
    return NewSpectrum1(float64(nClear) / float64(integrator.nSamples))
}

func (integrator *DiffusePRTIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {
    bbox := scene.WorldBound()
    p := bbox.pMin.Add(CreateVectorFromPoint(&bbox.pMax)).Scale(0.5)
    rng := NewRNG(1)
    SHProjectIncidentDirectRadiance(p, 0.0, camera.ShutterOpen(), nil, scene, false, integrator.lmax, rng, integrator.c_in)
}

func (integrator *DiffusePRTIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *DiffusePRTIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    L := NewSpectrum1(0.0)
    wo := ray.dir.Negate()
    // Compute emitted light if ray hit an area light source
    L = L.Add(isect.Le(wo))

    // Evaluate BSDF at hit point
    bsdf := isect.GetBSDF(ray, arena)
    p := bsdf.dgShading.p
    n := bsdf.dgShading.nn
    // Compute reflected radiance using diffuse PRT

    // Project diffuse transfer function at point to SH
    c_transfer := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHComputeDiffuseTransfer(p, FaceforwardNormalVector(n, wo), isect.rayEpsilon,
                             scene, rng, integrator.nSamples, integrator.lmax, c_transfer)

    // Compute integral of product of incident radiance and transfer function
    Kd := bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION, 6).InvScale(math.Pi)
    Lo := NewSpectrum1(0.0)
    for i := 0; i < SHTerms(integrator.lmax); i++ {
        Lo = Lo.Add(integrator.c_in[i].Mult(&c_transfer[i]))
    }    
    return L.Add(Kd.Mult(Lo.Clamp(0.0, INFINITY)))
}

func (*DirectLightingIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer) {}
func (integrator *DirectLightingIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
	if integrator.strategy == SAMPLE_ALL_UNIFORM {
		// Allocate and request samples for sampling all lights
		nLights := len(scene.lights)
		integrator.lightSampleOffsets = make([]LightSampleOffsets, nLights, nLights)
		integrator.bsdfSampleOffsets = make([]BSDFSampleOffsets, nLights, nLights)
		for i := 0; i < nLights; i++ {
			light := scene.lights[i]
			nSamples := light.NumSamples()
			if sampler != nil {
				nSamples = sampler.RoundSize(nSamples)
			}
			integrator.lightSampleOffsets[i] = *CreateLightSampleOffsets(nSamples, sample)
			integrator.bsdfSampleOffsets[i] = *CreateBSDFSampleOffsets(nSamples, sample)
		}
		integrator.lightNumOffset = -1
	} else {
		// Allocate and request samples for sampling one light
		integrator.lightSampleOffsets = make([]LightSampleOffsets, 1, 1)
		integrator.lightSampleOffsets[0] = *CreateLightSampleOffsets(1, sample)
		integrator.lightNumOffset = sample.Add1D(1)
		integrator.bsdfSampleOffsets = make([]BSDFSampleOffsets, 1, 1)
		integrator.bsdfSampleOffsets[0] = *CreateBSDFSampleOffsets(1, sample)
	}
}

func (integrator *DirectLightingIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	// Evaluate BSDF at hit point
	bsdf := isect.GetBSDF(ray, arena)
	Assert(bsdf != nil)
	wo := ray.dir.Negate()
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	// Compute emitted light if ray hit an area light source
	L = L.Add(isect.Le(wo))

	// Compute direct lighting for _DirectLightingIntegrator_ integrator
	if len(scene.lights) > 0 {
		// Apply direct lighting strategy
		switch integrator.strategy {
		case SAMPLE_ALL_UNIFORM:
			L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n, wo,
				isect.rayEpsilon, ray.time, bsdf, sample, rng,
				integrator.lightSampleOffsets, integrator.bsdfSampleOffsets))
		case SAMPLE_ONE_UNIFORM:
			L = L.Add(UniformSampleOneLight(scene, renderer, arena, p, n, wo,
				isect.rayEpsilon, ray.time, bsdf, sample, rng,
				integrator.lightNumOffset, &integrator.lightSampleOffsets[0], &integrator.bsdfSampleOffsets[0]))
		}
	}
	if ray.depth+1 < integrator.maxDepth {
		// Trace rays for specular reflection and refraction
		L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
		L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
	}
	return L
}

func NewGlossyPRTIntegrator(Kd, Ks Spectrum, roughness float64, lmax, ns int) *GlossyPRTIntegrator {
	return &GlossyPRTIntegrator{Kd, Ks, roughness, lmax, ns, nil, nil}
}

func (integrator *GlossyPRTIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {
    // Project direct lighting into SH for _GlossyPRTIntegrator_
    bbox := scene.WorldBound()
    p := bbox.pMin.Add(CreateVectorFromPoint(&bbox.pMax)).Scale(0.5)
    rng := NewRNG(42)
    integrator.c_in = make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHProjectIncidentDirectRadiance(p, 0.0, camera.ShutterOpen(), nil,
        scene, false, integrator.lmax, rng, integrator.c_in)

    // Compute glossy BSDF matrix for PRT
    integrator.B = make([]Spectrum, SHTerms(integrator.lmax)*SHTerms(integrator.lmax), SHTerms(integrator.lmax)*SHTerms(integrator.lmax))
    SHComputeBSDFMatrix(&integrator.Kd, &integrator.Ks, integrator.roughness, rng, 1024, integrator.lmax, integrator.B)
}

func (*GlossyPRTIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}

func (integrator *GlossyPRTIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    L := NewSpectrum1(0.0)
    wo := ray.dir.Negate()
    // Compute emitted light if ray hit an area light source
    L = L.Add(isect.Le(wo))

    // Evaluate BSDF at hit point
    bsdf := isect.GetBSDF(ray, arena)
    p := bsdf.dgShading.p
    // Compute reflected radiance with glossy PRT at point

    // Compute SH radiance transfer matrix at point and SH coefficients
    c_t := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    T := make([]Spectrum, SHTerms(integrator.lmax)*SHTerms(integrator.lmax), SHTerms(integrator.lmax)*SHTerms(integrator.lmax))
    SHComputeTransferMatrix(p, isect.rayEpsilon, scene, rng, integrator.nSamples, integrator.lmax, T)
    SHMatrixVectorMultiply(T, integrator.c_in, c_t, integrator.lmax)

    // Rotate incident SH lighting to local coordinate frame
    r1 := bsdf.LocalToWorld(CreateVector(1,0,0))
    r2 := bsdf.LocalToWorld(CreateVector(0,1,0))
    nl := CreateNormalFromVector(bsdf.LocalToWorld(CreateVector(0,0,1)))
    rot := NewMatrix4x4(r1.x, r2.x, nl.x, 0,
                  r1.y, r2.y, nl.y, 0,
                  r1.z, r2.z, nl.z, 0,
                     0,    0,    0, 1)
    c_l := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHRotate(c_t, c_l, rot, integrator.lmax, arena)

    // Compute final coefficients _c\_o_ using BSDF matrix
    c_o := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHMatrixVectorMultiply(integrator.B, c_l, c_o, integrator.lmax)

    // Evaluate outgoing radiance function for $\wo$ and add to _L_
    woLocal := bsdf.WorldToLocal(wo)
    Ylm := make([]float64, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHEvaluate(woLocal, integrator.lmax, Ylm)
    Li := NewSpectrum1(0.0)
    for i := 0; i < SHTerms(integrator.lmax); i++ {
        Li = Li.Add(c_o[i].Scale(Ylm[i]))
    }    
    L = L.Add(Li.Clamp(0.0, INFINITY))

    return L
}

func NewIrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing, maxAngle float64,
        maxSpecularDepth, maxIndirectDepth, nSamples int) *IrradianceCacheIntegrator {
	Unimplemented()
	return nil
}
func (integrator *IrradianceCacheIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (integrator *IrradianceCacheIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *IrradianceCacheIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	Unimplemented()
	return nil
}

func NewPathIntegrator(maxDepth int) *PathIntegrator {
	integrator := new(PathIntegrator)
	integrator.maxDepth = maxDepth
	return integrator
}
func (*PathIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (integrator *PathIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
    for i := 0; i < SAMPLE_DEPTH; i++ {
        integrator.lightSampleOffsets[i] = *CreateLightSampleOffsets(1, sample)
        integrator.lightNumOffset[i] = sample.Add1D(1)
        integrator.bsdfSampleOffsets[i] = *CreateBSDFSampleOffsets(1, sample)
        integrator.pathSampleOffsets[i] = *CreateBSDFSampleOffsets(1, sample)
    }
}
func (integrator *PathIntegrator) Li(scene *Scene, renderer Renderer, r *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    // Declare common path integration variables
    pathThroughput, L := NewSpectrum1(1.0), NewSpectrum1(0.0)
    ray := &RayDifferential{Ray{r.origin, r.dir, r.mint, r.maxt, r.time, r.depth}, r.hasDifferentials, r.rxOrigin, r.ryOrigin, r.rxDirection, r.ryDirection}
    specularBounce := false
    localIsect := new(Intersection)
    isectp := isect
    for bounces := 0; ; bounces++ {
        // Possibly add emitted light at path vertex
        if bounces == 0 || specularBounce {
            L = L.Add(pathThroughput.Mult(isectp.Le(ray.dir.Negate())))
		}
        // Sample illumination from lights to find path contribution
        bsdf := isectp.GetBSDF(ray, arena)
        p := bsdf.dgShading.p
        n := bsdf.dgShading.nn
        wo := ray.dir.Negate()
        if bounces < SAMPLE_DEPTH {
            L = L.Add(pathThroughput.Mult(UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp.rayEpsilon, ray.time, bsdf, sample, rng,
                     integrator.lightNumOffset[bounces], &integrator.lightSampleOffsets[bounces],
                     &integrator.bsdfSampleOffsets[bounces])))
        } else {
            L = L.Add(pathThroughput.Mult(UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp.rayEpsilon, ray.time, bsdf, sample, rng, 0, nil, nil)))
		}
        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        var outgoingBSDFSample *BSDFSample
        if bounces < SAMPLE_DEPTH {
            outgoingBSDFSample = CreateBSDFSample(sample, &integrator.pathSampleOffsets[bounces], 0)
        } else {
            outgoingBSDFSample = CreateRandomBSDFSample(rng)
        } 
           
        f, wi, pdf, flags := bsdf.Sample_f(wo, outgoingBSDFSample, BSDF_ALL)
        if f.IsBlack() || pdf == 0.0 {
            break
        }    
        specularBounce = (flags & BSDF_SPECULAR) != 0
        pathThroughput = pathThroughput.Mult(f.Scale(AbsDotVectorNormal(wi, n) / pdf))
        ray = CreateChildRayDifferential(p, wi, CreateRayFromRayDifferential(ray), isectp.rayEpsilon, INFINITY)

        // Possibly terminate the path
        if bounces > 3 {
            continueProbability := math.Min(0.5, pathThroughput.Y())
            if rng.RandomFloat() > continueProbability {
                break
            }    
            pathThroughput = pathThroughput.InvScale(continueProbability)
        }
        if bounces == integrator.maxDepth {
            break
		}
        // Find next vertex of path
        var ok bool
        if ok, localIsect = scene.Intersect(ray); ok {
            if specularBounce {
                for i := 0; i < len(scene.lights); i++ {
                   L = L.Add(pathThroughput.Mult(scene.lights[i].Le(ray)))
                }
            }       
            break
        }
        pathThroughput = pathThroughput.Mult(renderer.Transmittance(scene, ray, nil, rng, arena))
        isectp = localIsect
    }
    return L
}

func NewPhotonIntegrator(nCaustic, nIndirect, nUsed, maxSpecularDepth, maxPhotonDepth int, maxDist float64, finalGather bool, gatherSamples int, gatherAngle float64) *PhotonIntegrator {
	Unimplemented()
	return nil
}
func (integrator *PhotonIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (integrator *PhotonIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *PhotonIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	Unimplemented()
	return nil
}

func NewUseRadianceProbes(filename string) *UseRadianceProbes {
	probes := new(UseRadianceProbes)
    probes.lightSampleOffsets = nil
    probes.bsdfSampleOffsets = nil
    // Read precomputed radiance probe values from file
    f, err := os.Open(filename)
    defer f.Close()
    if err == nil {
		_, err = fmt.Fscanf(f, "%d %d %d", &probes.lmax, &probes.includeDirectInProbes,&probes.includeIndirectInProbes)
		if err != nil {
            Severe("Error reading data from radiance probe file \"%s\" Error: %s", filename, err)
		}	
		_, err = fmt.Fscanf(f, "%d %d %d", &probes.nProbes[0], &probes.nProbes[1], &probes.nProbes[2])
		if err != nil {
            Severe("Error reading data from radiance probe file \"%s\" Error: %s", filename, err)			
		}
        _, err = fmt.Fscanf(f, "%f %f %f %f %f %f", &probes.bbox.pMin.x, &probes.bbox.pMin.y, &probes.bbox.pMin.z,
                   &probes.bbox.pMax.x, &probes.bbox.pMax.y, &probes.bbox.pMax.z)
		if err != nil {
            Severe("Error reading data from radiance probe file \"%s\" Error: %s", filename, err)			
		}
    
        probes.c_in = make([]Spectrum, SHTerms(probes.lmax) * probes.nProbes[0] * probes.nProbes[1] * probes.nProbes[2], SHTerms(probes.lmax) * probes.nProbes[0] * probes.nProbes[1] * probes.nProbes[2])
        offset := 0
        for i := 0; i < probes.nProbes[0] * probes.nProbes[1] * probes.nProbes[2]; i++ {
            for j := 0; j < SHTerms(probes.lmax); j++ {
                if !probes.c_in[offset].Read(f) {
                    Severe("Error reading data from radiance probe file \"%s\"", filename)
                }
                offset++
            }     
        }
    } else {
        Error("Unable to read saved radiance volume values from file \"%s\"", filename)
    }
    return probes    
}

func (*UseRadianceProbes) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}

func (integrator *UseRadianceProbes) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
    // Allocate and request samples for sampling all lights
    nLights := len(scene.lights)
    integrator.lightSampleOffsets = make([]LightSampleOffsets, nLights, nLights)
    integrator.bsdfSampleOffsets = make([]BSDFSampleOffsets, nLights, nLights)
    for i := 0; i < nLights; i++ {
        light := scene.lights[i]
        nSamples := light.NumSamples()
        if sampler != nil { nSamples = sampler.RoundSize(nSamples) }
        integrator.lightSampleOffsets[i] = *CreateLightSampleOffsets(nSamples, sample)
        integrator.bsdfSampleOffsets[i] = *CreateBSDFSampleOffsets(nSamples, sample)
    }
}

func (integrator *UseRadianceProbes) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    L := NewSpectrum1(0.0)
    wo := ray.dir.Negate()
    // Compute emitted light if ray hit an area light source
    L = L.Add(isect.Le(wo))

    // Evaluate BSDF at hit point
    bsdf := isect.GetBSDF(ray, arena)
    p := bsdf.dgShading.p
    n := bsdf.dgShading.nn
    // Compute reflection for radiance probes integrator
    if integrator.includeDirectInProbes != 0 {
        L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n,
                wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
                integrator.lightSampleOffsets, integrator.bsdfSampleOffsets))
	}
    // Compute reflected lighting using radiance probes

    // Compute probe coordinates and offsets for lookup point
    offset := integrator.bbox.Offset(p)
     voxx := (offset.x * float64(integrator.nProbes[0])) - 0.5
     voxy := (offset.y * float64(integrator.nProbes[1])) - 0.5
     voxz := (offset.z * float64(integrator.nProbes[2])) - 0.5
     vx, vy, vz := Floor2Int(voxx), Floor2Int(voxy), Floor2Int(voxz)
     dx, dy, dz := voxx - float64(vx), voxy - float64(vy), voxz - float64(vz)

    // Get radiance probe coefficients around lookup point
    b000 := integrator.c_inXYZ(integrator.lmax, vx,   vy,   vz)
    b100 := integrator.c_inXYZ(integrator.lmax, vx+1, vy,   vz)
    b010 := integrator.c_inXYZ(integrator.lmax, vx,   vy+1, vz)
    b110 := integrator.c_inXYZ(integrator.lmax, vx+1, vy+1, vz)
    b001 := integrator.c_inXYZ(integrator.lmax, vx,   vy,   vz+1)
    b101 := integrator.c_inXYZ(integrator.lmax, vx+1, vy,   vz+1)
    b011 := integrator.c_inXYZ(integrator.lmax, vx,   vy+1, vz+1)
    b111 := integrator.c_inXYZ(integrator.lmax, vx+1, vy+1, vz+1)

    // Compute incident radiance from radiance probe coefficients
    c_inp := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    for i := 0; i < SHTerms(integrator.lmax); i++ {
        // Do trilinear interpolation to compute SH coefficients at point
        c00 := LerpSpectrum(dx, &b000[i], &b100[i])
        c10 := LerpSpectrum(dx, &b010[i], &b110[i])
        c01 := LerpSpectrum(dx, &b001[i], &b101[i])
        c11 := LerpSpectrum(dx, &b011[i], &b111[i])
        c0 := LerpSpectrum(dy, c00, c10)
        c1 := LerpSpectrum(dy, c01, c11)
        c_inp[i] = *LerpSpectrum(dz, c0, c1)
    }

    // Convolve incident radiance to compute irradiance function
    c_E := make([]Spectrum, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHConvolveCosTheta(integrator.lmax, c_inp, c_E)

    // Evaluate irradiance function and accumulate reflection
    rho := bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION, 6)
    Ylm := make([]float64, SHTerms(integrator.lmax), SHTerms(integrator.lmax))
    SHEvaluate(CreateVectorFromNormal(FaceforwardNormalVector(n, wo)), integrator.lmax, Ylm)
    E := NewSpectrum1(0.0)
    for i := 0; i < SHTerms(integrator.lmax); i++ {
        E = E.Add(c_E[i].Scale(Ylm[i]))
    }    
    L = L.Add(rho.Mult(E.Clamp(0.0, INFINITY).Scale(INV_PI)))
    
    return L
}

func (integrator *UseRadianceProbes) c_inXYZ(lmax, vx, vy, vz int) []Spectrum {
	vx = Clampi(vx, 0, integrator.nProbes[0]-1)
	vy = Clampi(vy, 0, integrator.nProbes[1]-1)
	vz = Clampi(vz, 0, integrator.nProbes[2]-1)
	offset := vx + vy * integrator.nProbes[0] + vz * integrator.nProbes[0] * integrator.nProbes[1]
	return integrator.c_in[SHTerms(lmax) * offset:]
}

func (*WhittedIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {}
func (*WhittedIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {}
func (integrator *WhittedIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	// Compute emitted and reflected light at ray intersection point

	// Evaluate BSDF at hit point
	bsdf := isect.GetBSDF(ray, arena)

	// Initialize common variables for Whitted integrator
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	wo := ray.dir.Negate()

	// Compute emitted light if ray hit an area light source
	L = L.Add(isect.Le(wo))

	// Add contribution of each light source
	for _, light := range scene.lights {
		var visibility VisibilityTester
		Li, wi, pdf := light.Sample_L(p, isect.rayEpsilon, CreateLightSampleRandom(rng), ray.time, &visibility)
		if Li.IsBlack() || pdf == 0.0 {
			continue
		}
		f := bsdf.f(wo, wi, BSDF_ALL)
		if !f.IsBlack() && visibility.Unoccluded(scene) {
			L = L.Add(f.Mult(Li.Scale(AbsDotVectorNormal(wi, n)).Mult(visibility.Transmittance(scene, renderer, sample, rng, arena).Scale(1.0 / pdf))))
		}
	}
	if ray.depth+1 < integrator.maxDepth {
		// Trace rays for specular reflection and refraction
		L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
		L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
	}
	return L
}

func CreateAmbientOcclusionIntegrator(params *ParamSet) *AmbientOcclusionIntegrator {
     nSamples := params.FindIntParam("nsamples", 2048)
     maxDist := params.FindFloatParam("maxdist", INFINITY)
    if options.QuickRender { nSamples = Maxi(1, nSamples / 4) }
    return &AmbientOcclusionIntegrator{nSamples, maxDist}
}

func CreateDiffusePRTIntegratorSurfaceIntegrator(params *ParamSet) *DiffusePRTIntegrator {
    lmax := params.FindIntParam("lmax", 4)
    ns := params.FindIntParam("nsamples", 4096)
    integrator := &DiffusePRTIntegrator{lmax, int(RoundUpPow2(uint32(ns))), nil}
    integrator.c_in = make([]Spectrum, SHTerms(lmax), SHTerms(lmax))
    return integrator
}

func CreateGlossyPRTIntegratorSurfaceIntegrator(params *ParamSet) *GlossyPRTIntegrator {
     lmax := params.FindIntParam("lmax", 4)
     ns := params.FindIntParam("nsamples", 4096)
     Kd := params.FindSpectrumParam("Kd", *NewSpectrum1(0.5))
     Ks := params.FindSpectrumParam("Ks", *NewSpectrum1(0.25))
     roughness := params.FindFloatParam("roughness", 0.1)
    return NewGlossyPRTIntegrator(Kd, Ks, roughness, lmax, ns)
}
	
func CreateIrradianceCacheIntegrator(params *ParamSet) *IrradianceCacheIntegrator {
     minWeight := params.FindFloatParam("minweight", 0.5)
     minSpacing := params.FindFloatParam("minpixelspacing", 2.5)
     maxSpacing := params.FindFloatParam("maxpixelspacing", 15.0)
     maxAngle := params.FindFloatParam("maxangledifference", 10.0)
     maxSpecularDepth := params.FindIntParam("maxspeculardepth", 5)
     maxIndirectDepth := params.FindIntParam("maxindirectdepth", 3)
     nSamples := params.FindIntParam("nsamples", 4096)
    if options.QuickRender { nSamples = Maxi(1, nSamples / 16) }
    return NewIrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing, maxAngle,
        maxSpecularDepth, maxIndirectDepth, nSamples)
}

func CreatePathSurfaceIntegrator(params *ParamSet) *PathIntegrator        { 
	maxDepth := params.FindIntParam("maxdepth", 5)
    return NewPathIntegrator(maxDepth)
}
	
func CreatePhotonMapSurfaceIntegrator(params *ParamSet) *PhotonIntegrator { 
    nCaustic := params.FindIntParam("causticphotons", 20000)
    nIndirect := params.FindIntParam("indirectphotons", 100000)
    nUsed := params.FindIntParam("nused", 50)
    if options.QuickRender { 
		nCaustic = nCaustic / 10
		nIndirect = nIndirect / 10
		nUsed = Maxi(1, nUsed / 10) 
    }
    maxSpecularDepth := params.FindIntParam("maxspeculardepth", 5)
    maxPhotonDepth := params.FindIntParam("maxphotondepth", 5)
    finalGather := params.FindBoolParam("finalgather", true)
    gatherSamples := params.FindIntParam("finalgathersamples", 32)
    if options.QuickRender { gatherSamples = Maxi(1, gatherSamples / 4) }
    maxDist := params.FindFloatParam("maxdist", 0.1)
    gatherAngle := params.FindFloatParam("gatherangle", 10.0)
    return NewPhotonIntegrator(nCaustic, nIndirect, nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples, gatherAngle)
}
	
func CreateRadianceProbesSurfaceIntegrator(params *ParamSet) *UseRadianceProbes {
    filename := params.FindFilenameParam("filename", "probes.out")
    return NewUseRadianceProbes(filename)
}

func CreateDirectLightingIntegrator(params *ParamSet) *DirectLightingIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	var strategy LightStrategy
	st := params.FindStringParam("strategy", "all")
	if strings.Compare(st, "one") == 0 {
		strategy = SAMPLE_ONE_UNIFORM
	} else if strings.Compare(st, "all") == 0 {
		strategy = SAMPLE_ALL_UNIFORM
	} else {
		Warning("Strategy \"%s\" for direct lighting unknown.  Using \"all\".", st)
		strategy = SAMPLE_ALL_UNIFORM
	}
	return &DirectLightingIntegrator{strategy, maxDepth, nil, nil, 0}
}

func CreateWhittedSurfaceIntegrator(params *ParamSet) *WhittedIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	return &WhittedIntegrator{maxDepth}
}

func UniformSampleAllLights(scene *Scene, renderer Renderer,
	arena *MemoryArena, p *Point, n *Normal, wo *Vector,
	rayEpsilon, time float64, bsdf *BSDF, sample *Sample, rng *RNG,
	lightSampleOffsets []LightSampleOffsets,
	bsdfSampleOffsets []BSDFSampleOffsets) *Spectrum {

	L := NewSpectrum1(0.0)
	for i, light := range scene.lights {
		nSamples := 1
		if lightSampleOffsets != nil {
			nSamples = lightSampleOffsets[i].nSamples
		}
		// Estimate direct lighting from _light_ samples
		Ld := NewSpectrum1(0.0)
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
				BxDFType(BSDF_ALL^BSDF_SPECULAR)))
		}
		L = L.Add(Ld.InvScale(1.0 / float64(nSamples)))
	}
	return L
}

func UniformSampleOneLight(scene *Scene, renderer Renderer,
	arena *MemoryArena, p *Point, n *Normal, wo *Vector,
	rayEpsilon, time float64, bsdf *BSDF,
	sample *Sample, rng *RNG, lightNumOffset int,
	lightSampleOffsets *LightSampleOffsets,
	bsdfSampleOffsets *BSDFSampleOffsets) *Spectrum {

	// Randomly choose a single light to sample, _light_
	nLights := len(scene.lights)
	if nLights == 0 {
		return NewSpectrum1(0.0)
	}
	var lightNum int
	if lightNumOffset != -1 {
		lightNum = Floor2Int(sample.oneD[lightNumOffset][0] * float64(nLights))
	} else {
		lightNum = Floor2Int(rng.RandomFloat() * float64(nLights))
	}
	lightNum = Mini(lightNum, nLights-1)
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
		bsdfSample, BxDFType(BSDF_ALL^BSDF_SPECULAR)).Scale(float64(nLights))
}

func EstimateDirect(scene *Scene, renderer Renderer,
	arena *MemoryArena, light Light, p *Point,
	n *Normal, wo *Vector, rayEpsilon, time float64, bsdf *BSDF,
	rng *RNG, lightSample *LightSample, bsdfSample *BSDFSample,
	flags BxDFType) *Spectrum {
	Ld := NewSpectrum1(0.0)

	// Sample light source with multiple importance sampling
	var visibility VisibilityTester
	Li, wi, lightPdf := light.Sample_L(p, rayEpsilon, lightSample, time, &visibility)
	if lightPdf > 0.0 && !Li.IsBlack() {
		f := bsdf.f(wo, wi, flags)
		if !f.IsBlack() && visibility.Unoccluded(scene) {
			// Add light's contribution to reflected radiance
			Li = Li.Mult(visibility.Transmittance(scene, renderer, nil, rng, arena))
			if light.IsDeltaLight() {
				Ld = Ld.Add(f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) / lightPdf)))
			} else {
				bsdfPdf := bsdf.Pdf(wo, wi, flags)
				weight := PowerHeuristic(1, lightPdf, 1, bsdfPdf)
				Ld = Ld.Add(f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) * weight / lightPdf)))
			}
		}
	}

	// Sample BSDF with multiple importance sampling
	if !light.IsDeltaLight() {
		f, wi, bsdfPdf, sampledType := bsdf.Sample_f(wo, bsdfSample, flags)
		if !f.IsBlack() && bsdfPdf > 0.0 {
			weight := 1.0
			if sampledType&BSDF_SPECULAR != 0 {
				lightPdf = light.Pdf(p, wi)
				if lightPdf == 0.0 {
					return Ld
				}
				weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf)
			}
			// Add light contribution from BSDF sampling
			Li := NewSpectrum1(0.0)
			ray := CreateRayDifferential(p, wi, rayEpsilon, INFINITY, time, 0)
			if hit, lightIsect := scene.Intersect(ray); hit {
				if lightIsect.primitive.GetAreaLight() == light {
					Li = lightIsect.Le(wi.Negate())
				}
			} else {
				Li = light.Le(ray)
			}
			if !Li.IsBlack() {
				Li = Li.Mult(renderer.Transmittance(scene, ray, nil, rng, arena))
				Ld = Ld.Add(f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) * weight / bsdfPdf)))
			}
		}
	}

	return Ld
}

func SpecularReflect(ray *RayDifferential, bsdf *BSDF, rng *RNG,
	isect *Intersection, renderer Renderer, scene *Scene,
	sample *Sample, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)

	wo := ray.dir.Negate()
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	f, wi, pdf, _ := bsdf.Sample_f(wo, CreateRandomBSDFSample(rng), BxDFType(BSDF_REFLECTION|BSDF_SPECULAR))
	if pdf > 0.0 && !f.IsBlack() && AbsDotVectorNormal(wi, n) != 0.0 {
		// Compute ray differential _rd_ for specular reflection
		rd := CreateChildRayDifferential(p, wi, CreateRayFromRayDifferential(ray), isect.rayEpsilon, INFINITY)
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = *p.Add(isect.dg.dpdx)
			rd.ryOrigin = *p.Add(isect.dg.dpdy)
			// Compute differential reflected directions
			dndx := bsdf.dgShading.dndu.Scale(bsdf.dgShading.dudx).Add(bsdf.dgShading.dndv.Scale(bsdf.dgShading.dvdx))
			dndy := bsdf.dgShading.dndu.Scale(bsdf.dgShading.dudy).Add(bsdf.dgShading.dndv.Scale(bsdf.dgShading.dvdy))
			dwodx, dwody := (ray.rxDirection.Negate()).Sub(wo), (ray.ryDirection.Negate()).Sub(wo)
			dDNdx := DotVectorNormal(dwodx, n) + DotVectorNormal(wo, dndx)
			dDNdy := DotVectorNormal(dwody, n) + DotVectorNormal(wo, dndy)
			rd.rxDirection = *wi.Sub(dwodx).Add(CreateVectorFromNormal(dndx.Scale(DotVectorNormal(wo, n)).Add(n.Scale(dDNdx))).Scale(2))
			rd.ryDirection = *wi.Sub(dwody).Add(CreateVectorFromNormal(dndy.Scale(DotVectorNormal(wo, n)).Add(n.Scale(dDNdy))).Scale(2))
		}
		//PBRT_STARTED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd))
		Li, _, _ := renderer.Li(scene, rd, sample, rng, arena)
		L = f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) / pdf))
		//PBRT_FINISHED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd))
	}

	return L
}

func SpecularTransmit(ray *RayDifferential, bsdf *BSDF, rng *RNG,
	isect *Intersection, renderer Renderer, scene *Scene,
	sample *Sample, arena *MemoryArena) *Spectrum {

	L := NewSpectrum1(0.0)

	wo := ray.dir.Negate()
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	f, wi, pdf, _ := bsdf.Sample_f(wo, CreateRandomBSDFSample(rng), BxDFType(BSDF_TRANSMISSION|BSDF_SPECULAR))
	if pdf > 0.0 && !f.IsBlack() && AbsDotVectorNormal(wi, n) != 0.0 {
		// Compute ray differential _rd_ for specular transmission
		rd := CreateChildRayDifferential(p, wi, CreateRayFromRayDifferential(ray), isect.rayEpsilon, INFINITY)
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = *p.Add(isect.dg.dpdx)
			rd.ryOrigin = *p.Add(isect.dg.dpdy)

			eta := bsdf.eta
			w := wo.Negate()
			if DotVectorNormal(wo, n) < 0.0 {
				eta = 1.0 / eta
			}

			dndx := bsdf.dgShading.dndu.Scale(bsdf.dgShading.dudx).Add(bsdf.dgShading.dndv.Scale(bsdf.dgShading.dvdx))
			dndy := bsdf.dgShading.dndu.Scale(bsdf.dgShading.dudy).Add(bsdf.dgShading.dndv.Scale(bsdf.dgShading.dvdy))
			dwodx, dwody := (ray.rxDirection.Negate()).Sub(wo), (ray.ryDirection.Negate()).Sub(wo)
			dDNdx := DotVectorNormal(dwodx, n) + DotVectorNormal(wo, dndx)
			dDNdy := DotVectorNormal(dwody, n) + DotVectorNormal(wo, dndy)

			mu := eta*DotVectorNormal(w, n) - DotVectorNormal(wi, n)
			dmudx := (eta - (eta*eta*DotVectorNormal(w, n))/DotVectorNormal(wi, n)) * dDNdx
			dmudy := (eta - (eta*eta*DotVectorNormal(w, n))/DotVectorNormal(wi, n)) * dDNdy

			rd.rxDirection = *wi.Add(dwodx.Scale(eta)).Sub(CreateVectorFromNormal(dndx.Scale(mu).Add(n.Scale(dmudx))))
			rd.ryDirection = *wi.Add(dwody.Scale(eta)).Sub(CreateVectorFromNormal(dndy.Scale(mu).Add(n.Scale(dmudy))))
		}
		//PBRT_STARTED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
		Li, _, _ := renderer.Li(scene, rd, sample, rng, arena)
		L = f.Mult(Li.Scale(AbsDotVectorNormal(wi, n) / pdf))
		//PBRT_FINISHED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
	}

	return L
}

func ComputeLightSamplingCDF(scene *Scene) *Distribution1D {
	nLights := len(scene.lights)
	lightPower := make([]float64, nLights, nLights)
	for i, light := range scene.lights {
		lightPower[i] = float64(light.Power(scene).Y())
	}
	return NewDistribution1D(lightPower)
}
