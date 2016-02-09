/*
	gopbrt

	Port of pbrt v2.0.0 by Matt Pharr and Greg Humphreys to the go language.
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

	The MIT License (MIT)
	Copyright (c) 2016 Rick Weyrauch

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
package core

import (
	"math"
)

type (

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
		causticMap, indirectMap                           *KdTree // Photons
		radianceMap                                       *KdTree // RadiancePhoton
	}
	
	photon struct {
		p Point
		alpha Spectrum
		wi Vector
	}
	
	radiancePhoton struct {
		p Point
		n Normal
		Lo Spectrum
	}
	
	closePhoton struct {
		photon *photon
		distanceSquared float64
	}
	
	photonShootingTask struct {
    	taskNum int
     	time float64
     	integrator *PhotonIntegrator
     	progress *ProgressReporter
     	abortTasks bool
     	nDirectPaths int
    	directPhotons, indirectPhotons, causticPhotons *[]*photon
    	radiancePhotons *[]*radiancePhoton
    	rpReflectances, rpTransmittances *[]Spectrum
    	nshot int
    	lightDistribution *Distribution1D
    	scene *Scene
    	renderer Renderer	
	}
	
	computeRadianceTask struct {
    	progress *ProgressReporter
    	taskNum, numTasks int
 	   	radiancePhotons *[]radiancePhoton
   		rpReflectances, rpTransmittances *[]Spectrum
     	nLookup int
    	maxDistSquared float64
    	nDirectPaths, nIndirectPaths, nCausticPaths int
    	directMap, indirectMap, causticMap *KdTree		
	}
	
	photonProcess struct {
		photons []closePhoton
		nLookup, nFound int
	}
	
	radiancePhotonProcess struct {
		n *Normal
		photon *radiancePhoton
	}
)


func NewPhotonIntegrator(nCaustic, nIndirect, nUsed, maxSpecularDepth, maxPhotonDepth int, maxDist float64, finalGather bool, gatherSamples int, gatherAngle float64) *PhotonIntegrator {
	integrator := new(PhotonIntegrator)
    integrator.nCausticPhotonsWanted = nCaustic
    integrator.nIndirectPhotonsWanted = nIndirect
    integrator.nLookup = nUsed
    integrator.maxSpecularDepth = maxSpecularDepth
    integrator.maxPhotonDepth = maxPhotonDepth
    integrator.maxDistSquared = maxDist * maxDist
    integrator.finalGather = finalGather
    integrator.cosGatherAngle = math.Cos(Radians(gatherAngle))
    integrator.gatherSamples = gatherSamples
    integrator.nCausticPaths = 0
    integrator.nIndirectPaths = 0
    integrator.causticMap = nil
    integrator.indirectMap = nil
    integrator.radianceMap = nil
    integrator.lightSampleOffsets = nil
    integrator.bsdfSampleOffsets = nil
	
	return integrator
}

func (integrator *PhotonIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer)    {
    if len(scene.lights) == 0 { return }
    // Declare shared variables for photon shooting
    nDirectPaths := 0
    directPhotons := make([]*photon, 0, 0)
    radiancePhotons := make([]*radiancePhoton, 0, 0)
    abortTasks := false
    causticPhotons := make([]*photon, 0, integrator.nCausticPhotonsWanted)
    indirectPhotons := make([]*photon, 0, integrator.nIndirectPhotonsWanted)
    nshot := 0
    rpReflectances := make([]Spectrum, 0, 0)
    rpTransmittances := make([]Spectrum, 0, 0)

    // Compute light power CDF for photon shooting
    lightDistribution := ComputeLightSamplingCDF(scene)

	shutterOpen := 0.0
	if camera != nil { shutterOpen = camera.ShutterOpen() }

    // Run parallel tasks for photon shooting
    progress := NewProgressReporter(integrator.nCausticPhotonsWanted+integrator.nIndirectPhotonsWanted, "Shooting photons", TerminalWidth())
    nTasks := NumSystemCores()
    for i := 0; i < nTasks; i++ {
        shooter := newPhotonShootingTask(
            i, shutterOpen, integrator, progress, abortTasks, nDirectPaths,
            directPhotons, indirectPhotons, causticPhotons, radiancePhotons,
            rpReflectances, rpTransmittances,
            nshot, lightDistribution, scene, renderer)
    	shooter.run()
    }		
    progress.Done()

    // Build kd-trees for indirect and caustic photons
    var directMap *KdTree = nil
    if len(directPhotons) > 0 {
    	directNodeData := make([]NodeData, len(directPhotons), len(directPhotons))
    	for i, p := range directPhotons {
    		directNodeData[i] = p
    	}
        directMap = NewKdTree(directNodeData)
    }
    if len(causticPhotons) > 0 {
    	causticNodeData := make([]NodeData, len(causticPhotons), len(causticPhotons))
    	for i, p := range causticPhotons {
    		causticNodeData[i] = p
    	}
        integrator.causticMap = NewKdTree(causticNodeData)
    }    
    if len(indirectPhotons) > 0 {
    	indirectNodeData := make([]NodeData, len(indirectPhotons), len(indirectPhotons))
    	for i, p := range indirectPhotons {
    		indirectNodeData[i] = p
    	}
        integrator.indirectMap = NewKdTree(indirectNodeData)
	}
    
    // Precompute radiance at a subset of the photons
    if integrator.finalGather && len(radiancePhotons) > 0 {
        // Launch tasks to compute photon radiances
        numTasks := 64
        progRadiance := NewProgressReporter(numTasks, "Computing photon radiances", TerminalWidth())
        for i := 0; i < numTasks; i++ {
            compRad := newComputeRadianceTask(progRadiance,
                i, numTasks, radiancePhotons, rpReflectances, rpTransmittances,
                integrator.nLookup, integrator.maxDistSquared, nDirectPaths, directMap,
                integrator.nIndirectPaths, integrator.indirectMap,
                integrator.nCausticPaths, integrator.causticMap)
            compRad.run()
        }    
        progRadiance.Done()
        radianceNodeData := make([]NodeData, len(radiancePhotons), len(radiancePhotons))
        for i, p := range radiancePhotons {
        	radianceNodeData[i] = p
        }
        integrator.radianceMap = NewKdTree(radianceNodeData)
    }
}

func (integrator *PhotonIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
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

    // Request samples for final gathering
    if integrator.finalGather {
        integrator.gatherSamples = Maxi(1, integrator.gatherSamples/2)
        if sampler != nil { integrator.gatherSamples = sampler.RoundSize(integrator.gatherSamples) }
        integrator.bsdfGatherSampleOffsets = *CreateBSDFSampleOffsets(integrator.gatherSamples, sample)
        integrator.indirGatherSampleOffsets = *CreateBSDFSampleOffsets(integrator.gatherSamples, sample)
    }
}

const (
	debugging = false
    nIndirSamplePhotons = 50	
)

func (integrator *PhotonIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
    L := NewSpectrum1(0.0)
    wo := ray.Dir.Negate()
    // Compute emitted light if ray hit an area light source
    L = L.Add(isect.Le(wo))

    // Evaluate BSDF at hit point
    bsdf := isect.GetBSDF(ray, arena)
    p := bsdf.dgShading.p
    n := bsdf.dgShading.nn
    L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.Time, bsdf, sample, rng,
        integrator.lightSampleOffsets, integrator.bsdfSampleOffsets))
    // Compute caustic lighting for photon map integrator
    lookupBuf := make([]closePhoton, integrator.nLookup, integrator.nLookup)
    L = L.Add(LPhoton(integrator.causticMap, integrator.nCausticPaths, integrator.nLookup, lookupBuf, bsdf,
                 rng, isect, wo, integrator.maxDistSquared))

    // Compute indirect lighting for photon map integrator
    if integrator.finalGather && integrator.indirectMap != nil {
    	if !debugging {
        // Do one-bounce final gather for photon map
        nonSpecular := BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY)
        if bsdf.NumComponentsMatching(nonSpecular) > 0 {
            // Find indirect photons around point for importance sampling
            proc := photonProcess{make([]closePhoton, nIndirSamplePhotons, nIndirSamplePhotons), nIndirSamplePhotons, 0}
            searchDist2 := integrator.maxDistSquared
            for proc.nFound < nIndirSamplePhotons {
                md2 := searchDist2
                proc.nFound = 0
                integrator.indirectMap.Lookup(p, proc.photonProc, &md2)
                searchDist2 *= 2.0
            }

            // Copy photon directions to local array
            photonDirs := make([]Vector, nIndirSamplePhotons, nIndirSamplePhotons) 
            for i := 0; i < nIndirSamplePhotons; i++ {
                photonDirs[i] = proc.photons[i].photon.wi
			}
            // Use BSDF to do final gathering
            Li := NewSpectrum1(0.0)
            for i := 0; i < integrator.gatherSamples; i++ {
                // Sample random direction from BSDF for final gather ray
                bsdfSample := CreateBSDFSample(sample, &integrator.bsdfGatherSampleOffsets, i)
                fr, wi, pdf, _ := bsdf.Sample_f(wo, bsdfSample, BxDFType(BSDF_ALL ^ BSDF_SPECULAR))
                if fr.IsBlack() || pdf == 0.0 { continue }
                Assert(pdf >= 0.0)

                // Trace BSDF final gather ray and accumulate radiance
                bounceRay := CreateChildRayDifferentialFromRayDifferential(p, wi, ray, isect.rayEpsilon, INFINITY)
                if ok, gatherIsect := scene.Intersect(bounceRay); ok {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Lindir := NewSpectrum1(0.0)
                    nGather := gatherIsect.dg.nn
                    nGather = FaceforwardNormalVector(nGather, bounceRay.Dir.Negate())
                    proc := radiancePhotonProcess{nGather, nil}
                    md2 := INFINITY
                    integrator.radianceMap.Lookup(gatherIsect.dg.p, proc.radiancePhotonProc, &md2)
                    if proc.photon != nil {
                        Lindir = &proc.photon.Lo
				    }                        
                    Lindir = Lindir.Mult(renderer.Transmittance(scene, bounceRay, nil, rng, arena))

                    // Compute MIS weight for BSDF-sampled gather ray

                    // Compute PDF for photon-sampling of direction _wi_
                    photonPdf := 0.0
                    conePdf := UniformConePdf(integrator.cosGatherAngle)
                    for j := 0; j < nIndirSamplePhotons; j++ {
                        if DotVector(&photonDirs[j], wi) > 0.999 * integrator.cosGatherAngle {
                            photonPdf += conePdf
                        }
                    }        
                    photonPdf /= nIndirSamplePhotons
                    wt := PowerHeuristic(integrator.gatherSamples, pdf, integrator.gatherSamples, photonPdf)
                    Li = Li.Add(fr.Mult(Lindir.Scale(AbsDotVectorNormal(wi, n) * wt / pdf)))
                }
            }
            L = L.Add(Li.InvScale(float64(integrator.gatherSamples)))

            // Use nearby photons to do final gathering
            Li = NewSpectrum1(0.0)
            for i := 0; i < integrator.gatherSamples; i++ {
                // Sample random direction using photons for final gather ray
                gatherSample := CreateBSDFSample(sample, &integrator.indirGatherSampleOffsets, i)
                photonNum := Mini(nIndirSamplePhotons - 1,
                    Floor2Int(gatherSample.uComponent * nIndirSamplePhotons))

                // Sample gather ray direction from _photonNum_
                vx, vy := CoordinateSystem(&photonDirs[photonNum])
                wi := UniformSampleConeVector(gatherSample.uDir[0], gatherSample.uDir[1],
                                              integrator.cosGatherAngle, vx, vy, &photonDirs[photonNum])

                // Trace photon-sampled final gather ray and accumulate radiance
                fr := bsdf.f(wo, wi, BSDF_ALL)
                if fr.IsBlack() { continue }
                bounceRay := CreateChildRayDifferentialFromRayDifferential(p, wi, ray, isect.rayEpsilon, INFINITY)
                //PBRT_PHOTON_MAP_STARTED_GATHER_RAY(&bounceRay);
                if ok, gatherIsect := scene.Intersect(bounceRay); ok {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Lindir := NewSpectrum1(0.0)
                    nGather := gatherIsect.dg.nn
                    nGather = FaceforwardNormal(nGather, CreateNormalFromVector(bounceRay.Dir.Negate()))
                    proc := radiancePhotonProcess{nGather, nil}
                    md2 := INFINITY
                    integrator.radianceMap.Lookup(gatherIsect.dg.p, proc.radiancePhotonProc, &md2)
                    if proc.photon != nil {
                        Lindir = &proc.photon.Lo
                    }    
                    Lindir = Lindir.Mult(renderer.Transmittance(scene, bounceRay, nil, rng, arena))

                    // Compute PDF for photon-sampling of direction _wi_
                    photonPdf := 0.0
                    conePdf := UniformConePdf(integrator.cosGatherAngle)
                    for j := 0; j < nIndirSamplePhotons; j++ {
                        if DotVector(&photonDirs[j], wi) > 0.999 * integrator.cosGatherAngle {
                            photonPdf += conePdf
                         }   
                    }        
                    photonPdf /= nIndirSamplePhotons

                    // Compute MIS weight for photon-sampled gather ray
                    bsdfPdf := bsdf.Pdf(wo, wi, BSDF_ALL)
                    wt := PowerHeuristic(integrator.gatherSamples, photonPdf, integrator.gatherSamples, bsdfPdf)
                    Li = Li.Add(fr.Mult(Lindir.Scale(AbsDotVectorNormal(wi, n) * wt / photonPdf)))
                }
                //PBRT_PHOTON_MAP_FINISHED_GATHER_RAY(&bounceRay);
            }
            L = L.Add(Li.InvScale(float64(integrator.gatherSamples)))
        }
    } else {
        // for debugging / examples: use the photon map directly
        nn := FaceforwardNormal(n, CreateNormalFromVector(ray.Dir.Negate()))
        proc := radiancePhotonProcess{nn, nil}
        md2 := INFINITY
        integrator.radianceMap.Lookup(p, proc.radiancePhotonProc, &md2)
        if proc.photon != nil {
            L = L.Add(&proc.photon.Lo)
        }    
    }
    } else {
        L = L.Add(LPhoton(integrator.indirectMap, integrator.nIndirectPaths, integrator.nLookup, lookupBuf,
                     bsdf, rng, isect, wo, integrator.maxDistSquared))
    }    
    if ray.Depth+1 < integrator.maxSpecularDepth {
        // Trace rays for specular reflection and refraction
        L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
        L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
    }
    return L
}
	
func (shooter *photonShootingTask) run() {
    // Declare local variables for _PhotonShootingTask_
    var arena *MemoryArena
    rng := NewRNG(int64(31 * shooter.taskNum))
    var localDirectPhotons, localIndirectPhotons, localCausticPhotons []photon
    var localRadiancePhotons []radiancePhoton
    totalPaths := 0
    causticDone := (shooter.integrator.nCausticPhotonsWanted == 0)
    indirectDone := (shooter.integrator.nIndirectPhotonsWanted == 0)
    halton := NewPermutedHalton(6, rng)
    var localRpReflectances, localRpTransmittances []Spectrum
    for {
        // Follow photon paths for a block of samples
        blockSize := 4096
        for i := 0; i < blockSize; i++ {
            u := halton.Sample(totalPaths)
            totalPaths++
            // Choose light to shoot photon from
            lightNum, lightPdf := shooter.lightDistribution.SampleDiscrete(u[0])
            light := shooter.scene.lights[lightNum]

            // Generate _photonRay_ from light source and initialize _alpha_
            ls := LightSample{[2]float64{u[1], u[2]}, u[3]}
            Le, pRay, Nl, pdf := light.Sample_L2(shooter.scene, &ls, u[4], u[5], shooter.time)
            photonRay := CreateRayDifferentialFromRay(pRay)
            if pdf == 0.0 || Le.IsBlack() { continue }
            alpha := Le.Scale(AbsDotNormalVector(Nl, &photonRay.Dir)).InvScale(pdf * lightPdf)
            if !alpha.IsBlack() {
                // Follow photon path through scene and record intersections
                //PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);
                specularPath := true
                nIntersections := 0
                for {
                	var ok bool
                	var photonIsect *Intersection
                	if ok, photonIsect = shooter.scene.Intersect(photonRay); !ok { break }
                    nIntersections++
                    // Handle photon/surface intersection
                    alpha = alpha.Mult(shooter.renderer.Transmittance(shooter.scene, photonRay, nil, rng, arena))
                    photonBSDF := photonIsect.GetBSDF(photonRay, arena)
                    specularType := BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)
                    hasNonSpecular := (photonBSDF.NumComponents() > photonBSDF.NumComponentsMatching(specularType))
                    wo := photonRay.Dir.Negate()
                    if hasNonSpecular {
                        // Deposit photon at surface
                        photon := photon{*photonIsect.dg.p, *alpha, *wo}
                        depositedPhoton := false
                        if specularPath && nIntersections > 1 {
                            if !causticDone {
                                //PBRT_PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true
                                localCausticPhotons = append(localCausticPhotons, photon)
                            }
                        } else {
                            // Deposit either direct or indirect photon
                            // stop depositing direct photons once indirectDone is true; don't
                            // want to waste memory storing too many if we're going a long time
                            // trying to get enough caustic photons desposited.
                            if nIntersections == 1 && !indirectDone && shooter.integrator.finalGather {
                                //PBRT_PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(&photonIsect.dg, &alpha, &wo)
                                depositedPhoton = true
                                localDirectPhotons = append(localDirectPhotons, photon)
                            } else if nIntersections > 1 && !indirectDone {
                                //PBRT_PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true
                                localIndirectPhotons = append(localIndirectPhotons, photon)
                            }
                        }

                        // Possibly create radiance photon at photon intersection point
                        if depositedPhoton && shooter.integrator.finalGather &&
                                rng.RandomFloat() < 0.125 {
                          	n := photonIsect.dg.nn
                            n = FaceforwardNormal(n, CreateNormalFromVector(photonRay.Dir.Negate()))
                            localRadiancePhotons = append(localRadiancePhotons, radiancePhoton{*photonIsect.dg.p, *n, *NewSpectrum1(0.0)})
                            rho_r := photonBSDF.rho(rng, BSDF_ALL_REFLECTION, DEFAULT_BSDF_SQRT_SAMPLES)
                            localRpReflectances = append(localRpReflectances, *rho_r)
                            rho_t := photonBSDF.rho(rng, BSDF_ALL_TRANSMISSION, DEFAULT_BSDF_SQRT_SAMPLES)
                            localRpTransmittances = append(localRpTransmittances, *rho_t)
                        }
                    }
                    if nIntersections >= shooter.integrator.maxPhotonDepth { break }

                    // Sample new photon ray direction
                    fr, wi, pdf, flags := photonBSDF.Sample_f(wo, CreateRandomBSDFSample(rng), BSDF_ALL)
                    if fr.IsBlack() || pdf == 0.0 { break }
                    anew := alpha.Mult(fr.Scale(AbsDotVectorNormal(wi, photonBSDF.dgShading.nn) / pdf))

                    // Possibly terminate photon path with Russian roulette
                    continueProb := math.Min(1.0, anew.Y() / alpha.Y())
                    if rng.RandomFloat() > continueProb {
                        break
                    }    
                    alpha = anew.InvScale(continueProb)
                    specularPath = specularPath && ((flags & BSDF_SPECULAR) != 0)
                    
                    if indirectDone && !specularPath { break }
                    photonRay = CreateChildRayDifferentialFromRayDifferential(photonIsect.dg.p, wi, photonRay,
                                                photonIsect.rayEpsilon, INFINITY)
                }
                //PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
            //arena.FreeAll();
        }

        // Merge local photon data with data in _PhotonIntegrator_
        { 

        // Give up if we're not storing enough photons
        if shooter.abortTasks {
            return
        }    
        if shooter.nshot > 500000 &&
            (unsuccessful(shooter.integrator.nCausticPhotonsWanted,
                                      len(*shooter.causticPhotons), blockSize) ||
             unsuccessful(shooter.integrator.nIndirectPhotonsWanted,
                                      len(*shooter.indirectPhotons), blockSize)) {
            Error("Unable to store enough photons.  Giving up.\n")
            *shooter.causticPhotons = make([]*photon, 0, 0)
            *shooter.indirectPhotons = make([]*photon, 0, 0)
            *shooter.radiancePhotons = make([]*radiancePhoton, 0, 0)
            shooter.abortTasks = true
            return
        }
        shooter.progress.Update(len(localIndirectPhotons) + len(localCausticPhotons))
        shooter.nshot += blockSize

        // Merge indirect photons into shared array
        if !indirectDone {
            shooter.integrator.nIndirectPaths += blockSize
            for i := 0; i < len(localIndirectPhotons); i++ {
                *shooter.indirectPhotons = append(*shooter.indirectPhotons, &localIndirectPhotons[i])
            }    
            localIndirectPhotons = make([]photon, 0, 0)
            if len(*shooter.indirectPhotons) >= shooter.integrator.nIndirectPhotonsWanted {
                indirectDone = true
            }    
            shooter.nDirectPaths += blockSize
            for i := 0; i < len(localDirectPhotons); i++ {
                *shooter.directPhotons = append(*shooter.directPhotons, &localDirectPhotons[i])
            }    
            localDirectPhotons = make([]photon, 0, 0)
        }

        // Merge direct, caustic, and radiance photons into shared array
        if !causticDone {
            shooter.integrator.nCausticPaths += blockSize
            for i := 0; i < len(localCausticPhotons); i++ {
                *shooter.causticPhotons = append(*shooter.causticPhotons, &localCausticPhotons[i])
            }    
            localCausticPhotons = make([]photon, 0, 0)
            if len(*shooter.causticPhotons) >= shooter.integrator.nCausticPhotonsWanted {
                causticDone = true
            }    
        }
        
        for i := 0; i < len(localRadiancePhotons); i++ {
            *shooter.radiancePhotons = append(*shooter.radiancePhotons, &localRadiancePhotons[i])
        }    
        localRadiancePhotons = make([]radiancePhoton, 0, 0)
        for i := 0; i < len(localRpReflectances); i++ {
            *shooter.rpReflectances = append(*shooter.rpReflectances, localRpReflectances[i])
        }    
        localRpReflectances = make([]Spectrum, 0, 0)
        for i := 0; i < len(localRpTransmittances); i++ {
            *shooter.rpTransmittances = append(*shooter.rpTransmittances, localRpTransmittances[i])
        }    
        localRpTransmittances = make([]Spectrum, 0, 0)
        }

        // Exit task if enough photons have been found
        if indirectDone && causticDone {
            break
        }    
    }
}

func (compRad *computeRadianceTask) run() {
    // Compute range of radiance photons to process in task
    taskSize := len(*compRad.radiancePhotons) / compRad.numTasks
    excess := len(*compRad.radiancePhotons) % compRad.numTasks
    rpStart := Mini(compRad.taskNum, excess) * (taskSize+1) +
                       Maxi(0, compRad.taskNum-excess) * taskSize
    rpEnd := rpStart + taskSize
    if compRad.taskNum < excess { rpEnd++ }
    if compRad.taskNum == compRad.numTasks-1 { Assert(rpEnd == len(*compRad.radiancePhotons)) }
    lookupBuf := make([]closePhoton, compRad.nLookup, compRad.nLookup)
    for  i := rpStart; i < rpEnd; i++ {
        // Compute radiance for radiance photon _i_
        rp := (*compRad.radiancePhotons)[i]
        rho_r := &(*compRad.rpReflectances)[i]
        rho_t := &(*compRad.rpTransmittances)[i]
        if !rho_r.IsBlack() {
            // Accumulate outgoing radiance due to reflected irradiance
            E := EPhoton(compRad.directMap, compRad.nDirectPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, &rp.n).Add(
                         EPhoton(compRad.indirectMap, compRad.nIndirectPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, &rp.n)).Add(
                         EPhoton(compRad.causticMap, compRad.nCausticPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, &rp.n))
            rp.Lo = *rp.Lo.Add(E.Mult(rho_r.Scale(INV_PI)))
        }
        if !rho_t.IsBlack() {
            // Accumulate outgoing radiance due to transmitted irradiance
             E := EPhoton(compRad.directMap, compRad.nDirectPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, rp.n.Negate()).Add(
                         EPhoton(compRad.indirectMap, compRad.nIndirectPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, rp.n.Negate())).Add(
                         EPhoton(compRad.causticMap, compRad.nCausticPaths, compRad.nLookup, lookupBuf,
                                 compRad.maxDistSquared, &rp.p, rp.n.Negate()))
            rp.Lo = *rp.Lo.Add(E.Mult(rho_t.Scale(INV_PI)))
        }
    }
    compRad.progress.Update(1)
}

func unsuccessful(needed, found, shot int) bool {
    return (found < needed && (found == 0 || found < shot / 1024))
}

func kernel(photon *photon, p *Point, maxDist2 float64) float64 {
    s := (1.0 - DistanceSquaredPoint(&photon.p, p) / maxDist2)
    return 3.0 * INV_PI * s * s
}


func LPhoton(pmap *KdTree, nPaths, nLookup int,
      lookupBuf []closePhoton, bsdf *BSDF, rng *RNG,
      isect *Intersection, wo *Vector, maxDist2 float64) *Spectrum {
    L := NewSpectrum1(0.0)
    nonSpecular := BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY)
    if pmap != nil && bsdf.NumComponentsMatching(nonSpecular) > 0 {
        //PBRT_PHOTON_MAP_STARTED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg));
        // Do photon map lookup at intersection point
        proc := photonProcess{lookupBuf, nLookup, 0}
        pmap.Lookup(isect.dg.p, proc.photonProc, &maxDist2)

        // Estimate reflected radiance due to incident photons
        photons := proc.photons
        nFound := proc.nFound
        Nf := FaceforwardNormalVector(bsdf.dgShading.nn, wo)
        if bsdf.NumComponentsMatching(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0 {
            // Compute exitant radiance from photons for glossy surface
            for i := 0; i < nFound; i++ {
                p := photons[i].photon
                k := kernel(p, isect.dg.p, maxDist2)
                L = L.Add(bsdf.f(wo, &p.wi, BSDF_ALL).Scale(k / (float64(nPaths) * maxDist2)).Mult(&p.alpha))
            }
        } else {
            // Compute exitant radiance from photons for diffuse surface
            Lr := NewSpectrum1(0.0)
            Lt := NewSpectrum1(0.0)
            for i := 0; i < nFound; i++ {
                if DotNormalVector(Nf, &photons[i].photon.wi) > 0.0 {
                    k := kernel(photons[i].photon, isect.dg.p, maxDist2)
                    Lr = Lr.Add(photons[i].photon.alpha.Scale(k / (float64(nPaths) * maxDist2)))
                } else {
                    k := kernel(photons[i].photon, isect.dg.p, maxDist2)
                    Lt = Lt.Add(photons[i].photon.alpha.Scale(k / (float64(nPaths) * maxDist2)))
                }
            }
            L = L.Add(Lr.Mult(bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION, DEFAULT_BSDF_SQRT_SAMPLES).Scale(INV_PI)).Add(
                 Lt.Mult(bsdf.rho2(wo, rng, BSDF_ALL_TRANSMISSION, DEFAULT_BSDF_SQRT_SAMPLES).Scale(INV_PI))))
        }
        //PBRT_PHOTON_MAP_FINISHED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg), proc.nFound, proc.nLookup, &L)
    }
    return L;
}


func EPhoton(pmap *KdTree, count, nLookup int, lookupBuf []closePhoton, maxDist2 float64, p *Point, n *Normal) *Spectrum {
    if pmap == nil { return NewSpectrum1(0.0) }
    // Lookup nearby photons at irradiance computation point
    proc := photonProcess{lookupBuf, nLookup, 0}
    md2 := maxDist2
    pmap.Lookup(p, proc.photonProc, &md2)
    Assert(md2 > 0.0)

    // Accumulate irradiance value from nearby photons
    if proc.nFound == 0 { return NewSpectrum1(0.0) }
    photons := proc.photons
    E := NewSpectrum1(0.0)
    for i := 0; i < proc.nFound; i++ {
        if DotNormalVector(n, &photons[i].photon.wi) > 0.0 {
            E = E.Add(&photons[i].photon.alpha)
        }    
    }        
    return E.InvScale(float64(count) * md2 * math.Pi)
}

func (p *photon) Location() *Point {
	return &p.p
}
func (p *radiancePhoton) Location() *Point {
	return &p.p
}

func (proc *photonProcess) photonProc(p *Point, nodeData NodeData, dist2 float64, maxDistSquared *float64) {
	Unimplemented()
/*	
    if (proc.nFound < proc.nLookup) {
        // Add photon to unordered array of photons
        photons[nFound++] = ClosePhoton(&photon, distSquared)
        if proc.nFound == proc.nLookup {
            std::make_heap(&photons[0], &photons[nLookup]);
            maxDistSquared = photons[0].distanceSquared;
        }
    } else {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&photons[0], &photons[nLookup]);
        photons[nLookup-1] = ClosePhoton(&photon, distSquared);
        std::push_heap(&photons[0], &photons[nLookup]);
        maxDistSquared = photons[0].distanceSquared;
    }
*/    
}

func (proc *radiancePhotonProcess) radiancePhotonProc(p *Point, nodeData NodeData, dist2 float64, maxDistSquared *float64) {
	rp, ok := nodeData.(*radiancePhoton)
	if ok {
        if DotNormal(&rp.n, proc.n) > 0 {
            proc.photon = rp
            *maxDistSquared = dist2
        }
    }    
}
	
func newPhotonShootingTask(taskNum int, time float64, integrator *PhotonIntegrator, prog *ProgressReporter, abortTasks bool,
            nDirectPaths int, directPhotons, indirectPhotons, causticPhotons []*photon, radiancePhotons []*radiancePhoton,
            rpReflectances, rpTransmittances []Spectrum,
            nshot int, lightDistribution *Distribution1D, scene *Scene, renderer Renderer) *photonShootingTask {
    shooter := new(photonShootingTask)
    shooter.taskNum = taskNum
    shooter.time = time
    shooter.integrator = integrator
    shooter.progress = prog
    shooter.abortTasks = abortTasks
    shooter.nDirectPaths = nDirectPaths
    shooter.directPhotons = &directPhotons
    shooter.indirectPhotons = &indirectPhotons
    shooter.causticPhotons = &causticPhotons
    shooter.radiancePhotons = &radiancePhotons
    shooter.nshot = nshot
    shooter.lightDistribution = lightDistribution
    shooter.scene = scene
    shooter.renderer = renderer
    return shooter     	
}
            
func newComputeRadianceTask(progRadiance *ProgressReporter,
                taskNum, numTasks int , radiancePhotons []*radiancePhoton, rpReflectances, rpTransmittances []Spectrum,
                nLookup int, maxDistSquared float64, nDirectPaths int, directMap *KdTree,
                nIndirectPaths int, indirectMap *KdTree,
                nCausticPaths int, causticMap *KdTree) *computeRadianceTask {
                	return nil
                }            