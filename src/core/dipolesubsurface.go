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
	DipoleSubsurfaceIntegrator struct {
		maxSpecularDepth        int
		maxError, minSampleDist float64
		filename                string
		irradiancePoints        []IrradiancePoint
		octreeBounds            BBox
		octree                  []SubsurfaceOctreeNode
		octreeArena             *MemoryArena

		// Declare sample parameters for light source sampling
		lightSampleOffsets []LightSampleOffsets
		bsdfSampleOffsets  []BSDFSampleOffsets
	}

	// DipoleSubsurfaceIntegrator Helper Declarations
	IrradiancePoint struct {
		p                Point
		n                Normal
		E                Spectrum
		area, rayEpsilon float64
	}

	SubsurfaceOctreeNode struct {
		p       Point
		isLeaf  bool
		E       Spectrum
		sumArea float64
		// 'union'
		children [8]*SubsurfaceOctreeNode // isLeaf == false
		ips      [8]*IrradiancePoint      // isLeaf == true
	}

	DiffusionReflectance struct {
		zpos, zneg, sigmap_t, sigma_tr, alphap Spectrum
		A                                      float64
	}
)

func NewSubsurfaceOctreeNode() *SubsurfaceOctreeNode {
	node := new(SubsurfaceOctreeNode)
	node.isLeaf = true
	node.sumArea = 0.0
	return node
}

func (node *SubsurfaceOctreeNode) Insert(nodeBound *BBox, ip *IrradiancePoint, arena *MemoryArena) {
	pMid := nodeBound.pMin.AddPoint(&nodeBound.pMax).Scale(0.5)
	if node.isLeaf {
		// Add _IrradiancePoint_ to leaf octree node
		for i := 0; i < len(node.ips); i++ {
			if node.ips[i] == nil {
				node.ips[i] = ip
				return
			}
		}

		// Convert leaf node to interior node, redistribute points
		node.isLeaf = false
		var localIps [8]*IrradiancePoint
		for i := 0; i < len(node.ips); i++ {
			localIps[i] = node.ips[i]
			node.children[i] = nil
		}
		for _, ip := range localIps {
			// Add _IrradiancePoint_ _ip_ to interior octree node
			child := 0
			if ip.p.x > pMid.x {
				child += 4
			}
			if ip.p.y > pMid.y {
				child += 2
			}
			if ip.p.z > pMid.z {
				child += 1
			}
			if node.children[child] == nil {
				node.children[child] = NewSubsurfaceOctreeNode()
			}
			childBound := octreeChildBound(child, nodeBound, pMid)
			node.children[child].Insert(&childBound, ip, arena)
		}
		/* fall through to interior case to insert the new point... */
	}
	// Add _IrradiancePoint_ _ip_ to interior octree node
	child := 0
	if ip.p.x > pMid.x {
		child += 4
	}
	if ip.p.y > pMid.y {
		child += 2
	}
	if ip.p.z > pMid.z {
		child += 1
	}
	if node.children[child] == nil {
		node.children[child] = NewSubsurfaceOctreeNode()
	}
	childBound := octreeChildBound(child, nodeBound, pMid)
	node.children[child].Insert(&childBound, ip, arena)
}

func (node *SubsurfaceOctreeNode) InitHierarchy() {
	if node.isLeaf {
		// Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
		sumWt := 0.0
		i := 0
		for i = 0; i < len(node.ips); i++ {
			if node.ips[i] == nil {
				break
			}
			wt := node.ips[i].E.Y()
			node.E = *node.E.Add(&node.ips[i].E)
			node.p = *node.p.AddPoint(node.ips[i].p.Scale(wt))
			sumWt += wt
			node.sumArea += node.ips[i].area
		}
		if sumWt > 0.0 {
			node.p = *node.p.InvScale(sumWt)
		}
		node.E = *node.E.InvScale(float64(i))
	} else {
		// Init interior _SubsurfaceOctreeNode_
		sumWt := 0.0
		nChildren := 0
		for i := 0; i < len(node.children); i++ {
			if node.children[i] == nil {
				continue
			}
			nChildren++
			node.children[i].InitHierarchy()
			wt := node.children[i].E.Y()
			node.E = *node.E.Add(&node.children[i].E)
			node.p = *node.p.AddPoint(node.children[i].p.Scale(wt))
			sumWt += wt
			node.sumArea += node.children[i].sumArea
		}
		if sumWt > 0.0 {
			node.p = *node.p.InvScale(sumWt)
		}
		node.E = *node.E.InvScale(float64(nChildren))
	}
}

func (node *SubsurfaceOctreeNode) Mo(nodeBound *BBox, pt *Point, Rd *DiffusionReflectance, maxError float64) *Spectrum {
	// Compute $M_\roman{o}$ at node if error is low enough
	dw := node.sumArea / DistanceSquaredPoint(pt, &node.p)
	if dw < maxError && !nodeBound.Inside(pt) {
		//PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<SubsurfaceOctreeNode *>(this))
		return Rd.Evaluate(DistanceSquaredPoint(pt, &node.p)).Mult(node.E.Scale(node.sumArea))
	}

	// Otherwise compute $M_\roman{o}$ from points in leaf or recursively visit children
	mo := NewSpectrum1(0.0)
	if node.isLeaf {
		// Accumulate $M_\roman{o}$ from leaf node
		for i := 0; i < len(node.ips); i++ {
			if node.ips[i] == nil {
				break
			}
			//PBRT_SUBSURFACE_ADDED_POINT_CONTRIBUTION(const_cast<IrradiancePoint *>(ips[i]))
			mo = mo.Add(Rd.Evaluate(DistanceSquaredPoint(pt, &node.ips[i].p)).Mult(node.ips[i].E.Scale(node.ips[i].area)))
		}
	} else {
		// Recursively visit children nodes to compute $M_\roman{o}$
		pMid := nodeBound.pMin.AddPoint(&nodeBound.pMax).Scale(0.5)
		for child := 0; child < len(node.children); child++ {
			if node.children[child] == nil {
				continue
			}
			childBound := octreeChildBound(child, nodeBound, pMid)
			mo = mo.Add(node.children[child].Mo(&childBound, pt, Rd, maxError))
		}
	}
	return mo
}

func NewDiffusionReflectance(sigma_a, sigmap_s *Spectrum, eta float64) *DiffusionReflectance {
	dr := new(DiffusionReflectance)
	dr.A = (1.0 + Fdr(eta)) / (1.0 - Fdr(eta))
	dr.sigmap_t = *sigma_a.Add(sigmap_s)
	dr.sigma_tr = *SqrtSpectrum(sigma_a.Mult(&dr.sigmap_t).Scale(3.0))
	dr.alphap = *sigmap_s.Divide(&dr.sigmap_t)
	dr.zpos = *NewSpectrum1(1.0).Divide(&dr.sigmap_t)
	dr.zneg = *(dr.zpos.Scale(1.0 + (4.0/3.0)*dr.A)).Negate()
	return dr
}

func (dr *DiffusionReflectance) Evaluate(d2 float64) *Spectrum {
	dpos := SqrtSpectrum(NewSpectrum1(d2).Add(dr.zpos.Mult(&dr.zpos)))
	dneg := SqrtSpectrum(NewSpectrum1(d2).Add(dr.zneg.Mult(&dr.zneg)))
	Rd := (dr.alphap.InvScale(4.0 * math.Pi)).Mult((dr.zpos.Mult(dpos.Mult(&dr.sigma_tr).Add(NewSpectrum1(1.0))).Mult(ExpSpectrum(dr.sigma_tr.Negate().Mult(dpos))).Divide(dpos.Mult(dpos.Mult(dpos)))).Sub(dr.zneg.Mult(dneg.Mult(&dr.sigma_tr).Add(NewSpectrum1(1.0))).Mult(ExpSpectrum(dr.sigma_tr.Negate().Mult(dneg))).Divide(dneg.Mult(dneg.Mult(dneg)))))
	return Rd.Clamp(0.0, INFINITY)
}

func NewIrradiancePoint(sp *SurfacePoint, ee *Spectrum) *IrradiancePoint {
	return &IrradiancePoint{p: sp.p, n: sp.n, E: *ee, area: sp.area, rayEpsilon: sp.rayEpsilon}
}

func NewDipoleSubsurfaceIntegrator(maxDepth int, maxError, minDist float64, pointsfile string) *DipoleSubsurfaceIntegrator {
	dipole := new(DipoleSubsurfaceIntegrator)
	dipole.maxSpecularDepth = maxDepth
	dipole.maxError = maxError
	dipole.minSampleDist = minDist
	dipole.filename = pointsfile

	return dipole
}

func (integrator *DipoleSubsurfaceIntegrator) Preprocess(scene *Scene, camera Camera, renderer Renderer) {
	if len(scene.lights) == 0 {
		return
	}
	pts := make([]SurfacePoint, 0, 16)
	// Get _SurfacePoint_s for translucent objects in scene
	if len(integrator.filename) != 0 {
		// Initialize _SurfacePoint_s from file
		var fpts []float64
		var ok bool
		if ok, fpts = ReadFloatFile(integrator.filename); ok {
			if (len(fpts) % 8) != 0 {
				Error("Excess values (%d) in points file \"%s\"", len(fpts)%8, integrator.filename)
			}
			for i := 0; i < len(fpts); i = i + 8 {
				pts = append(pts, SurfacePoint{Point{fpts[i], fpts[i+1], fpts[i+2]},
					Normal{fpts[i+3], fpts[i+4], fpts[i+5]},
					fpts[i+6], fpts[i+7]})
			}
		}
	}
	if len(pts) == 0 {
		pCamera := PointAnimatedTransform(camera.CameraToWorld(), camera.ShutterOpen(), CreatePoint(0, 0, 0))
		FindPoissonPointDistribution(pCamera, camera.ShutterOpen(), integrator.minSampleDist, scene, &pts)
	}

	// Compute irradiance values at sample points
	rng := NewRNG(10)
	var arena *MemoryArena
	//PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
	progress := NewProgressReporter(len(pts), "Computing Irradiances", TerminalWidth())
	for i := 0; i < len(pts); i++ {
		sp := &pts[i]
		E := NewSpectrum1(0.0)
		for j := 0; j < len(scene.lights); j++ {
			// Add irradiance from light at point
			light := scene.lights[j]
			Elight := NewSpectrum1(0.0)
			nSamples := RoundUpPow2(uint32(light.NumSamples()))
			scramble := [2]uint32{rng.RandomUInt(), rng.RandomUInt()}
			compScramble := rng.RandomUInt()
			var s uint32
			for s = 0; s < nSamples; s++ {
				lpos := [2]float64{0.0, 0.0}
				lpos[0], lpos[1] = Sample02(s, scramble)
				lcomp := VanDerCorput(s, compScramble)
				ls := &LightSample{lpos, lcomp}
				var vis VisibilityTester
				Li, wi, lightPdf := light.Sample_L(&sp.p, sp.rayEpsilon, ls, camera.ShutterOpen(), &vis)
				if DotVectorNormal(wi, &sp.n) <= 0.0 {
					continue
				}
				if Li.IsBlack() || lightPdf == 0.0 {
					continue
				}
				Li = Li.Mult(vis.Transmittance(scene, renderer, nil, rng, arena))
				if vis.Unoccluded(scene) {
					Elight = Elight.Add(Li.Scale(AbsDotVectorNormal(wi, &sp.n) / lightPdf))
				}
			}
			E = E.Add(Elight.InvScale(float64(nSamples)))
		}
		integrator.irradiancePoints = append(integrator.irradiancePoints, *NewIrradiancePoint(sp, E))
		//PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
		//arena.FreeAll();
		progress.Update(1)
	}
	progress.Done()
	//PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

	// Create octree of clustered irradiance samples
	integrator.octree = make([]SubsurfaceOctreeNode, 1, 1)
	integrator.octree[0] = *NewSubsurfaceOctreeNode()
	for _, ip := range integrator.irradiancePoints {
		integrator.octreeBounds = *UnionBBoxPoint(&integrator.octreeBounds, &ip.p)
	}
	for _, ip := range integrator.irradiancePoints {
		integrator.octree[0].Insert(&integrator.octreeBounds, &ip, nil)
	}
	integrator.octree[0].InitHierarchy()
}

func (integrator *DipoleSubsurfaceIntegrator) RequestSamples(sampler Sampler, sample *Sample, scene *Scene) {
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
}

func (integrator *DipoleSubsurfaceIntegrator) Li(scene *Scene, renderer Renderer, ray *RayDifferential, isect *Intersection,
	sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	L := NewSpectrum1(0.0)
	wo := ray.dir.Negate()
	// Compute emitted light if ray hit an area light source
	L = L.Add(isect.Le(wo))

	// Evaluate BSDF at hit point
	bsdf := isect.GetBSDF(ray, arena)
	p := bsdf.dgShading.p
	n := bsdf.dgShading.nn
	// Evaluate BSSRDF and possibly compute subsurface scattering
	bssrdf := isect.GetBSSRDF(ray, arena)
	if bssrdf != nil && integrator.octree != nil {
		sigma_a := bssrdf.sigma_a
		sigmap_s := bssrdf.sigma_prime_s
		sigmap_t := sigmap_s.Add(sigma_a)
		if !sigmap_t.IsBlack() {
			// Use hierarchical integration to evaluate reflection from dipole model
			//PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
			Rd := NewDiffusionReflectance(sigma_a, sigmap_s, bssrdf.eta)
			Mo := integrator.octree[0].Mo(&integrator.octreeBounds, p, Rd, integrator.maxError)
			fresnel := &FresnelDielectric{1.0, bssrdf.eta}
			Ft := NewSpectrum1(1.0).Sub(fresnel.Evaluate(AbsDotVectorNormal(wo, n)))
			Fdt := 1.0 - Fdr(bssrdf.eta)
			L = L.Add(Ft.InvScale(math.Pi).Mult(Mo.Scale(Fdt)))
			//PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
		}
	}
	L = L.Add(UniformSampleAllLights(scene, renderer, arena, p, n,
		wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, integrator.lightSampleOffsets,
		integrator.bsdfSampleOffsets))
	if ray.depth < integrator.maxSpecularDepth {
		// Trace rays for specular reflection and refraction
		L = L.Add(SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample, arena))
		L = L.Add(SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample, arena))
	}
	return L
}

func CreateDipoleSubsurfaceIntegrator(params *ParamSet) *DipoleSubsurfaceIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	maxError := params.FindFloatParam("maxerror", 0.05)
	minDist := params.FindFloatParam("minsampledistance", 0.25)
	pointsfile := params.FindFilenameParam("pointsfile", "")
	if options.QuickRender {
		maxError *= 4.0
		minDist *= 4.0
	}
	return NewDipoleSubsurfaceIntegrator(maxDepth, maxError, minDist, pointsfile)
}
