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
	"fmt"
	"os"
)

type (
	CreateRadianceProbes struct {
		surfaceIntegrator                              SurfaceIntegrator
		volumeIntegrator                               VolumeIntegrator
		camera                                         Camera
		lmax, nIndirSamples                            int
		bbox                                           BBox
		includeDirectInProbes, includeIndirectInProbes bool
		time, probeSpacing                             float64
		filename                                       string
	}

	createRadProbeTask struct {
		pointNum                                       int
		nProbes                                        [3]int
		bbox                                           BBox
		lmax, nIndirSamples                            int
		time                                           float64
		progress                                       *ProgressReporter
		includeDirectInProbes, includeIndirectInProbes bool
		origSample                                     *Sample
		renderer                                       Renderer
		scene                                          *Scene
		surfacePoints                                  []Point
		c_in                                           []Spectrum
	}
)

func NewCreateRadianceProbes(surf SurfaceIntegrator, vol VolumeIntegrator, camera Camera, lmax int, probeSpacing float64, bbox *BBox, nindir int,
	includeDirect, includeIndirect bool, time float64, filename string) *CreateRadianceProbes {
	probes := new(CreateRadianceProbes)
	probes.surfaceIntegrator = surf
	probes.volumeIntegrator = vol
	probes.camera = camera
	probes.lmax = lmax
	probes.nIndirSamples = nindir
	probes.bbox = *bbox
	probes.includeDirectInProbes = includeDirect
	probes.includeIndirectInProbes = includeIndirect
	probes.time = time
	probes.probeSpacing = probeSpacing
	probes.filename = filename
	return probes
}

func (probes *CreateRadianceProbes) Render(scene *Scene) {
	// Compute scene bounds and initialize probe integrators
	if probes.bbox.PMin.X > probes.bbox.PMax.X {
		probes.bbox = *scene.WorldBound()
	}
	probes.surfaceIntegrator.Preprocess(scene, probes.camera, probes)
	probes.volumeIntegrator.Preprocess(scene, probes.camera, probes)
	origSample := NewSample(nil, probes.surfaceIntegrator, probes.volumeIntegrator, scene)

	// Compute sampling rate in each dimension
	delta := probes.bbox.PMax.Sub(&probes.bbox.PMin)
	var nProbes [3]int
	for i := 0; i < 3; i++ {
		nProbes[i] = Maxi(1, Ceil2Int(delta.At(i)/probes.probeSpacing))
	}
	// Allocate SH coefficient vector pointers for sample points
	count := nProbes[0] * nProbes[1] * nProbes[2]
	c_in := make([][]Spectrum, count, count)
	for i := 0; i < count; i++ {
		c_in[i] = make([]Spectrum, SHTerms(probes.lmax), SHTerms(probes.lmax))
	}
	// Compute random points on surfaces of scene

	// Create scene bounding sphere to catch rays that leave the scene
	sceneCenter, sceneRadius := scene.WorldBound().BoundingSphere()
	ObjectToWorld := TranslateTransform(sceneCenter.Sub(CreatePoint(0, 0, 0)))
	WorldToObject := InverseTransform(ObjectToWorld)
	sph := CreateSphere(ObjectToWorld, WorldToObject, true, sceneRadius, -sceneRadius, sceneRadius, 360.0)
	sphere := NewGeometricPrimitive(sph, nil, nil)
	nPoints := 32768
	maxDepth := 32
	surfacePoints := make([]Point, 0, nPoints+maxDepth)
	pCamera := PointAnimatedTransform(probes.camera.CameraToWorld(), probes.camera.ShutterOpen(), CreatePoint(0, 0, 0))
	surfacePoints = append(surfacePoints, *pCamera)
	rng := NewRNG(13)
	for len(surfacePoints) < nPoints {
		// Generate random path from camera and deposit surface points
		pray := pCamera
		dir := UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat())
		rayEpsilon := 0.0
		for i := 0; i < maxDepth; i++ {
			ray := CreateRayDifferential(pray, dir, rayEpsilon, INFINITY, probes.time, 0)

			var isect *Intersection
			var hit bool
			if hit, isect = scene.Intersect(ray); !hit {
				if hit, isect = sphere.Intersect(ray); !hit {
					break
				}
			}
			surfacePoints = append(surfacePoints, *ray.PointAt(ray.Maxt()))

			hitGeometry := isect.dg
			pray = isect.dg.p
			rayEpsilon = isect.rayEpsilon
			hitGeometry.nn = FaceforwardNormalVector(hitGeometry.nn, ray.Dir().Negate())

			dir = UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat())
			dir = FaceforwardVectorNormal(dir, hitGeometry.nn)
		}
	}

	nWorkers := Maxi(NumSystemCores(), 1)
	jobs := make(chan *createRadProbeTask, count)
	completed := make(chan bool, count)

	for w := 0; w < nWorkers; w++ {
		go radProbeWorker(jobs, completed)
	}

	// Launch tasks to compute radiance probes at sample points
	prog := NewProgressReporter(count, "Radiance Probes", -1)
	for i := 0; i < count; i++ {
		task := newCreateRadProbeTask(i, nProbes, probes.time,
			&probes.bbox, probes.lmax, probes.includeDirectInProbes,
			probes.includeIndirectInProbes, probes.nIndirSamples,
			prog, origSample, surfacePoints,
			scene, probes, &c_in[i])
		jobs <- task
	}
	close(jobs)
	numCompleted := 0
	for numCompleted < nWorkers {
		select {
		case done := <-completed:
			if done {
				numCompleted++
			}
		default:
			break
		}
	}
	prog.Done()

	// Write radiance probe coefficients to file
	f, err := os.Open(probes.filename)
	defer f.Close()
	if err == nil {
		incDir := 0
		if probes.includeDirectInProbes {
			incDir = 1
		}
		incIndir := 0
		if probes.includeIndirectInProbes {
			incIndir = 1
		}
		_, err = fmt.Fprintf(f, "%d %d %d\n", probes.lmax, incDir, incIndir)
		_, err = fmt.Fprintf(f, "%d %d %d\n", nProbes[0], nProbes[1], nProbes[2])
		_, err = fmt.Fprintf(f, "%f %f %f %f %f %f\n", probes.bbox.PMin.X, probes.bbox.PMin.Y, probes.bbox.PMin.Z,
			probes.bbox.PMax.X, probes.bbox.PMax.Y, probes.bbox.PMax.Z)
		if err != nil {
			Severe("Error writing radiance file \"%s\" (%s)", probes.filename, err)
		}

		for i := 0; i < nProbes[0]*nProbes[1]*nProbes[2]; i++ {
			for j := 0; j < SHTerms(probes.lmax); j++ {
				_, err = fmt.Fprintf(f, "  ")
				if c_in[i][j].Write(f) == false {
					Severe("Error writing radiance file \"%s\" (%s)", probes.filename, err)
				}
				fmt.Fprintf(f, "\n")
			}
			fmt.Fprintf(f, "\n")
		}
	}
}

func (probes *CreateRadianceProbes) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (*Spectrum, *Intersection, *Spectrum) {
	Assert(ray.Time() == sample.time)
	Assert(!ray.HasNaNs())
	Lo := NewSpectrum1(0.0)
	var hit bool
	var isect *Intersection
	if hit, isect = scene.Intersect(ray); hit {
		Lo = probes.surfaceIntegrator.Li(scene, probes, ray, isect, sample, rng, arena)
	} else {
		for i := 0; i < len(scene.lights); i++ {
			Lo = Lo.Add(scene.lights[i].Le(ray))
		}
	}
	Lv, T := probes.volumeIntegrator.Li(scene, probes, ray, sample, rng, arena)
	li := T.Mult(Lo).Add(Lv)
	return li, isect, T
}

func (probes *CreateRadianceProbes) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return probes.volumeIntegrator.Transmittance(scene, probes, ray, sample, rng, arena)
}

func CreateRadianceProbesRenderer(camera Camera, surf SurfaceIntegrator, vol VolumeIntegrator, params *ParamSet) *CreateRadianceProbes {
	includeDirect := params.FindBoolParam("directlighting", true)
	includeIndirect := params.FindBoolParam("indirectlighting", true)
	lmax := params.FindIntParam("lmax", 4)
	nindir := params.FindIntParam("indirectsamples", 512)

	var bounds *BBox
	b := params.FindFloatArrayParam("bounds")
	if b != nil {
		if len(b) != 6 {
			Warning("Expecting six values [x0 y0 z0 x1 y1 z1] for bounds")
		} else {
			bounds = CreateBBoxFromPoints(CreatePoint(b[0], b[1], b[2]), CreatePoint(b[3], b[4], b[5]))
		}
	}
	probeSpacing := params.FindFloatParam("samplespacing", 1.0)
	time := params.FindFloatParam("time", 0.0)
	filename := params.FindFilenameParam("filename", "probes.out")

	return NewCreateRadianceProbes(surf, vol, camera, lmax, probeSpacing,
		bounds, nindir, includeDirect, includeIndirect, time, filename)
}

func newCreateRadProbeTask(taskNum int, nProbes [3]int, time float64, bbox *BBox, lmax int, incDir, incIndir bool, nIndir int,
	progress *ProgressReporter, sample *Sample, surfPoints []Point, scene *Scene, renderer Renderer, c_in *[]Spectrum) *createRadProbeTask {
	return nil
}

func radProbeWorker(workQueue <-chan *createRadProbeTask, completed chan<- bool) {

	for t := range workQueue {
		// Compute region in which to compute incident radiance probes
		sx := t.pointNum % t.nProbes[0]
		sy := (t.pointNum / t.nProbes[0]) % t.nProbes[1]
		sz := t.pointNum / (t.nProbes[0] * t.nProbes[1])
		Assert(sx >= 0 && sx < t.nProbes[0])
		Assert(sy >= 0 && sy < t.nProbes[1])
		Assert(sz >= 0 && sz < t.nProbes[2])
		tx0, tx1 := float64(sx)/float64(t.nProbes[0]), float64(sx+1)/float64(t.nProbes[0])
		ty0, ty1 := float64(sy)/float64(t.nProbes[1]), float64(sy+1)/float64(t.nProbes[1])
		tz0, tz1 := float64(sz)/float64(t.nProbes[2]), float64(sz+1)/float64(t.nProbes[2])
		b := CreateBBoxFromPoints(t.bbox.Lerp(tx0, ty0, tz0), t.bbox.Lerp(tx1, ty1, tz1))

		// Initialize common variables for _CreateRadProbeTask::Run()_
		rng := NewRNG(int64(t.pointNum))
		c_probe := make([]Spectrum, SHTerms(t.lmax), SHTerms(t.lmax))
		var arena *MemoryArena
		nFound := 0
		lastVisibleOffset := 0
		for i := 0; i < 256; i++ {
			if nFound == 32 {
				break
			}
			// Try to compute radiance probe contribution at _i_th sample point

			// Compute _i_th candidate point _p_ in cell's bounding box
			dx := RadicalInverse(i+1, 2)
			dy := RadicalInverse(i+1, 3)
			dz := RadicalInverse(i+1, 5)
			p := b.Lerp(dx, dy, dz)

			// Skip point _p_ if not indirectly visible from camera
			if t.scene.IntersectP(CreateRay(&t.surfacePoints[lastVisibleOffset], p.Sub(&t.surfacePoints[lastVisibleOffset]), 1.0e-4, 1.0, t.time, 0)) {
				// See if point is visible to any element of _surfacePoints_
				var j int
				for j = 0; j < len(t.surfacePoints); j++ {
					if !t.scene.IntersectP(CreateRay(&t.surfacePoints[j], p.Sub(&t.surfacePoints[j]), 1.0e-4, 1.0, t.time, 0)) {
						lastVisibleOffset = j
						break
					}
				}
				if j == len(t.surfacePoints) {
					continue
				}
			}
			nFound++

			// Compute SH coefficients of incident radiance at point _p_
			if t.includeDirectInProbes {
				for i := 0; i < SHTerms(t.lmax); i++ {
					c_probe[i] = *NewSpectrum1(0.0)
				}
				SHProjectIncidentDirectRadiance(p, 0.0, t.time, arena, t.scene, true, t.lmax, rng, c_probe)
				for i := 0; i < SHTerms(t.lmax); i++ {
					t.c_in[i] = *t.c_in[i].Add(&c_probe[i])
				}
			}

			if t.includeIndirectInProbes {
				for i := 0; i < SHTerms(t.lmax); i++ {
					c_probe[i] = *NewSpectrum1(0.0)
				}
				SHProjectIncidentIndirectRadiance(p, 0.0, t.time, t.renderer, t.origSample, t.scene, t.lmax, rng, t.nIndirSamples, c_probe)
				for i := 0; i < SHTerms(t.lmax); i++ {
					t.c_in[i] = *t.c_in[i].Add(&c_probe[i])
				}
			}
		}
		// Compute final average value for probe and cleanup
		if nFound > 0 {
			for i := 0; i < SHTerms(t.lmax); i++ {
				t.c_in[i] = *t.c_in[i].InvScale(float64(nFound))
			}
		}

		t.progress.Update(1)
	}

	completed <- true
}
