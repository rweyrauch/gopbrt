package core

import (
	"os"
	"fmt"
	"math"
)

type (
    
   SurfacePoint struct {
       p Point
       n Normal
       area, rayEpsilon float64
   }
   
   SurfacePointsRenderer struct {
       minDist, time float64
       pCamera Point
       filename string
       points []SurfacePoint
   }
	
	surfacePointTask struct {
		taskNum int
		scene *Scene
		origin Point
		time float64
		minSampleDist float64
		maxFails int
		
    	repeatedFails, maxRepeatedFails int
    	totalPathsTraced, totalRaysTraced, numPointsAdded int
		
   	 	sphere *GeometricPrimitive
    	octree *Octree
    	surfacePoints []SurfacePoint
		prog *ProgressReporter
	}
	
	poissonCheck struct {
		maxDist2 float64
		failed bool
		p Point
	}
)

func NewSurfacePointsRenderer(minDist float64, pCamera *Point, time float64, filename string) *SurfacePointsRenderer {
	renderer := new(SurfacePointsRenderer)
	
	return renderer
}

func (renderer *SurfacePointsRenderer) Render(scene *Scene) {
    // Declare shared variables for Poisson point generation
    octBounds := scene.WorldBound()
    octBounds.Expand(.0010 * math.Pow(octBounds.Volume(), 1.0/3.0))
    pointOctree := NewOctree(octBounds, 16)

    // Create scene bounding sphere to catch rays that leave the scene
    sceneCenter, sceneRadius := scene.WorldBound().BoundingSphere()
    ObjectToWorld := TranslateTransform(sceneCenter.Sub(CreatePoint(0,0,0)))
    WorldToObject := InverseTransform(ObjectToWorld)
    sph := CreateSphere(ObjectToWorld, WorldToObject,
        true, sceneRadius, -sceneRadius, sceneRadius, 360.0)
    sphere := NewGeometricPrimitive(sph, nil, nil)
    maxFails, repeatedFails, maxRepeatedFails := 2000, 0, 0
    if options.QuickRender { maxFails = Maxi(10, maxFails / 10) }
    totalPathsTraced, totalRaysTraced, numPointsAdded := 0, 0, 0
    prog := NewProgressReporter(maxFails, "Depositing samples", TerminalWidth())
    // Launch tasks to trace rays to find Poisson points
    //PBRT_SUBSURFACE_STARTED_RAYS_FOR_POINTS();
    //vector<Task *> tasks;
    //RWMutex *mutex = RWMutex::Create();
    nTasks := NumSystemCores()
    for i := 0; i < nTasks; i++ {
        spt := newSurfacePointTask(scene, &renderer.pCamera, renderer.time, i,
            renderer.minDist, maxFails, repeatedFails, maxRepeatedFails,
            totalPathsTraced, totalRaysTraced, numPointsAdded, sphere, pointOctree,
            renderer.points, prog)
        spt.Run()
    }   
    //EnqueueTasks(tasks);
    //WaitForAllTasks();
    //for (uint32_t i = 0; i < tasks.size(); ++i)
    //    delete tasks[i];
    //RWMutex::Destroy(mutex);
    prog.Done()
    //PBRT_SUBSURFACE_FINISHED_RAYS_FOR_POINTS(totalRaysTraced, numPointsAdded);
    if len(renderer.filename) != 0 {
        // Write surface points to file
        fi, err := os.Open(renderer.filename)
        defer fi.Close()
        if err != nil {
            Error("Unable to open output file \"%s\" (%s)", renderer.filename, err)
            return
        }

        fmt.Fprintf(fi, "# points generated by SurfacePointsRenderer\n")
        fmt.Fprintf(fi, "# position (x,y,z), normal (x,y,z), area, rayEpsilon\n")
        for _, sp := range renderer.points {
            fmt.Fprintf(fi, "%g %g %g %g %g %g %g %g\n", sp.p.x, sp.p.y, sp.p.z,
                sp.n.x, sp.n.y, sp.n.z, sp.area, sp.rayEpsilon)
        }
    }	
}

func (*SurfacePointsRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (li *Spectrum, isect *Intersection, T *Spectrum) { return NewSpectrum1(0.0), nil, nil }
func (*SurfacePointsRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return NewSpectrum1(0.0) }

func newSurfacePointTask(scene *Scene, pCamera *Point, time float64, taskNum int,
            minDist float64, maxFails, repeatedFails, maxRepeatedFails,
            totalPathsTraced, totalRaysTraced, numPointsAdded int, sphere Primitive, octree *Octree, points []SurfacePoint, prog *ProgressReporter) *surfacePointTask {
 	return nil          	
}
            
func (spt *surfacePointTask) Run() {
    // Declare common variables for _SurfacePointTask::Run()_
    rng := NewRNG(int64(37 * spt.taskNum))
    var arena *MemoryArena
    candidates := make([]SurfacePoint, 0, 100)
    for {
        pathsTraced, raysTraced := 0, 0
        for pathsTraced = 0; pathsTraced < 20000; pathsTraced++ {
            // Follow ray path and attempt to deposit candidate sample points
            dir := UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat())
            ray := CreateRay(&spt.origin, dir, 0.0, INFINITY, spt.time, 0)
            for ray.depth < 30 {
                // Find ray intersection with scene geometry or bounding sphere
                raysTraced++
                
                var isect *Intersection
                var hit bool
                
                hitOnSphere := false
                if hit, isect = spt.scene.Intersect(CreateRayDifferentialFromRay(ray)); !hit {
                    if hit, isect = spt.sphere.Intersect(CreateRayDifferentialFromRay(ray)); !hit {
                        break
                    }    
                    hitOnSphere = true
                }
                hitGeometry := isect.dg
                hitGeometry.nn = FaceforwardNormalVector(hitGeometry.nn, ray.dir.Negate())

                // Store candidate sample point at ray intersection if appropriate
                if !hitOnSphere && ray.depth >= 3 &&
                    isect.GetBSSRDF(CreateRayDifferentialFromRay(ray), arena) != nil {
                    area := math.Pi * (spt.minSampleDist / 2.0) * (spt.minSampleDist / 2.0)
                    candidates = append(candidates, SurfacePoint{*hitGeometry.p, *hitGeometry.nn,
                                                      area, isect.rayEpsilon})
                }

                // Generate random ray from intersection point
                dir := UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat())
                dir = FaceforwardVectorNormal(dir, hitGeometry.nn)
                ray = CreateChildRay(hitGeometry.p, dir, ray, isect.rayEpsilon, INFINITY)
            }
            //arena.FreeAll()
        }
        // Make first pass through candidate points with reader lock
        candidateRejected := make([]bool, 0, len(candidates))
        //RWMutexLock lock(mutex, READ);
                
        for _, c := range candidates {
            check := poissonCheck{spt.minSampleDist, false, c.p}
	        checkPoisson := func(obj Object) bool {
	        	sp, ok := obj.(*SurfacePoint)
	        	if ok {
	        		if DistanceSquaredPoint(&sp.p, &check.p) < check.maxDist2 {
	        			check.failed = true
	        			return false
	        		}
	        	}
	        	return true
	        }
            
            spt.octree.Lookup(&c.p, checkPoisson)
            candidateRejected = append(candidateRejected, check.failed)
        }

        // Make second pass through points with writer lock and update octree
        //lock.UpgradeToWrite();
        if spt.repeatedFails >= spt.maxFails {
            return
		}            
        spt.totalPathsTraced += pathsTraced
        spt.totalRaysTraced += raysTraced 
        oldMaxRepeatedFails := spt.maxRepeatedFails;
        for i, sp := range candidates {
            if candidateRejected[i] {
                // Update for rejected candidate point
                spt.repeatedFails++
                spt.maxRepeatedFails = Maxi(spt.maxRepeatedFails, spt.repeatedFails)
                if spt.repeatedFails >= spt.maxFails {
                    return
                }
            } else {
                // Recheck candidate point and possibly add to octree
                check := poissonCheck{spt.minSampleDist, false, sp.p}
                
                checkPoisson := func(obj Object) bool {
                	sp, ok := obj.(*SurfacePoint)
                	if ok {
                		if DistanceSquaredPoint(&sp.p, &check.p) < check.maxDist2 {
                			check.failed = true
                			return false
                		}
                	}
                	return true
                }
                
                spt.octree.Lookup(&sp.p, checkPoisson)
                if check.failed {
                    // Update for rejected candidate point
                    spt.repeatedFails++
                    spt.maxRepeatedFails = Maxi(spt.maxRepeatedFails, spt.repeatedFails)
                    if spt.repeatedFails >= spt.maxFails {
                        return
                    }    
                } else {
                    spt.numPointsAdded++
                    spt.repeatedFails = 0
                    delta := CreateVector(spt.minSampleDist, spt.minSampleDist, spt.minSampleDist)
                    spt.octree.Add(sp, CreateBBoxFromPoints(sp.p.SubVector(delta), sp.p.Add(delta)))
                    //PBRT_SUBSURFACE_ADDED_POINT_TO_OCTREE(&sp, minSampleDist);
                    spt.surfacePoints = append(spt.surfacePoints, sp)
                }
            }
        }

        // Stop following paths if not finding new points
        if spt.repeatedFails > oldMaxRepeatedFails {
            delta := spt.repeatedFails - oldMaxRepeatedFails
            spt.prog.Update(delta)
        }
        if spt.totalPathsTraced > 50000 && spt.numPointsAdded == 0 {
            Warning("There don't seem to be any objects with BSSRDFs in this scene.  Giving up.");
            return
        }
        candidates = candidates[0:0]
    }
}

func CreateSurfacePointsRenderer(params *ParamSet, pCamera *Point, time float64) *SurfacePointsRenderer { 
     minDist := params.FindFloatParam("minsampledistance", 0.25)
     filename := params.FindFilenameParam("filename", "")
    if options.QuickRender { minDist *= 4.0 }
    return NewSurfacePointsRenderer(minDist, pCamera, time, filename)
}

func FindPoissonPointDistribution(pCamera *Point, time, minDist float64, scene *Scene, points *[]SurfacePoint) {
    sp := NewSurfacePointsRenderer(minDist, pCamera, time, "")
    sp.Render(scene)
    *points = make([]SurfacePoint, len(sp.points), len(sp.points))
    copy(*points, sp.points)
}