package core

type AggregateTest struct {
	nIterations int
	primitives []Primitive
	bboxes []BBox
}

func NewAggregateTest(niterations int, primitives []Primitive) *AggregateTest {
	aggtest := &AggregateTest{niterations, primitives, nil}
	
    for _, p := range aggtest.primitives {
        aggtest.primitives = p.FullyRefine(aggtest.primitives)
    }
    aggtest.bboxes = make([]BBox, len(aggtest.primitives),len(aggtest.primitives)) 
    for i, p := range aggtest.primitives {
        aggtest.bboxes[i] = *p.WorldBound()
	}
	return aggtest
}

func (aggtest *AggregateTest) Render(scene *Scene) {
    rng := NewRNG(42)
    	prog := NewProgressReporter(aggtest.nIterations, "Aggregate Test", -1)
    	
    // Compute bounding box of region used to generate random rays
    bbox := scene.WorldBound()
    bbox.Expand(bbox.pMax.At(bbox.MaximumExtent()) - bbox.pMin.At(bbox.MaximumExtent()))
    var lastHit Point
    lastEps := 0.0
    for i := 0; i < aggtest.nIterations; i++ {
        // Choose random rays, _rayAccel_ and _rayAll_ for testing

        // Choose ray origin for testing accelerator
        org := *CreatePoint(Lerp(rng.RandomFloat(), bbox.pMin.x, bbox.pMax.x),
                  Lerp(rng.RandomFloat(), bbox.pMin.y, bbox.pMax.y),
                  Lerp(rng.RandomFloat(), bbox.pMin.z, bbox.pMax.z))
        if (rng.RandomUInt() % 4) == 0 { org = lastHit }

        // Choose ray direction for testing accelerator
        dir := *UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat())
        if (rng.RandomUInt() % 32) == 0 { 
        	dir.x, dir.y = 0.0, 0.0 
       	} else if (rng.RandomUInt() % 32) == 0 {
       		dir.x, dir.z = 0.0, 0.0
        }  else if (rng.RandomUInt() % 32) == 0 {
          	dir.y, dir.z = 0.0, 0.0
		}
        // Choose ray epsilon for testing accelerator
        eps := 0.0
        if rng.RandomFloat() < 0.25 { 
        	eps = lastEps 
       	} else if rng.RandomFloat() < 0.25 {
       		eps = 1.0e-3
		}       		
        rayAccel := CreateRayDifferential(&org, &dir, eps, INFINITY, 0.0, 1)
        rayAll := CreateRayFromRayDifferential(rayAccel)
		rayDiffAll := CreateRayDifferential(&org, &dir, eps, INFINITY, 0.0, 1)
		
        // Compute intersections using accelerator and exhaustive testing
        var isectAll Intersection 
        hitAccel, _ := scene.Intersect(rayAccel)
        hitAll := false
        inconsistentBounds := false
        for j, p := range aggtest.primitives {
        	hitP, _, _ := aggtest.bboxes[j].IntersectP(rayAll)
            if hitP {
                hit2, isect := p.Intersect(rayDiffAll)
                hitAll = hitAll || hit2
                if hit2 { isectAll = *isect }
            } else {
            	hit, _ := p.Intersect(rayDiffAll)
            	if hit { inconsistentBounds = true }
            }    
        }

        // Report any inconsistencies between intersections
        if !inconsistentBounds &&
            ((hitAccel != hitAll) || (rayAccel.maxt != rayDiffAll.maxt)) {
            Warning("Disagreement: t accel %.16g t exhaustive %.16g\nRay: org [%g, %g, %g], dir [%g, %g, %g], mint = %g",
                    rayAccel.maxt, rayDiffAll.maxt,
                    rayDiffAll.origin.x, rayDiffAll.origin.y, rayDiffAll.origin.z,
                    rayDiffAll.dir.x, rayDiffAll.dir.y, rayDiffAll.dir.z, rayDiffAll.mint)
		}            
        if hitAll {
            lastHit = *rayAll.PointAt(rayAll.maxt)
            lastEps = isectAll.rayEpsilon
        }
        prog.Update(1)
    }
    prog.Done()
	
}

func (r *AggregateTest) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) (*Spectrum, *Intersection, *Spectrum) {
	return nil, nil, nil
}

func (r *AggregateTest) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum {
	return nil
}

func CreateAggregateTestRenderer(param *ParamSet, primitives []Primitive) *AggregateTest { 
    niters := param.FindIntParam("niters", 100000)
    return NewAggregateTest(niters, primitives)
	}
