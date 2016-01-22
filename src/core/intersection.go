package core

type Intersection struct {
	dg                           *DifferentialGeometry
	primitive                    Primitive
	WorldToObject, ObjectToWorld *Transform
	shapeId, primitiveId         uint32
	rayEpsilon                   float64
}

func (isect *Intersection) GetBSDF(ray *RayDifferential, arena *MemoryArena) *BSDF {
    //PBRT_STARTED_BSDF_SHADING(const_cast<RayDifferential *>(&ray));
    isect.dg.ComputeDifferentials(ray)
    bsdf := isect.primitive.GetBSDF(isect.dg, isect.ObjectToWorld, arena)
    //PBRT_FINISHED_BSDF_SHADING(const_cast<RayDifferential *>(&ray), bsdf);
    return bsdf;
}
func (isect *Intersection) GetBSSRDF(ray *RayDifferential, arena *MemoryArena) *BSSRDF {
    //PBRT_STARTED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray));
    isect.dg.ComputeDifferentials(ray)
    bssrdf := isect.primitive.GetBSSRDF(isect.dg, isect.ObjectToWorld, arena)
    //PBRT_FINISHED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray), bssrdf);
    return bssrdf;
}
func (isect *Intersection) Le(wo *Vector) *Spectrum {
    area := isect.primitive.GetAreaLight()
    if area != nil {
		return area.L(isect.dg.p, isect.dg.nn, wo)
	} else {
		return NewSpectrum1(0.0)
	}
}
