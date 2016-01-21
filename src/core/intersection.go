package core

type Intersection struct {
	dg                           *DifferentialGeometry
	primitive                    Primitive
	WorldToObject, ObjectToWorld *Transform
	shapeId, primitiveId         uint32
	rayEpsilon                   float64
}

func (isect *Intersection) GetBSDF(ray *RayDifferential, arena *MemoryArena) *BSDF {
	// TODO: implement this
	return nil
}
func (isect *Intersection) GetBSSRDF(ray *RayDifferential, arena *MemoryArena) *BSSRDF {
	// TODO: implement this
	return nil
}
func (isect *Intersection) Le(wo *Vector) *Spectrum {
	// TODO: implement this
	return nil
}
