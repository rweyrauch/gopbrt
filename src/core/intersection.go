package pbrt

type Intersection struct {
	dg                           *DifferentialGeometry
	primitive                    *Primitive
	WorldToObject, ObjectToWorld *Transform
	shapeId, primitiveId         uint32
	rayEpsilon                   float64
}

func GetBSDF(ray *RayDifferential, arena *MemoryArena) *BSDF {
	// TODO: implement this
	return nil
}
func GetBSSRDF(ray *RayDifferential, arena *MemoryArena) *BSSRDF {
	// TODO: implement this
	return nil
}
func Le(wo *Vector) *Spectrum {
	// TODO: implement this
	return nil
}
