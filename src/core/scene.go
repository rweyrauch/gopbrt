package pbrt

type Scene struct {
	aggregate    Primitive
	lights       []Light
	volumeRegion VolumeRegion
	bound        BBox
}

func (*Scene) Intersect(ray *Ray, isect *Intersection) bool {
	return false
}

func (*Scene) IntersectP(ray *Ray) bool {
	return false
}
