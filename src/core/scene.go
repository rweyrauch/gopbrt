package core

type Scene struct {
	aggregate    Primitive
	lights       []Light
	volumeRegion VolumeRegion
	bound        BBox
}

func CreateScene(accel Primitive, lights []Light, vr VolumeRegion) *Scene {
    scene := new(Scene)
    scene.aggregate = accel
    scene.lights = lights
    scene.volumeRegion = vr
    scene.bound = *accel.WorldBound()
    if scene.volumeRegion != nil {
        scene.bound = *UnionBBoxes(&scene.bound, scene.volumeRegion.WorldBound())
    }
    return scene
}
func (s *Scene) Intersect(ray *Ray) (hit bool,  isect *Intersection) {
	return s.aggregate.Intersect(ray)
}

func (s *Scene) IntersectP(ray *Ray) bool {
	return s.aggregate.IntersectP(ray) 
}
