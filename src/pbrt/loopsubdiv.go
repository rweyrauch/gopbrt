package pbrt

import ()

type LoopSubdiv struct {
	ShapeData
}

func CreateLoopSubdivShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *LoopSubdiv {
	return nil
}

func (c *LoopSubdiv) ObjectBound() *BBox {
	return nil
}

func (c *LoopSubdiv) WorldBound() *BBox {
	return nil
}
func (c *LoopSubdiv) CanIntersect() bool {
	return false
}
func (c *LoopSubdiv) Refine(refined []Shape) []Shape {
	return refined
}
func (c *LoopSubdiv) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *LoopSubdiv) IntersectP(ray *Ray) bool {
	return false
}
func (c *LoopSubdiv) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *LoopSubdiv) Area() float64 {
	return 0.0
}
func (c *LoopSubdiv) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *LoopSubdiv) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *LoopSubdiv) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *LoopSubdiv) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *LoopSubdiv) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *LoopSubdiv) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *LoopSubdiv) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *LoopSubdiv) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *LoopSubdiv) ShapeId() uint32 {
	return c.shapeId
}
