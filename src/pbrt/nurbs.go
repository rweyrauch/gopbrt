package pbrt

import ()

type NURBS struct {
	ShapeData
}

func CreateNURBSShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *NURBS {
	return nil
}

func (c *NURBS) ObjectBound() *BBox {
	return nil
}

func (c *NURBS) WorldBound() *BBox {
	return nil
}
func (c *NURBS) CanIntersect() bool {
	return false
}
func (c *NURBS) Refine() (refined []Shape) {
	return nil
}
func (c *NURBS) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *NURBS) IntersectP(ray *Ray) bool {
	return false
}
func (c *NURBS) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *NURBS) Area() float64 {
	return 0.0
}
func (c *NURBS) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *NURBS) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *NURBS) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *NURBS) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *NURBS) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *NURBS) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *NURBS) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *NURBS) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *NURBS) ShapeId() uint32 {
	return c.shapeId
}
