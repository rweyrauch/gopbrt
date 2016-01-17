package core

import ()

type Paraboloid struct {
	ShapeData
}

func CreateParaboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Paraboloid {
	return nil
}

func (c *Paraboloid) ObjectBound() *BBox {
	return nil
}

func (c *Paraboloid) WorldBound() *BBox {
	return nil
}
func (c *Paraboloid) CanIntersect() bool {
	return false
}
func (c *Paraboloid) Refine(refined []Shape) []Shape {
	return refined
}
func (c *Paraboloid) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Paraboloid) IntersectP(ray *Ray) bool {
	return false
}
func (c *Paraboloid) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Paraboloid) Area() float64 {
	return 0.0
}
func (c *Paraboloid) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Paraboloid) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Paraboloid) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Paraboloid) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Paraboloid) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Paraboloid) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Paraboloid) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Paraboloid) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Paraboloid) ShapeId() uint32 {
	return c.shapeId
}
