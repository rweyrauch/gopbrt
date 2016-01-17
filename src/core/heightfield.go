package pbrt

import ()

type Heightfield struct {
	ShapeData
}

func CreateHeightfieldShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Heightfield {
	return nil
}

func (c *Heightfield) ObjectBound() *BBox {
	return nil
}

func (c *Heightfield) WorldBound() *BBox {
	return nil
}
func (c *Heightfield) CanIntersect() bool {
	return false
}
func (c *Heightfield) Refine(refined []Shape) []Shape {
	return refined
}
func (c *Heightfield) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Heightfield) IntersectP(ray *Ray) bool {
	return false
}
func (c *Heightfield) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Heightfield) Area() float64 {
	return 0.0
}
func (c *Heightfield) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Heightfield) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Heightfield) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Heightfield) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Heightfield) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Heightfield) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Heightfield) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Heightfield) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Heightfield) ShapeId() uint32 {
	return c.shapeId
}
