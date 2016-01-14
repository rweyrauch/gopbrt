package pbrt

import ()

type Cone struct {
	ShapeData
}

func CreateConeShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cone {
	return nil
}

func (c *Cone) ObjectBound() *BBox {
	return nil
}

func (c *Cone) WorldBound() *BBox {
	return nil
}
func (c *Cone) CanIntersect() bool {
	return false
}
func (c *Cone) Refine() (refined []*Shape) {
	return nil
}
func (c *Cone) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Cone) IntersectP(ray *Ray) bool {
	return false
}
func (c *Cone) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Cone) Area() float64 {
	return 0.0
}
func (c *Cone) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Cone) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Cone) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Cone) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Cone) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Cone) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Cone) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Cone) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Cone) ShapeId() uint32 {
	return c.shapeId
}
