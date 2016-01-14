package pbrt

import ()

type Hyperboloid struct {
	ShapeData
}

func CreateHyperboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Hyperboloid {
	return nil
}

func (c *Hyperboloid) ObjectBound() *BBox {
	return nil
}

func (c *Hyperboloid) WorldBound() *BBox {
	return nil
}
func (c *Hyperboloid) CanIntersect() bool {
	return false
}
func (c *Hyperboloid) Refine() (refined []*Shape) {
	return nil
}
func (c *Hyperboloid) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Hyperboloid) IntersectP(ray *Ray) bool {
	return false
}
func (c *Hyperboloid) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Hyperboloid) Area() float64 {
	return 0.0
}
func (c *Hyperboloid) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Hyperboloid) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Hyperboloid) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Hyperboloid) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Hyperboloid) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Hyperboloid) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Hyperboloid) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Hyperboloid) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Hyperboloid) ShapeId() uint32 {
	return c.shapeId
}
