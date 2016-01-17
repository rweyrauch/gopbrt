package core

type Cylinder struct {
	ShapeData
}

func CreateCylinderShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cylinder {
	return nil
}

func (c *Cylinder) ObjectBound() *BBox {
	return nil
}

func (c *Cylinder) WorldBound() *BBox {
	return nil
}
func (c *Cylinder) CanIntersect() bool {
	return false
}
func (c *Cylinder) Refine(refined []Shape) []Shape {
	return refined
}
func (c *Cylinder) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Cylinder) IntersectP(ray *Ray) bool {
	return false
}
func (c *Cylinder) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Cylinder) Area() float64 {
	return 0.0
}
func (c *Cylinder) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Cylinder) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Cylinder) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Cylinder) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Cylinder) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Cylinder) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Cylinder) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Cylinder) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Cylinder) ShapeId() uint32 {
	return c.shapeId
}
