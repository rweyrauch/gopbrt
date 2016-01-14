package pbrt

import ()

type Disk struct {
	ShapeData
}

func CreateDiskShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Disk {
	return nil
}

func (c *Disk) ObjectBound() *BBox {
	return nil
}

func (c *Disk) WorldBound() *BBox {
	return nil
}
func (c *Disk) CanIntersect() bool {
	return false
}
func (c *Disk) Refine() (refined []*Shape) {
	return nil
}
func (c *Disk) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *Disk) IntersectP(ray *Ray) bool {
	return false
}
func (c *Disk) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *Disk) Area() float64 {
	return 0.0
}
func (c *Disk) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Disk) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *Disk) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *Disk) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *Disk) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *Disk) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *Disk) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *Disk) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *Disk) ShapeId() uint32 {
	return c.shapeId
}
