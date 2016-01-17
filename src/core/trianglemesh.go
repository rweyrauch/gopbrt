package core

import ()

type TriangleMesh struct {
	ShapeData
}

func CreateTriangleMeshShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet, textures map[string]TextureFloat) *TriangleMesh {
	return nil
}

func (c *TriangleMesh) ObjectBound() *BBox {
	return nil
}

func (c *TriangleMesh) WorldBound() *BBox {
	return nil
}
func (c *TriangleMesh) CanIntersect() bool {
	return false
}
func (c *TriangleMesh) Refine(refined []Shape) []Shape {
	return refined
}
func (c *TriangleMesh) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *TriangleMesh) IntersectP(ray *Ray) bool {
	return false
}
func (c *TriangleMesh) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *TriangleMesh) Area() float64 {
	return 0.0
}
func (c *TriangleMesh) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *TriangleMesh) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *TriangleMesh) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *TriangleMesh) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *TriangleMesh) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *TriangleMesh) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *TriangleMesh) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *TriangleMesh) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *TriangleMesh) ShapeId() uint32 {
	return c.shapeId
}
