package pbrt

type Shape interface {
    ObjectBound() *BBox
    WorldBound() *BBox
    CanIntersect() bool
    Refine() (refined []*Shape)
    Intersect(ray *Ray) (bool, tHit, rayEpsilon float64, dg *DifferentialGeometry)
    IntersectP(ray *Ray) bool
    GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry
    Area() float64
    Sample(u1, u2 float64) (*Point, *Normal)
    Pdf(pshape *Point) float64
    SampleAt(p *Point, u1, u2 float64) (*Point, *Normal)
    Pdf2(p *Point, wi *Vector) float64	// TODO: create better name for this function
    ReverseOrientation() bool
    TransformSwapsHandedness() bool
}