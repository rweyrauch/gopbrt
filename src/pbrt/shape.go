package pbrt

import (
	"math"
)

var (
	nextShapeId uint32 = 0
)

func GenerateShapeId() uint32 {
	nextShapeId++
	return nextShapeId
}

type Shape interface {
	ObjectBound() *BBox
	WorldBound() *BBox
	CanIntersect() bool
	Refine() (refined []*Shape)
	Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry)
	IntersectP(ray *Ray) bool
	GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry
	Area() float64
	Sample(u1, u2 float64) (*Point, *Normal)
	Pdf(pshape *Point) float64
	SampleAt(p *Point, u1, u2 float64) (*Point, *Normal)
	Pdf2(p *Point, wi *Vector) float64 // TODO: create better name for this function
	ObjectToWorld() *Transform
	WorldToObject() *Transform
	ReverseOrientation() bool
	TransformSwapsHandedness() bool
	ShapeId() uint32
}

type ShapeData struct {
	objectToWorld, worldToObject                 *Transform
	reverseOrientation, transformSwapsHandedness bool
	shapeId                                      uint32
}

func ShapePdf(shape Shape, p *Point, wi *Vector) float64 {
	// Intersect sample ray with area light geometry
	ray := CreateRay(p, wi, 1.0e-3, INFINITY, 0.0, 0)
	ray.depth = -1 // temporary hack to ignore alpha mask

	var thit float64
	var dgLight *DifferentialGeometry
	var ok bool
	if ok, thit, _, dgLight = shape.Intersect(ray); !ok {
		return 0.0
	}

	// Convert light sample weight to solid angle measure
	pdf := DistanceSquaredPoint(p, ray.PointAt(thit)) / (AbsDotNormalVector(dgLight.nn, wi.Negate()) * shape.Area())
	if math.IsInf(pdf, 0) {
		pdf = 0.0
	}
	return pdf
}
