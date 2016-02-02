/*
	gopbrt

	Port of pbrt v2.0.0 by Matt Pharr and Greg Humphreys to the go language.
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

	The MIT License (MIT)
	Copyright (c) 2016 Rick Weyrauch

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
package core

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
	Refine(refined []Shape) []Shape
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
