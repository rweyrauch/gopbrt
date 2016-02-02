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

import ()

type LoopSubdiv struct {
	ShapeData
}

func CreateLoopSubdivShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *LoopSubdiv {
	return nil
}

func (c *LoopSubdiv) ObjectBound() *BBox {
	return nil
}

func (c *LoopSubdiv) WorldBound() *BBox {
	return nil
}
func (c *LoopSubdiv) CanIntersect() bool {
	return false
}
func (c *LoopSubdiv) Refine(refined []Shape) []Shape {
	return refined
}
func (c *LoopSubdiv) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (c *LoopSubdiv) IntersectP(ray *Ray) bool {
	return false
}
func (c *LoopSubdiv) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (c *LoopSubdiv) Area() float64 {
	return 0.0
}
func (c *LoopSubdiv) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *LoopSubdiv) Pdf(pshape *Point) float64 {
	return 0.0
}
func (c *LoopSubdiv) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (c *LoopSubdiv) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (c *LoopSubdiv) ObjectToWorld() *Transform {
	return c.objectToWorld
}
func (c *LoopSubdiv) WorldToObject() *Transform {
	return c.worldToObject
}
func (c *LoopSubdiv) ReverseOrientation() bool {
	return c.reverseOrientation
}
func (c *LoopSubdiv) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}
func (c *LoopSubdiv) ShapeId() uint32 {
	return c.shapeId
}
