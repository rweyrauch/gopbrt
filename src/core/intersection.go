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

type Intersection struct {
	dg                           *DifferentialGeometry
	primitive                    Primitive
	WorldToObject, ObjectToWorld *Transform
	shapeId, primitiveId         uint32
	rayEpsilon                   float64
}

func NewIntersection() *Intersection {
	isect := new(Intersection)
	isect.dg = CreateDiffGeometry(CreatePoint(0, 0, 0), CreateVector(0, 0, 0), CreateVector(0, 0, 0), CreateNormal(0, 0, 0), CreateNormal(0, 0, 0), 0.0, 0.0, nil)
	isect.primitive = nil
	isect.WorldToObject = NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4())
	isect.ObjectToWorld = InverseTransform(isect.WorldToObject)
	isect.shapeId = 0
	isect.primitiveId = 0
	isect.rayEpsilon = 0
	return isect
}
func (isect *Intersection) GetBSDF(ray *RayDifferential, arena *MemoryArena) *BSDF {
	isect.dg.ComputeDifferentials(ray)
	bsdf := isect.primitive.GetBSDF(isect.dg, isect.ObjectToWorld, arena)
	return bsdf
}
func (isect *Intersection) GetBSSRDF(ray *RayDifferential, arena *MemoryArena) *BSSRDF {
	isect.dg.ComputeDifferentials(ray)
	bssrdf := isect.primitive.GetBSSRDF(isect.dg, isect.ObjectToWorld, arena)
	return bssrdf
}
func (isect *Intersection) Le(wo *Vector) *Spectrum {
	area := isect.primitive.GetAreaLight()
	if area != nil {
		return area.L(isect.dg.p, isect.dg.nn, wo)
	} else {
		return NewSpectrum1(0.0)
	}
}
