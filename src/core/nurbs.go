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

type NURBS struct {
	ShapeData
	nu, uorder, nv, vorder int
	umin, umax, vmin, vmax float64
	uknot, vknot           []float64
	isHomogeneous          bool
	P                      []float64
}

type Homogeneous3 struct {
	x, y, z, w float64
}

func NewNURBS(o2w, w2o *Transform, ro bool, nu, uorder int, uknot []float64, umin, umax float64,
	nv, vorder int, vknot []float64, vmin, vmax float64, P []float64, isHomogeneous bool) *NURBS {
	nurbs := new(NURBS)
	nurbs.objectToWorld = o2w
	nurbs.worldToObject = w2o
	nurbs.reverseOrientation = ro
	nurbs.transformSwapsHandedness = SwapsHandednessTransform(nurbs.objectToWorld)
	nurbs.shapeId = GenerateShapeId()

	nurbs.nu = nu
	nurbs.uorder = uorder
	nurbs.umin = umin
	nurbs.umax = umax

	nurbs.nv = nv
	nurbs.vorder = vorder
	nurbs.vmin = vmin
	nurbs.vmax = vmax
	nurbs.isHomogeneous = isHomogeneous
	if nurbs.isHomogeneous {
		Assert(len(P) == 4*nu*nv)
	} else {
		Assert(len(P) == 3*nu*nv)
	}
	nurbs.P = make([]float64, len(P), len(P))
	copy(nurbs.P, P)
	Assert(len(uknot) == nurbs.nu+nurbs.uorder)
	nurbs.uknot = make([]float64, len(uknot), len(uknot))
	copy(nurbs.uknot, uknot)
	Assert(len(vknot) == nurbs.nv+nurbs.vorder)
	nurbs.vknot = make([]float64, len(vknot), len(vknot))
	copy(nurbs.vknot, vknot)

	return nurbs
}

func (nurbs *NURBS) ObjectBound() *BBox {
	if !nurbs.isHomogeneous {
		// Compute object-space bound of non-homogeneous NURBS
		bound := CreateEmptyBBox()
		pi := 0
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			bound = UnionBBoxPoint(bound, CreatePoint(nurbs.P[pi+0], nurbs.P[pi+1], nurbs.P[pi+2]))
			pi += 3
		}
		return bound
	} else {
		// Compute object-space bound of homogeneous NURBS
		bound := CreateEmptyBBox()
		pi := 0
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			bound = UnionBBoxPoint(bound, CreatePoint(nurbs.P[pi+0]/nurbs.P[pi+3], nurbs.P[pi+1]/nurbs.P[pi+3], nurbs.P[pi+2]/nurbs.P[pi+3]))
			pi += 4
		}
		return bound
	}
}

func (nurbs *NURBS) WorldBound() *BBox {
	if !nurbs.isHomogeneous {
		// Compute world-space bound of non-homogeneous NURBS
		bound := CreateEmptyBBox()
		pi := 0
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			pt := PointTransform(nurbs.objectToWorld, CreatePoint(nurbs.P[pi+0], nurbs.P[pi+1], nurbs.P[pi+2]))
			bound = UnionBBoxPoint(bound, pt)
			pi += 3
		}
		return bound
	} else {
		// Compute world-space bound of homogeneous NURBS
		bound := CreateEmptyBBox()
		pi := 0
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			pt := PointTransform(nurbs.objectToWorld, CreatePoint(nurbs.P[pi+0]/nurbs.P[pi+3], nurbs.P[pi+1]/nurbs.P[pi+3], nurbs.P[pi+2]/nurbs.P[pi+3]))
			bound = UnionBBoxPoint(bound, pt)
			pi += 4
		}
		return bound
	}
}
func (nurbs *NURBS) CanIntersect() bool {
	return false
}
func (nurbs *NURBS) Refine(refined []Shape) []Shape {
	// Compute NURBS dicing rates
	diceu, dicev := 30, 30
	ueval := make([]float64, diceu, diceu)
	veval := make([]float64, dicev, dicev)
	evalPs := make([]Point, diceu*dicev, diceu*dicev)
	evalNs := make([]Normal, diceu*dicev, diceu*dicev)

	for i := 0; i < diceu; i++ {
		ueval[i] = Lerp(float64(i)/float64(diceu-1), nurbs.umin, nurbs.umax)
	}
	for i := 0; i < dicev; i++ {
		veval[i] = Lerp(float64(i)/float64(dicev-1), nurbs.vmin, nurbs.vmax)
	}
	// Evaluate NURBS over grid of points
	uvs := make([]float64, 2*diceu*dicev, 2*diceu*dicev)
	// Turn NURBS into triangles
	Pw := make([]Homogeneous3, nurbs.nu*nurbs.nv, nurbs.nu*nurbs.nv)
	if !nurbs.isHomogeneous {
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			Pw[i].x = nurbs.P[3*i]
			Pw[i].y = nurbs.P[3*i+1]
			Pw[i].z = nurbs.P[3*i+2]
			Pw[i].w = 1.0
		}
	} else {
		for i := 0; i < nurbs.nu*nurbs.nv; i++ {
			Pw[i].x = nurbs.P[4*i]
			Pw[i].y = nurbs.P[4*i+1]
			Pw[i].z = nurbs.P[4*i+2]
			Pw[i].w = nurbs.P[4*i+3]
		}
	}
	for v := 0; v < dicev; v++ {
		for u := 0; u < diceu; u++ {
			uvs[2*(v*diceu+u)] = ueval[u]
			uvs[2*(v*diceu+u)+1] = veval[v]

			pt, dPdu, dPdv := NURBSEvaluateSurface(nurbs.uorder, nurbs.uknot, nurbs.nu, ueval[u],
				nurbs.vorder, nurbs.vknot, nurbs.nv, veval[v], Pw)
			evalPs[v*diceu+u].X = pt.X
			evalPs[v*diceu+u].Y = pt.Y
			evalPs[v*diceu+u].Z = pt.Z
			evalNs[v*diceu+u] = *CreateNormalFromVector(NormalizeVector(CrossVector(dPdu, dPdv)))
		}
	}
	// Generate points-polygons mesh
	nTris := 2 * (diceu - 1) * (dicev - 1)
	vertices := make([]int, 3*nTris, 3*nTris)
	vi := 0
	VN := func(u, v int) int { return v*diceu + u }
	// Compute the vertex offset numbers for the triangles
	for v := 0; v < dicev-1; v++ {
		for u := 0; u < diceu-1; u++ {
			vertices[vi] = VN(u, v)
			vi++
			vertices[vi] = VN(u+1, v)
			vi++
			vertices[vi] = VN(u+1, v+1)
			vi++

			vertices[vi] = VN(u, v)
			vi++
			vertices[vi] = VN(u+1, v+1)
			vi++
			vertices[vi] = VN(u, v+1)
			vi++
		}
	}

	refined = append(refined, NewTriangleMesh(nurbs.objectToWorld,
		nurbs.worldToObject, nurbs.reverseOrientation, vertices, evalPs, evalNs, nil, uvs, nil))

	return refined
}
func (nurbs *NURBS) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (nurbs *NURBS) IntersectP(ray *Ray) bool {
	return false
}
func (nurbs *NURBS) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (nurbs *NURBS) Area() float64 {
	return 0.0
}
func (nurbs *NURBS) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (nurbs *NURBS) Pdf(pshape *Point) float64 {
	return 0.0
}
func (nurbs *NURBS) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (nurbs *NURBS) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (nurbs *NURBS) ObjectToWorld() *Transform {
	return nurbs.objectToWorld
}
func (nurbs *NURBS) WorldToObject() *Transform {
	return nurbs.worldToObject
}
func (nurbs *NURBS) ReverseOrientation() bool {
	return nurbs.reverseOrientation
}
func (nurbs *NURBS) TransformSwapsHandedness() bool {
	return nurbs.transformSwapsHandedness
}
func (nurbs *NURBS) ShapeId() uint32 {
	return nurbs.shapeId
}

// NURBS Evaluation Functions
func KnotOffset(knot []float64, order, np int, t float64) int {
	firstKnot := order - 1

	knotOffset := firstKnot
	for t > knot[knotOffset+1] {
		knotOffset++
	}
	Assert(knotOffset < np) // np == lastKnot
	Assert(t >= knot[knotOffset] && t <= knot[knotOffset+1])
	return knotOffset
}

func NURBSEvaluate(order int, knot []float64, cp []Homogeneous3, np int,
	cpStride int, t float64) (val Homogeneous3, deriv *Vector) {
	//    int nKnots = np + order;
	var alpha float64

	knotOffset := KnotOffset(knot, order, np, t)
	knot = knot[knotOffset:]

	cpOffset := knotOffset - order + 1
	Assert(cpOffset >= 0 && cpOffset < np)

	cpWork := make([]Homogeneous3, order, order)
	for i := 0; i < order; i++ {
		cpWork[i] = cp[(cpOffset+i)*cpStride]
	}
	for i := 0; i < order-2; i++ {
		for j := 0; j < order-1-i; j++ {
			alpha = (knot[1+j] - t) / (knot[1+j] - knot[j+2-order+i])
			Assert(alpha >= 0.0 && alpha <= 1.0)

			cpWork[j].x = cpWork[j].x*alpha + cpWork[j+1].x*(1-alpha)
			cpWork[j].y = cpWork[j].y*alpha + cpWork[j+1].y*(1-alpha)
			cpWork[j].z = cpWork[j].z*alpha + cpWork[j+1].z*(1-alpha)
			cpWork[j].w = cpWork[j].w*alpha + cpWork[j+1].w*(1-alpha)
		}
	}
	alpha = (knot[1] - t) / (knot[1] - knot[0])
	Assert(alpha >= 0.0 && alpha <= 1.0)

	val = Homogeneous3{cpWork[0].x*alpha + cpWork[1].x*(1-alpha),
		cpWork[0].y*alpha + cpWork[1].y*(1-alpha),
		cpWork[0].z*alpha + cpWork[1].z*(1-alpha),
		cpWork[0].w*alpha + cpWork[1].w*(1-alpha)}

	factor := float64(order-1) / (knot[1] - knot[0])
	delta := Homogeneous3{(cpWork[1].x - cpWork[0].x) * factor,
		(cpWork[1].y - cpWork[0].y) * factor,
		(cpWork[1].z - cpWork[0].z) * factor,
		(cpWork[1].w - cpWork[0].w) * factor}

	deriv.X = delta.x/val.w - (val.x * delta.w / (val.w * val.w))
	deriv.Y = delta.y/val.w - (val.y * delta.w / (val.w * val.w))
	deriv.Z = delta.z/val.w - (val.z * delta.w / (val.w * val.w))

	return val, deriv
}

func NURBSEvaluateSurface(uOrder int, uKnot []float64, ucp int, u float64,
	vOrder int, vKnot []float64, vcp int, v float64,
	cp []Homogeneous3) (Po *Point, dPdu, dPdv *Vector) {
	iso := make([]Homogeneous3, Maxi(uOrder, vOrder), Maxi(uOrder, vOrder))

	uOffset := KnotOffset(uKnot, uOrder, ucp, u)
	uFirstCp := uOffset - uOrder + 1
	Assert(uFirstCp >= 0 && uFirstCp+uOrder-1 < ucp)

	for i := 0; i < uOrder; i++ {
		iso[i], _ = NURBSEvaluate(vOrder, vKnot, cp[uFirstCp+i:], vcp, ucp, v)
	}
	vOffset := KnotOffset(vKnot, vOrder, vcp, v)
	vFirstCp := vOffset - vOrder + 1
	Assert(vFirstCp >= 0 && vFirstCp+vOrder-1 < vcp)

	var P Homogeneous3
	Unimplemented()
	//P, dPdu = NURBSEvaluate(uOrder, uKnot, iso - uFirstCp, ucp, 1, u)

	for i := 0; i < vOrder; i++ {
		iso[i], _ = NURBSEvaluate(uOrder, uKnot, cp[(vFirstCp+i)*ucp:], ucp, 1, u)
	}
	Unimplemented()
	//_, dPdv = NURBSEvaluate(vOrder, vKnot, iso - vFirstCp, vcp, 1, v)

	Po = CreatePoint(P.x/P.w, P.y/P.w, P.z/P.w)
	return Po, dPdu, dPdv
}
