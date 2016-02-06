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
	"fmt"
)

type TriangleMesh struct {
	ShapeData
	ntris, nverts int
	vertexIndex   []int
	p             []Point
	n             []Normal
	s             []Vector
	uvs           []float64
	alphaTexture  TextureFloat
}

type Triangle struct {
	ShapeData
	mesh *TriangleMesh
	v    []int
}

func CreateTriangleMeshShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet, textures map[string]TextureFloat) *TriangleMesh {
	vi := params.FindIntArrayParam("indices")
	P := params.FindPointArrayParam("P")
	uvs := params.FindFloatArrayParam("uv")
	if uvs == nil {
		uvs = params.FindFloatArrayParam("st")
	}
	discardDegnerateUVs := params.FindBoolParam("discarddegenerateUVs", false)

	// XXX should complain if uvs aren't an array of 2...
	if uvs != nil {
		if len(uvs) < 2*len(P) {
			Error("Not enough of \"uv\"s for triangle mesh.  Expencted %d, found %d.  Discarding.", 2*len(P), len(uvs))
			uvs = nil
		} else if len(uvs) > 2*len(P) {
			Warning("More \"uv\"s provided than will be used for triangle  mesh.  (%d expcted, %d found)", 2*len(P), len(uvs))
		}
	}
	if vi == nil || P == nil {
		return nil
	}

	S := params.FindVectorArrayParam("S")
	if S != nil && len(S) != len(P) {
		Error("Number of \"S\"s for triangle mesh must match \"P\"s")
		S = nil
	}
	N := params.FindNormalArrayParam("N")
	if N != nil && len(N) != len(P) {
		Error("Number of \"N\"s for triangle mesh must match \"P\"s")
		N = nil
	}
	if discardDegnerateUVs && uvs != nil && N != nil {
		// if there are normals, check for bad uv's that
		// give degenerate mappings; discard them if so
		for i := 0; i < len(vi); i = i + 3 {
			area := 0.5 * CrossVector(P[vi[i+0]].Sub(&P[vi[i+1]]), P[vi[i+2]].Sub(&P[vi[i+1]])).Length()
			if area < 1.0e-7 {
				continue
			} // ignore degenerate tris.
			if (uvs[2*vi[i+0]] == uvs[2*vi[i+1]] &&
				uvs[2*vi[i+0]+1] == uvs[2*vi[i+1]+1]) ||
				(uvs[2*vi[i+1]] == uvs[2*vi[i+2]] &&
					uvs[2*vi[i+1]+1] == uvs[2*vi[i+2]+1]) ||
				(uvs[2*vi[i+2]] == uvs[2*vi[i+0]] &&
					uvs[2*vi[i+2]+1] == uvs[2*vi[i+0]+1]) {
				Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.")
				uvs = nil
				break
			}
		}
	}

	for i := 0; i < len(vi); i++ {
		if vi[i] >= len(P) {
			Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
				vi[i], len(P))
			return nil
		}
	}

	var alphaTex TextureFloat
	alphaTexName := params.FindTextureParam("alpha")
	if len(alphaTexName) != 0 {
		alphaTex = textures[alphaTexName]
		if alphaTex == nil {
			Error("Couldn't find float texture \"%s\" for \"alpha\" parameter", alphaTexName)
		}
	} else if params.FindFloatParam("alpha", 1.0) == 0.0 {
		alphaTex = &ConstantTextureFloat{0.0}
	}

	return NewTriangleMesh(o2w, w2o, reverseOrientation, vi, P,
		N, S, uvs, alphaTex)
}

func NewTriangleMesh(o2w, w2o *Transform, reverseOrientation bool, vi []int, P []Point, N []Normal, S []Vector, uv []float64, atex TextureFloat) *TriangleMesh {
	t := new(TriangleMesh)
	t.objectToWorld = o2w
	t.worldToObject = w2o
	t.reverseOrientation = reverseOrientation
	t.transformSwapsHandedness = SwapsHandednessTransform(t.objectToWorld)
	t.shapeId = GenerateShapeId()
	t.alphaTexture = atex
	t.ntris = len(vi) / 3
	t.nverts = len(P)
	t.vertexIndex = make([]int, 3*t.ntris, 3*t.ntris)
	copy(t.vertexIndex, vi)

	if uv != nil && len(uv) == 2*t.nverts {
		t.uvs = make([]float64, 2*t.nverts, 2*t.nverts)
		copy(t.uvs, uv)
	} else {
		t.uvs = nil
	}
	t.p = make([]Point, t.nverts, t.nverts)
	copy(t.p, P)
	if N != nil {
		t.n = make([]Normal, t.nverts, t.nverts)
		copy(t.n, N)
	} else {
		t.n = nil
	}
	if S != nil {
		t.s = make([]Vector, t.nverts, t.nverts)
		copy(t.s, S)
	} else {
		t.s = nil
	}
	// Transform mesh vertices to world space
	for i, p := range P {
		t.p[i] = *PointTransform(t.objectToWorld, &p)
	}
	return t
}

func (t *TriangleMesh) ObjectBound() *BBox {
	objectBounds := CreateEmptyBBox()
	for i := 0; i < t.nverts; i++ {
		objectBounds = UnionBBoxPoint(objectBounds, PointTransform(t.worldToObject, &t.p[i]))
	}
	return objectBounds
}

func (t *TriangleMesh) WorldBound() *BBox {
	worldBounds := CreateEmptyBBox()
	if t == nil {
		panic("Null mesh arg.")
	}
	for i := 0; i < t.nverts; i++ {
		worldBounds = UnionBBoxPoint(worldBounds, &t.p[i])
	}
	return worldBounds
}

func (t *TriangleMesh) CanIntersect() bool {
	return false
}
func (t *TriangleMesh) Refine(refined []Shape) []Shape {
	for i := 0; i < t.ntris; i++ {
		refined = append(refined, NewTriangle(t.objectToWorld, t.worldToObject, t.reverseOrientation, t, i))
	}
	return refined
}

func (t *TriangleMesh) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (t *TriangleMesh) IntersectP(ray *Ray) bool {
	return false
}
func (t *TriangleMesh) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}
func (t *TriangleMesh) Area() float64 {
	return 0.0
}
func (t *TriangleMesh) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *TriangleMesh) Pdf(pshape *Point) float64 {
	return 0.0
}
func (t *TriangleMesh) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *TriangleMesh) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (t *TriangleMesh) ObjectToWorld() *Transform {
	return t.objectToWorld
}
func (t *TriangleMesh) WorldToObject() *Transform {
	return t.worldToObject
}
func (t *TriangleMesh) ReverseOrientation() bool {
	return t.reverseOrientation
}
func (t *TriangleMesh) TransformSwapsHandedness() bool {
	return t.transformSwapsHandedness
}
func (t *TriangleMesh) ShapeId() uint32 {
	return t.shapeId
}

func (t *TriangleMesh) String() string {
	return fmt.Sprintf("mesh[id: %d ntris: %d nverts: %d]", t.shapeId, t.ntris, t.nverts)
}

func NewTriangle(o2w, w2o *Transform, ro bool, m *TriangleMesh, n int) *Triangle {
	t := new(Triangle)
	t.objectToWorld = o2w
	t.worldToObject = w2o
	t.reverseOrientation = ro
	t.transformSwapsHandedness = SwapsHandednessTransform(t.objectToWorld)
	t.shapeId = GenerateShapeId()
	t.mesh = m
	t.v = m.vertexIndex[3*n : 3*(n+1)]
	return t
}

func (t *Triangle) ObjectBound() *BBox {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &t.mesh.p[t.v[0]]
	p2 := &t.mesh.p[t.v[1]]
	p3 := &t.mesh.p[t.v[2]]
	return UnionBBoxPoint(CreateBBoxFromPoints(PointTransform(t.worldToObject, p1), PointTransform(t.worldToObject, p2)), PointTransform(t.worldToObject, p3))
}

func (t *Triangle) WorldBound() *BBox {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &t.mesh.p[t.v[0]]
	p2 := &t.mesh.p[t.v[1]]
	p3 := &t.mesh.p[t.v[2]]
	return UnionBBoxPoint(CreateBBoxFromPoints(p1, p2), p3)
}

func (t *Triangle) CanIntersect() bool {
	return true
}

func (t *Triangle) Refine(refined []Shape) []Shape {
	return refined
}

func (tri *Triangle) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	//PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this))
	// Compute $\VEC{s}_1$

	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &tri.mesh.p[tri.v[0]]
	p2 := &tri.mesh.p[tri.v[1]]
	p3 := &tri.mesh.p[tri.v[2]]
	e1 := p2.Sub(p1)
	e2 := p3.Sub(p1)
	s1 := CrossVector(&ray.Dir, e2)
	divisor := DotVector(s1, e1)

	if divisor == 0.0 {
		return false, 0.0, 0.0, nil
	}

	invDivisor := 1.0 / divisor

	// Compute first barycentric coordinate
	s := ray.Origin.Sub(p1)
	b1 := DotVector(s, s1) * invDivisor
	if b1 < 0.0 || b1 > 1.0 {
		return false, 0.0, 0.0, nil
	}

	// Compute second barycentric coordinate
	s2 := CrossVector(s, e1)
	b2 := DotVector(&ray.Dir, s2) * invDivisor
	if b2 < 0.0 || b1+b2 > 1.0 {
		return false, 0.0, 0.0, nil
	}

	// Compute _t_ to intersection point
	t := DotVector(e2, s2) * invDivisor
	if t < ray.Mint || t > ray.Maxt {
		return false, 0.0, 0.0, nil
	}

	// Compute triangle partial derivatives
	var dpdu, dpdv *Vector
	uvs := tri.GetUVs()

	// Compute deltas for triangle partial derivatives
	du1 := uvs[0][0] - uvs[2][0]
	du2 := uvs[1][0] - uvs[2][0]
	dv1 := uvs[0][1] - uvs[2][1]
	dv2 := uvs[1][1] - uvs[2][1]
	dp1 := p1.Sub(p3)
	dp2 := p2.Sub(p3)
	determinant := du1*dv2 - dv1*du2
	if determinant == 0.0 {
		// Handle zero determinant for triangle partial derivative matrix
		dpdu, dpdv = CoordinateSystem(NormalizeVector(CrossVector(e2, e1)))
	} else {
		invdet := 1.0 / determinant
		dpdu = dp1.Scale(dv2).Sub(dp2.Scale(dv1)).Scale(invdet)
		dpdv = dp1.Scale(-du2).Add(dp2.Scale(du1)).Scale(invdet)
	}

	// Interpolate $(u,v)$ triangle parametric coordinates
	b0 := 1 - b1 - b2
	tu := b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0]
	tv := b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1]

	// Test intersection against alpha texture, if present
	if ray.Depth != -1 {
		if tri.mesh.alphaTexture != nil {
			dgLocal := CreateDiffGeometry(ray.PointAt(t), dpdu, dpdv,
				CreateNormal(0, 0, 0), CreateNormal(0, 0, 0),
				tu, tv, tri)
			if tri.mesh.alphaTexture.Evaluate(dgLocal) == 0.0 {
				return false, 0.0, 0.0, nil
			}
		}
	}

	// Fill in _DifferentialGeometry_ from triangle hit
	dg = CreateDiffGeometry(ray.PointAt(t), dpdu, dpdv,
		CreateNormal(0, 0, 0), CreateNormal(0, 0, 0),
		tu, tv, tri)
	tHit = t
	rayEpsilon = 1.0e-3 * tHit
	//PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
	return true, tHit, rayEpsilon, dg
}

func (tri *Triangle) IntersectP(ray *Ray) bool {
	//PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
	// Compute $\VEC{s}_1$

	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &tri.mesh.p[tri.v[0]]
	p2 := &tri.mesh.p[tri.v[1]]
	p3 := &tri.mesh.p[tri.v[2]]
	e1 := p2.Sub(p1)
	e2 := p3.Sub(p1)
	s1 := CrossVector(&ray.Dir, e2)
	divisor := DotVector(s1, e1)

	if divisor == 0.0 {
		return false
	}

	invDivisor := 1.0 / divisor

	// Compute first barycentric coordinate
	s := ray.Origin.Sub(p1)
	b1 := DotVector(s, s1) * invDivisor
	if b1 < 0.0 || b1 > 1.0 {
		return false
	}

	// Compute second barycentric coordinate
	s2 := CrossVector(s, e1)
	b2 := DotVector(&ray.Dir, s2) * invDivisor
	if b2 < 0.0 || b1+b2 > 1.0 {
		return false
	}

	// Compute _t_ to intersection point
	t := DotVector(e2, s2) * invDivisor
	if t < ray.Mint || t > ray.Maxt {
		return false
	}

	// Test shadow ray intersection against alpha texture, if present
	if ray.Depth != -1 && tri.mesh.alphaTexture != nil {
		// Compute triangle partial derivatives
		var dpdu, dpdv *Vector
		uvs := tri.GetUVs()

		// Compute deltas for triangle partial derivatives
		du1 := uvs[0][0] - uvs[2][0]
		du2 := uvs[1][0] - uvs[2][0]
		dv1 := uvs[0][1] - uvs[2][1]
		dv2 := uvs[1][1] - uvs[2][1]
		dp1 := p1.Sub(p3)
		dp2 := p2.Sub(p3)
		determinant := du1*dv2 - dv1*du2
		if determinant == 0.0 {
			// Handle zero determinant for triangle partial derivative matrix
			dpdu, dpdv = CoordinateSystem(NormalizeVector(CrossVector(e2, e1)))
		} else {
			invdet := 1.0 / determinant
			dpdu = dp1.Scale(dv2).Sub(dp2.Scale(dv1)).Scale(invdet)
			dpdv = dp1.Scale(-du2).Add(dp2.Scale(du1)).Scale(invdet)
		}

		// Interpolate $(u,v)$ triangle parametric coordinates
		b0 := 1 - b1 - b2
		tu := b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0]
		tv := b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1]
		dgLocal := CreateDiffGeometry(ray.PointAt(t), dpdu, dpdv,
			CreateNormal(0, 0, 0), CreateNormal(0, 0, 0),
			tu, tv, tri)
		if tri.mesh.alphaTexture.Evaluate(dgLocal) == 0.0 {
			return false
		}
	}
	//PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);
	return true
}

func (tri *Triangle) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	if tri.mesh.n == nil && tri.mesh.s == nil {
		return dg
	}
	// Initialize _Triangle_ shading geometry with _n_ and _s_

	// Compute barycentric coordinates for point
	var b [3]float64

	// Initialize _A_ and _C_ matrices for barycentrics
	uv := tri.GetUVs()

	A := [2][2]float64{{uv[1][0] - uv[0][0], uv[2][0] - uv[0][0]},
		{uv[1][1] - uv[0][1], uv[2][1] - uv[0][1]}}
	C := [2]float64{dg.u - uv[0][0], dg.v - uv[0][1]}
	var ok bool
	if ok, b[1], b[2] = SolveLinearSystem2x2(A, C); !ok {
		// Handle degenerate parametric mapping
		b[0], b[1], b[2] = 1.0/3.0, 1.0/3.0, 1.0/3.0
	} else {
		b[0] = 1.0 - b[1] - b[2]
	}

	// Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
	var ns *Normal
	var ss *Vector
	if tri.mesh.n != nil {
		ns = NormalizeNormal(NormalTransform(obj2world, tri.mesh.n[tri.v[0]].Scale(b[0]).Add(tri.mesh.n[tri.v[1]].Scale(b[1]).Add(tri.mesh.n[tri.v[2]].Scale(b[2])))))
	} else {
		ns = dg.nn
	}
	if tri.mesh.s != nil {
		ss = NormalizeVector(VectorTransform(obj2world, tri.mesh.s[tri.v[0]].Scale(b[0]).Add(tri.mesh.s[tri.v[1]].Scale(b[1]).Add(tri.mesh.s[tri.v[2]].Scale(b[2])))))
	} else {
		ss = NormalizeVector(dg.dpdu)
	}

	ts := CrossVectorNormal(ss, ns)
	if ts.LengthSquared() > 0.0 {
		ts = NormalizeVector(ts)
		ss = CrossVectorNormal(ts, ns)
	} else {
		ss, ts = CoordinateSystem(CreateVectorFromNormal(ns))
	}
	var dndu, dndv *Normal

	// Compute $\dndu$ and $\dndv$ for triangle shading geometry
	if tri.mesh.n != nil {
		uvs := tri.GetUVs()
		// Compute deltas for triangle partial derivatives of normal
		du1 := uvs[0][0] - uvs[2][0]
		du2 := uvs[1][0] - uvs[2][0]
		dv1 := uvs[0][1] - uvs[2][1]
		dv2 := uvs[1][1] - uvs[2][1]
		dn1 := tri.mesh.n[tri.v[0]].Sub(&tri.mesh.n[tri.v[2]])
		dn2 := tri.mesh.n[tri.v[1]].Sub(&tri.mesh.n[tri.v[2]])
		determinant := du1*dv2 - dv1*du2
		if determinant == 0.0 {
			dndu, dndv = CreateNormal(0, 0, 0), CreateNormal(0, 0, 0)
		} else {
			invdet := 1.0 / determinant
			dndu = dn1.Scale(dv2).Sub(dn2.Scale(dv1)).Scale(invdet)
			dndv = dn1.Scale(-du2).Add(dn2.Scale(du1)).Scale(invdet)
		}
	} else {
		dndu, dndv = CreateNormal(0, 0, 0), CreateNormal(0, 0, 0)
	}
	dgShading := CreateDiffGeometry(dg.p, ss, ts, NormalTransform(obj2world, dndu), NormalTransform(obj2world, dndv), dg.u, dg.v, dg.shape)
	dgShading.dudx = dg.dudx
	dgShading.dvdx = dg.dvdx
	dgShading.dudy = dg.dudy
	dgShading.dvdy = dg.dvdy
	dgShading.dpdx = dg.dpdx
	dgShading.dpdy = dg.dpdy
	return dgShading
}

func (t *Triangle) Area() float64 {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &t.mesh.p[t.v[0]]
	p2 := &t.mesh.p[t.v[1]]
	p3 := &t.mesh.p[t.v[2]]
	return 0.5 * CrossVector(p2.Sub(p1), p3.Sub(p1)).Length()
}

func (t *Triangle) Sample(u1, u2 float64) (ps *Point, ns *Normal) {
	b1, b2 := UniformSampleTriangle(u1, u2)
	// Get triangle vertices in _p1_, _p2_, and _p3_
	p1 := &t.mesh.p[t.v[0]]
	p2 := &t.mesh.p[t.v[1]]
	p3 := &t.mesh.p[t.v[2]]
	ps = p1.Scale(b1).Add(p2.Scale(b2).Sub(p3.Scale(b1 + b2 - 1.0)))
	n := CreateNormalFromVector(CrossVector(p2.Sub(p1), p3.Sub(p1)))
	ns = NormalizeNormal(n)
	if t.reverseOrientation {
		ns = ns.Scale(-1.0)
	}
	return ps, ns
}

func (t *Triangle) Pdf(pshape *Point) float64 {
	return 0.0
}

func (t *Triangle) SampleAt(p *Point, u1, u2 float64) (ps *Point, ns *Normal) {
	return t.Sample(u1, u2)
}

func (t *Triangle) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (t *Triangle) ObjectToWorld() *Transform {
	return t.objectToWorld
}
func (t *Triangle) WorldToObject() *Transform {
	return t.worldToObject
}
func (t *Triangle) ReverseOrientation() bool {
	return t.reverseOrientation
}
func (t *Triangle) TransformSwapsHandedness() bool {
	return t.transformSwapsHandedness
}
func (t *Triangle) ShapeId() uint32 {
	return t.mesh.shapeId
}

func (t *Triangle) GetUVs() (uv [3][2]float64) {
	if t.mesh.uvs != nil {
		uv[0][0] = t.mesh.uvs[2*t.v[0]]
		uv[0][1] = t.mesh.uvs[2*t.v[0]+1]
		uv[1][0] = t.mesh.uvs[2*t.v[1]]
		uv[1][1] = t.mesh.uvs[2*t.v[1]+1]
		uv[2][0] = t.mesh.uvs[2*t.v[2]]
		uv[2][1] = t.mesh.uvs[2*t.v[2]+1]
	} else {
		uv[0][0] = 0.0
		uv[0][1] = 0.0
		uv[1][0] = 1.0
		uv[1][1] = 0.0
		uv[2][0] = 1.0
		uv[2][1] = 1.0
	}
	return uv
}
