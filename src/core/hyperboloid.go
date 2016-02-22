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

type Hyperboloid struct {
	ShapeData
	p1, p2                         Point
	Zmin, Zmax, phiMax, rmax, a, c float64
}

func NewHyperboloid(o2w, w2o *Transform, ro bool, point1, point2 Point, phimax float64) *Hyperboloid {
	h := new(Hyperboloid)
	h.objectToWorld = o2w
	h.worldToObject = w2o
	h.reverseOrientation = ro
	h.transformSwapsHandedness = SwapsHandednessTransform(h.objectToWorld)
	h.shapeId = GenerateShapeId()
	h.p1 = point1
	h.p2 = point2
	h.phiMax = Radians(Clamp(phimax, 0.0, 360.0))

	radius1 := math.Sqrt(h.p1.X*h.p1.X + h.p1.Y*h.p1.Y)
	radius2 := math.Sqrt(h.p2.X*h.p2.X + h.p2.Y*h.p2.Y)
	h.rmax = math.Max(radius1, radius2)
	h.Zmin = math.Min(h.p1.Z, h.p2.Z)
	h.Zmax = math.Max(h.p1.Z, h.p2.Z)

	// Compute implicit function coefficients for hyperboloid
	if h.p2.Z == 0.0 {
		h.p1, h.p2 = h.p2, h.p1
	}
	pp := h.p1
	var xy1, xy2 float64
	for {
		pp = *pp.Add((h.p2.Sub(&h.p1)).Scale(2.0))
		xy1 = pp.X*pp.X + pp.Y*pp.Y
		xy2 = h.p2.Y*h.p2.Y + h.p2.Y*h.p2.Y
		h.a = (1.0/xy1 - (pp.Z*pp.Z)/(xy1*h.p2.Z*h.p2.Z)) / (1.0 - (xy2*pp.Z*pp.Z)/(xy1*h.p2.Z*h.p2.Z))
		h.c = (h.a*xy2 - 1.0) / (h.p2.Z * h.p2.Z)
		if !math.IsInf(h.a, 0) && !math.IsNaN(h.a) {
			break
		}
	}

	return h
}

func (h *Hyperboloid) ObjectBound() *BBox {
	return &BBox{Point{-h.rmax, -h.rmax, h.Zmin}, Point{h.rmax, h.rmax, h.Zmax}}
}

func (h *Hyperboloid) WorldBound() *BBox {
	return BBoxTransform(h.objectToWorld, h.ObjectBound())
}

func (h *Hyperboloid) CanIntersect() bool {
	return true
}

func (h *Hyperboloid) Refine(refined []Shape) []Shape {
	return refined
}

func (h *Hyperboloid) Intersect(r RayBase) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := r.Transform(h.worldToObject)

	// Compute quadratic hyperboloid coefficients
	A := h.a*ray.Dir().X*ray.Dir().X + h.a*ray.Dir().Y*ray.Dir().Y - h.c*ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (h.a*ray.Dir().X*ray.Origin().X + h.a*ray.Dir().Y*ray.Origin().Y - h.c*ray.Dir().Z*ray.Origin().Z)
	C := h.a*ray.Origin().X*ray.Origin().X + h.a*ray.Origin().Y*ray.Origin().Y - h.c*ray.Origin().Z*ray.Origin().Z - 1

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); !ok {
		return false, 0.0, 0.0, nil
	}

	// Compute intersection distance along ray
	if t0 > ray.Maxt() || t1 < ray.Mint() {
		return false, 0.0, 0.0, nil
	}

	thit := t0
	if t0 < ray.Mint() {
		thit = t1
		if thit > ray.Maxt() {
			return false, 0.0, 0.0, nil
		}
	}

	// Compute hyperboloid inverse mapping
	phit := ray.PointAt(thit)
	v := (phit.Z - h.p1.Z) / (h.p2.Z - h.p1.Z)
	pr := h.p1.Scale(1.0 - v).Sub(h.p2.Scale(-v)) // using Sub(-v) rather than Add(v)
	phi := math.Atan2(pr.X*phit.Y-phit.X*pr.Y, phit.X*pr.X+phit.Y*pr.Y)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test hyperboloid intersection against clipping parameters
	if phit.Z < h.Zmin || phit.Z > h.Zmax || phi > h.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false, 0.0, 0.0, nil
		}
		// Compute hyperboloid inverse mapping
		phit = ray.PointAt(thit)
		v = (phit.Z - h.p1.Z) / (h.p2.Z - h.p1.Z)
		pr := h.p1.Scale(1.0 - v).Sub(h.p2.Scale(-v))
		phi = math.Atan2(pr.X*phit.Y-phit.X*pr.Y, phit.X*pr.X+phit.Y*pr.Y)
		if phi < 0 {
			phi += 2 * math.Pi
		}
		if phit.Z < h.Zmin || phit.Z > h.Zmax || phi > h.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Compute parametric representation of hyperboloid hit
	u := phi / h.phiMax

	// Compute hyperboloid $\dpdu$ and $\dpdv$
	cosphi, sinphi := math.Cos(phi), math.Sin(phi)
	dpdu := CreateVector(-h.phiMax*phit.Y, h.phiMax*phit.X, 0.0)
	dpdv := CreateVector((h.p2.X-h.p1.X)*cosphi-(h.p2.Y-h.p1.Y)*sinphi,
		(h.p2.X-h.p1.X)*sinphi+(h.p2.Y-h.p1.Y)*cosphi,
		h.p2.Z-h.p1.Z)

	// Compute hyperboloid $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.X, phit.Y, 0).Scale(-h.phiMax * h.phiMax)
	d2Pduv := CreateVector(-dpdv.Y, dpdv.X, 0).Scale(h.phiMax)
	d2Pdvv := CreateVector(0, 0, 0)

	// Compute coefficients for fundamental forms
	E := DotVector(dpdu, dpdu)
	F := DotVector(dpdu, dpdv)
	G := DotVector(dpdv, dpdv)
	N := NormalizeVector(CrossVector(dpdu, dpdv))
	e := DotVector(N, d2Pduu)
	f := DotVector(N, d2Pduv)
	g := DotVector(N, d2Pdvv)

	// Compute $\dndu$ and $\dndv$ from fundamental form coefficients
	invEGF2 := 1.0 / (E*G - F*F)
	dndu := CreateNormalFromVector(dpdu.Scale((f*F - e*G) * invEGF2).Add(dpdv.Scale((e*F - f*E) * invEGF2)))
	dndv := CreateNormalFromVector(dpdu.Scale((g*F - f*G) * invEGF2).Add(dpdv.Scale((f*F - g*E) * invEGF2)))

	// Initialize _DifferentialGeometry_ from parametric information
	dg = CreateDiffGeometry(PointTransform(h.objectToWorld, phit), VectorTransform(h.objectToWorld, dpdu), VectorTransform(h.objectToWorld, dpdv),
		NormalTransform(h.objectToWorld, dndu), NormalTransform(h.objectToWorld, dndv), u, v, h)

	// Update _tHit_ for quadric intersection
	tHit = thit

	// Compute _rayEpsilon_ for quadric intersection
	rayEpsilon = 5.0e-4 * tHit

	return true, tHit, rayEpsilon, dg
}

func (h *Hyperboloid) IntersectP(r RayBase) bool {
	// Transform _Ray_ to object space
	ray := r.Transform(h.worldToObject)

	// Compute quadratic hyperboloid coefficients
	A := h.a*ray.Dir().X*ray.Dir().X + h.a*ray.Dir().Y*ray.Dir().Y - h.c*ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (h.a*ray.Dir().X*ray.Origin().X + h.a*ray.Dir().Y*ray.Origin().Y - h.c*ray.Dir().Z*ray.Origin().Z)
	C := h.a*ray.Origin().X*ray.Origin().X + h.a*ray.Origin().Y*ray.Origin().Y - h.c*ray.Origin().Z*ray.Origin().Z - 1

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); !ok {
		return false
	}

	// Compute intersection distance along ray
	if t0 > ray.Maxt() || t1 < ray.Mint() {
		return false
	}

	thit := t0
	if t0 < ray.Mint() {
		thit = t1
		if thit > ray.Maxt() {
			return false
		}
	}

	// Compute hyperboloid inverse mapping
	phit := ray.PointAt(thit)
	v := (phit.Z - h.p1.Z) / (h.p2.Z - h.p1.Z)
	pr := h.p1.Scale(1.0 - v).Sub(h.p2.Scale(-v)) // using Sub(-v) rather than Add(v)
	phi := math.Atan2(pr.X*phit.Y-phit.X*pr.Y, phit.X*pr.X+phit.Y*pr.Y)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test hyperboloid intersection against clipping parameters
	if phit.Z < h.Zmin || phit.Z > h.Zmax || phi > h.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false
		}
		// Compute hyperboloid inverse mapping
		phit = ray.PointAt(thit)
		v = (phit.Z - h.p1.Z) / (h.p2.Z - h.p1.Z)
		pr := h.p1.Scale(1.0 - v).Sub(h.p2.Scale(-v))
		phi = math.Atan2(pr.X*phit.Y-phit.X*pr.Y, phit.X*pr.X+phit.Y*pr.Y)
		if phi < 0 {
			phi += 2 * math.Pi
		}
		if phit.Z < h.Zmin || phit.Z > h.Zmax || phi > h.phiMax {
			return false
		}
	}
	return true
}

func (h *Hyperboloid) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (h *Hyperboloid) Area() float64 {
	SQR := func(a float64) float64 { return a * a }
	QUAD := func(a float64) float64 { return SQR(a) * SQR(a) }

	return h.phiMax / 6.0 *
		(2.0*QUAD(h.p1.X) - 2.0*h.p1.X*h.p1.X*h.p1.X*h.p2.X +
			2.0*QUAD(h.p2.X) +
			2.0*(h.p1.Y*h.p1.Y+h.p1.Y*h.p2.Y+h.p2.Y*h.p2.Y)*
				(SQR(h.p1.Y-h.p2.Y)+SQR(h.p1.Z-h.p2.Z)) +
			h.p2.X*h.p2.X*(5.0*h.p1.Y*h.p1.Y+2.0*h.p1.Y*h.p2.Y-
				4.0*h.p2.Y*h.p2.Y+2.0*SQR(h.p1.Z-h.p2.Z)) +
			h.p1.X*h.p1.X*(-4.0*h.p1.Y*h.p1.Y+2.0*h.p1.Y*h.p2.Y+
				5.0*h.p2.Y*h.p2.Y+2.0*SQR(h.p1.Z-h.p2.Z)) -
			2.0*h.p1.X*h.p2.X*(h.p2.X*h.p2.X-h.p1.Y*h.p1.Y+
				5.0*h.p1.Y*h.p2.Y-h.p2.Y*h.p2.Y-h.p1.Z*h.p1.Z+
				2.0*h.p1.Z*h.p2.Z-h.p2.Z*h.p2.Z))
}

func (h *Hyperboloid) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}

func (h *Hyperboloid) Pdf(pshape *Point) float64 {
	return 1.0 / h.Area()
}

func (h *Hyperboloid) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}

func (h *Hyperboloid) Pdf2(p *Point, wi *Vector) float64 {
	return ShapePdf(h, p, wi)
}

func (h *Hyperboloid) ObjectToWorld() *Transform {
	return h.objectToWorld
}

func (h *Hyperboloid) WorldToObject() *Transform {
	return h.worldToObject
}

func (h *Hyperboloid) ReverseOrientation() bool {
	return h.reverseOrientation
}

func (h *Hyperboloid) TransformSwapsHandedness() bool {
	return h.transformSwapsHandedness
}

func (h *Hyperboloid) ShapeId() uint32 {
	return h.shapeId
}
