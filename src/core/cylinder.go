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
	"math"
)

type Cylinder struct {
	ShapeData
	radius, Zmin, Zmax, phiMax float64
}

func NewCylinder(o2w, w2o *Transform, ro bool, rad, zmin, zmax, pmax float64) *Cylinder {
	c := new(Cylinder)
	c.objectToWorld = o2w
	c.worldToObject = w2o
	c.reverseOrientation = ro
	c.transformSwapsHandedness = SwapsHandednessTransform(c.objectToWorld)
	c.shapeId = GenerateShapeId()
	c.radius = rad
	c.Zmin = zmin
	c.Zmax = zmax
	c.phiMax = Radians(Clamp(pmax, 0.0, 360.0))

	return c
}

func (c *Cylinder) String() string {
	return fmt.Sprintf("cylinder[r: %f zmin: %f zmax: %f phimax: %f obj2world: %v]",
		c.radius, c.Zmin, c.Zmax, Degrees(c.phiMax), c.objectToWorld)
}

func (c *Cylinder) ObjectBound() *BBox {
	return &BBox{Point{-c.radius, -c.radius, c.Zmin}, Point{c.radius, c.radius, c.Zmax}}
}

func (c *Cylinder) WorldBound() *BBox {
	return BBoxTransform(c.objectToWorld, c.ObjectBound())
}

func (c *Cylinder) CanIntersect() bool {
	return true
}

func (c *Cylinder) Refine(refined []Shape) []Shape {
	return refined
}

func (c *Cylinder) Intersect(r RayBase) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := r.Transform(c.worldToObject)

	// Compute quadratic cylinder coefficients
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y)
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y - c.radius*c.radius

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

	// Compute cylinder hit point and $\phi$
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cylinder intersection against clipping parameters
	if phit.Z < c.Zmin || phit.Z > c.Zmax || phi > c.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false, 0.0, 0.0, nil
		}
		// Compute cylinder hit point and $\phi$
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < c.Zmin || phit.Z > c.Zmax || phi > c.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of cylinder hit
	u := phi / c.phiMax
	v := (phit.Z - c.Zmin) / (c.Zmax - c.Zmin)

	// Compute cylinder $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-c.phiMax*phit.Y, c.phiMax*phit.X, 0.0)
	dpdv := CreateVector(0.0, 0.0, c.Zmax-c.Zmin)

	// Compute cylinder $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.X, phit.Y, 0).Scale(-c.phiMax * c.phiMax)
	d2Pduv := CreateVector(0, 0, 0)
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
	dg = CreateDiffGeometry(PointTransform(c.objectToWorld, phit), VectorTransform(c.objectToWorld, dpdu), VectorTransform(c.objectToWorld, dpdv),
		NormalTransform(c.objectToWorld, dndu), NormalTransform(c.objectToWorld, dndv), u, v, c)

	// Update _tHit_ for quadric intersection
	tHit = thit

	// Compute _rayEpsilon_ for quadric intersection
	rayEpsilon = 5.0e-4 * tHit

	return true, tHit, rayEpsilon, dg
}

func (c *Cylinder) IntersectP(r RayBase) bool {
	// Transform _Ray_ to object space
	ray := r.Transform(c.worldToObject)

	// Compute quadratic cylinder coefficients
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y)
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y - c.radius*c.radius

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

	// Compute cylinder hit point and $\phi$
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cylinder intersection against clipping parameters
	if phit.Z < c.Zmin || phit.Z > c.Zmax || phi > c.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false
		}
		// Compute cylinder hit point and $\phi$
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < c.Zmin || phit.Z > c.Zmax || phi > c.phiMax {
			return false
		}
	}
	return true
}

func (c *Cylinder) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (c *Cylinder) Area() float64 {
	return (c.Zmax - c.Zmin) * c.phiMax * c.radius
}

func (c *Cylinder) Sample(u1, u2 float64) (*Point, *Normal) {
	z := Lerp(u1, c.Zmin, c.Zmax)
	t := u2 * c.phiMax
	p := CreatePoint(c.radius*math.Cos(t), c.radius*math.Sin(t), z)
	Ns := NormalizeNormal(NormalTransform(c.objectToWorld, CreateNormal(p.X, p.Y, 0.0)))
	if c.reverseOrientation {
		Ns = Ns.Negate()
	}
	return PointTransform(c.objectToWorld, p), Ns
}

func (c *Cylinder) Pdf(pshape *Point) float64 {
	return 1.0 / c.Area()
}

func (c *Cylinder) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return c.Sample(u1, u2)
}

func (c *Cylinder) Pdf2(p *Point, wi *Vector) float64 {
	return ShapePdf(c, p, wi)
}

func (c *Cylinder) ObjectToWorld() *Transform {
	return c.objectToWorld
}

func (c *Cylinder) WorldToObject() *Transform {
	return c.worldToObject
}

func (c *Cylinder) ReverseOrientation() bool {
	return c.reverseOrientation
}

func (c *Cylinder) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}

func (c *Cylinder) ShapeId() uint32 {
	return c.shapeId
}
