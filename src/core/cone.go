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

type Cone struct {
	ShapeData
	radius, height, phiMax float64
}

func NewCone(o2w, w2o *Transform, ro bool, rad, h, pm float64) *Cone {
	c := new(Cone)
	c.objectToWorld = o2w
	c.worldToObject = w2o
	c.reverseOrientation = ro
	c.transformSwapsHandedness = SwapsHandednessTransform(c.objectToWorld)
	c.shapeId = GenerateShapeId()
	c.radius = rad
	c.height = h
	c.phiMax = Radians(Clamp(pm, 0.0, 360.0))

	return c
}

func (c *Cone) String() string {
	return fmt.Sprintf("cone[r: %f h: %f phimax: %f obj2world: %v]", c.radius, c.height, Degrees(c.phiMax), c.objectToWorld)
}

func (c *Cone) ObjectBound() *BBox {
	return &BBox{Point{-c.radius, -c.radius, 0.0}, Point{c.radius, c.radius, c.height}}
}

func (c *Cone) WorldBound() *BBox {
	return BBoxTransform(c.objectToWorld, c.ObjectBound())
}

func (c *Cone) CanIntersect() bool {
	return true
}

func (c *Cone) Refine(refined []Shape) []Shape {
	return refined
}

func (c *Cone) Intersect(r RayBase) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := r.Transform(c.worldToObject)

	// Compute quadratic cone coefficients
	k := c.radius / c.height
	k = k * k
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y - k*ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y - k*ray.Dir().Z*(ray.Origin().Z-c.height))
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y - k*(ray.Origin().Z-c.height)*(ray.Origin().Z-c.height)

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

	// Compute cone inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cone intersection against clipping parameters
	if phit.Z < 0 || phit.Z > c.height || phi > c.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false, 0.0, 0.0, nil
		}
		// Compute cone inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < 0 || phit.Z > c.height || phi > c.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of cone hit
	u := phi / c.phiMax
	v := phit.Z / c.height

	// Compute cone $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-c.phiMax*phit.Y, c.phiMax*phit.X, 0.0)
	dpdv := CreateVector(-phit.X/(1.0-v), -phit.Y/(1.0-v), c.height)

	// Compute cone $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.X, phit.Y, 0.0).Scale(-c.phiMax * c.phiMax)
	d2Pduv := CreateVector(phit.Y, -phit.X, 0.0).Scale(c.phiMax / (1.0 - v))

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

func (c *Cone) IntersectP(r RayBase) bool {
	// Transform _Ray_ to object space
	ray := r.Transform(c.worldToObject)

	// Compute quadratic cone coefficients
	k := c.radius / c.height
	k = k * k
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y - k*ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y - k*ray.Dir().Z*(ray.Origin().Z-c.height))
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y - k*(ray.Origin().Z-c.height)*(ray.Origin().Z-c.height)

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

	// Compute cone inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cone intersection against clipping parameters
	if phit.Z < 0 || phit.Z > c.height || phi > c.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.Maxt() {
			return false
		}
		// Compute cone inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < 0 || phit.Z > c.height || phi > c.phiMax {
			return false
		}
	}
	return true
}

func (c *Cone) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (c *Cone) Area() float64 {
	return c.radius * math.Sqrt((c.height*c.height)+(c.radius*c.radius)) * c.phiMax / 2.0
}

func (c *Cone) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}

func (c *Cone) Pdf(pshape *Point) float64 {
	return 1.0 / c.Area()
}

func (c *Cone) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return c.Sample(u1, u2)
}

func (c *Cone) Pdf2(p *Point, wi *Vector) float64 {
	return ShapePdf(c, p, wi)
}

func (c *Cone) ObjectToWorld() *Transform {
	return c.objectToWorld
}

func (c *Cone) WorldToObject() *Transform {
	return c.worldToObject
}

func (c *Cone) ReverseOrientation() bool {
	return c.reverseOrientation
}

func (c *Cone) TransformSwapsHandedness() bool {
	return c.transformSwapsHandedness
}

func (c *Cone) ShapeId() uint32 {
	return c.shapeId
}
