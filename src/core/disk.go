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

type Disk struct {
	ShapeData
	height, radius, innerRadius, phiMax float64
}

func NewDisk(o2w, w2o *Transform, ro bool, height, radius, innerradius, phimax float64) *Disk {
	d := new(Disk)
	d.objectToWorld = o2w
	d.worldToObject = w2o
	d.reverseOrientation = ro
	d.transformSwapsHandedness = SwapsHandednessTransform(d.objectToWorld)
	d.shapeId = GenerateShapeId()
	d.height = height
	d.radius = radius
	d.innerRadius = innerradius
	d.phiMax = Radians(Clamp(phimax, 0.0, 360.0))

	return d
}

func (d *Disk) String() string {
	return fmt.Sprintf("disk[h: %f rad: %f, %f phimax: %f obj2world: %v]", d.height, d.radius, d.innerRadius, Degrees(d.phiMax), d.objectToWorld)
}

func (d *Disk) ObjectBound() *BBox {
	return &BBox{Point{-d.radius, -d.radius, d.height}, Point{d.radius, d.radius, d.height}}
}

func (d *Disk) WorldBound() *BBox {
	return BBoxTransform(d.objectToWorld, d.ObjectBound())
}

func (d *Disk) CanIntersect() bool {
	return true
}

func (d *Disk) Refine(refined []Shape) []Shape {
	return refined
}

func (d *Disk) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(d.worldToObject, r)

	// Compute plane intersection for disk
	if math.Abs(ray.Dir.Z) < 1.0e-7 {
		return false, 0.0, 0.0, nil
	}
	thit := (d.height - ray.Origin.Z) / ray.Dir.Z
	if thit < ray.Mint || thit > ray.Maxt {
		return false, 0.0, 0.0, nil
	}

	// See if hit point is inside disk radii and $\phimax$
	phit := ray.PointAt(thit)
	dist2 := phit.X*phit.X + phit.Y*phit.Y
	if dist2 > d.radius*d.radius || dist2 < d.innerRadius*d.innerRadius {
		return false, 0.0, 0.0, nil
	}

	// Test disk $\phi$ value against $\phimax$
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	if phi > d.phiMax {
		return false, 0.0, 0.0, nil
	}

	// Find parametric representation of disk hit
	u := phi / d.phiMax
	R := math.Sqrt(dist2)
	oneMinusV := ((R - d.innerRadius) / (d.radius - d.innerRadius))
	v := 1.0 - oneMinusV

	dpdu := CreateVector(-d.phiMax*phit.Y, d.phiMax*phit.X, 0.0)
	dpdv := CreateVector(phit.X, phit.Y, 0.0)
	dpdv.Scale((d.innerRadius - d.radius) / R)
	dndu := CreateNormal(0, 0, 0)
	dndv := CreateNormal(0, 0, 0)

	// Initialize _DifferentialGeometry_ from parametric information
	dg = CreateDiffGeometry(PointTransform(d.objectToWorld, phit), VectorTransform(d.objectToWorld, dpdu), VectorTransform(d.objectToWorld, dpdv),
		NormalTransform(d.objectToWorld, dndu), NormalTransform(d.objectToWorld, dndv), u, v, d)

	// Update _tHit_ for quadric intersection
	tHit = thit

	// Compute _rayEpsilon_ for quadric intersection
	rayEpsilon = 5.0e-4 * tHit

	return true, tHit, rayEpsilon, dg
}

func (d *Disk) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(d.worldToObject, r)

	// Compute plane intersection for disk
	if math.Abs(ray.Dir.Z) < 1.0e-7 {
		return false
	}
	thit := (d.height - ray.Origin.Z) / ray.Dir.Z
	if thit < ray.Mint || thit > ray.Maxt {
		return false
	}

	// See if hit point is inside disk radii and $\phimax$
	phit := ray.PointAt(thit)
	dist2 := phit.X*phit.X + phit.Y*phit.Y
	if dist2 > d.radius*d.radius || dist2 < d.innerRadius*d.innerRadius {
		return false
	}

	// Test disk $\phi$ value against $\phimax$
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	if phi > d.phiMax {
		return false
	}
	return true
}

func (d *Disk) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (d *Disk) Area() float64 {
	return d.phiMax * 0.5 * (d.radius*d.radius - d.innerRadius*d.innerRadius)
}

func (d *Disk) Sample(u1, u2 float64) (*Point, *Normal) {
	var p Point
	p.X, p.Y = ConcentricSampleDisk(u1, u2)
	p.X *= d.radius
	p.Y *= d.radius
	p.Z = d.height
	Ns := NormalizeNormal(NormalTransform(d.objectToWorld, CreateNormal(0.0, 0.0, 1.0)))
	if d.reverseOrientation {
		Ns = Ns.Negate()
	}
	return PointTransform(d.objectToWorld, &p), Ns
}

func (d *Disk) Pdf(pshape *Point) float64 {
	return 1.0 / d.Area()
}
func (d *Disk) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return d.Sample(u1, u2)
}

func (d *Disk) Pdf2(p *Point, wi *Vector) float64 {
	return ShapePdf(d, p, wi)
}

func (d *Disk) ObjectToWorld() *Transform {
	return d.objectToWorld
}
func (d *Disk) WorldToObject() *Transform {
	return d.worldToObject
}

func (d *Disk) ReverseOrientation() bool {
	return d.reverseOrientation
}

func (d *Disk) TransformSwapsHandedness() bool {
	return d.transformSwapsHandedness
}

func (d *Disk) ShapeId() uint32 {
	return d.shapeId
}
