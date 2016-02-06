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

type Paraboloid struct {
	ShapeData
	radius, Zmin, Zmax, phiMax float64
}

func NewParaboloid(o2w, w2o *Transform, ro bool, rad, zmin, zmax, pmax float64) *Paraboloid {
	p := new(Paraboloid)
	p.objectToWorld = o2w
	p.worldToObject = w2o
	p.reverseOrientation = ro
	p.transformSwapsHandedness = SwapsHandednessTransform(p.objectToWorld)
	p.shapeId = GenerateShapeId()
	p.radius = rad
	p.Zmin = zmin
	p.Zmax = zmax
	p.phiMax = Radians(Clamp(pmax, 0.0, 360.0))

	return p
}

func (p *Paraboloid) ObjectBound() *BBox {
	return &BBox{Point{-p.radius, -p.radius, p.Zmin}, Point{p.radius, p.radius, p.Zmax}}
}

func (p *Paraboloid) WorldBound() *BBox {
	return BBoxTransform(p.objectToWorld, p.ObjectBound())
}

func (p *Paraboloid) CanIntersect() bool {
	return true
}

func (p *Paraboloid) Refine(refined []Shape) []Shape {
	return refined
}

func (p *Paraboloid) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(p.worldToObject, r)

	// Compute quadratic paraboloid coefficients
	k := p.Zmax / (p.radius * p.radius)
	A := k * (ray.dir.X*ray.dir.X + ray.dir.Y*ray.dir.Y)
	B := 2*k*(ray.dir.X*ray.origin.X+ray.dir.Y*ray.origin.Y) - ray.dir.Z
	C := k*(ray.origin.X*ray.origin.X+ray.origin.Y*ray.origin.Y) - ray.origin.Z

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); !ok {
		return false, 0.0, 0.0, nil
	}

	// Compute intersection distance along ray
	if t0 > ray.maxt || t1 < ray.mint {
		return false, 0.0, 0.0, nil
	}

	thit := t0
	if t0 < ray.mint {
		thit = t1
		if thit > ray.maxt {
			return false, 0.0, 0.0, nil
		}
	}

	// Compute paraboloid inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test paraboloid intersection against clipping parameters
	if phit.Z < p.Zmin || phit.Z > p.Zmax || phi > p.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.maxt {
			return false, 0.0, 0.0, nil
		}
		// Compute paraboloid inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < p.Zmin || phit.Z > p.Zmax || phi > p.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of paraboloid hit
	u := phi / p.phiMax
	v := (phit.Z - p.Zmin) / (p.Zmax - p.Zmin)

	// Compute parabaloid $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-p.phiMax*phit.Y, p.phiMax*phit.X, 0.0)
	dpdv := CreateVector(phit.X/(2.0*phit.Z), phit.Y/(2.0*phit.Z), 1.0).Scale(p.Zmax - p.Zmin)

	// Compute parabaloid $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.X, phit.Y, 0).Scale(-p.phiMax * p.phiMax)
	d2Pduv := CreateVector(-phit.Y/(2.0*phit.Z), phit.X/(2.0*phit.Z), 0).Scale((p.Zmax - p.Zmin) * p.phiMax)
	d2Pdvv := CreateVector(phit.X/(4.0*phit.Z*phit.Z), phit.Y/(4.0*phit.Z*phit.Z), 0.0).Scale(-(p.Zmax - p.Zmin) * (p.Zmax - p.Zmin))

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
	dg = CreateDiffGeometry(PointTransform(p.objectToWorld, phit), VectorTransform(p.objectToWorld, dpdu), VectorTransform(p.objectToWorld, dpdv),
		NormalTransform(p.objectToWorld, dndu), NormalTransform(p.objectToWorld, dndv), u, v, p)

	// Update _tHit_ for quadric intersection
	tHit = thit

	// Compute _rayEpsilon_ for quadric intersection
	rayEpsilon = 5.0e-4 * tHit

	return true, tHit, rayEpsilon, dg
}

func (p *Paraboloid) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(p.worldToObject, r)

	// Compute quadratic paraboloid coefficients
	k := p.Zmax / (p.radius * p.radius)
	A := k * (ray.dir.X*ray.dir.X + ray.dir.Y*ray.dir.Y)
	B := 2*k*(ray.dir.X*ray.origin.X+ray.dir.Y*ray.origin.Y) - ray.dir.Z
	C := k*(ray.origin.X*ray.origin.X+ray.origin.Y*ray.origin.Y) - ray.origin.Z

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); !ok {
		return false
	}

	// Compute intersection distance along ray
	if t0 > ray.maxt || t1 < ray.mint {
		return false
	}

	thit := t0
	if t0 < ray.mint {
		thit = t1
		if thit > ray.maxt {
			return false
		}
	}

	// Compute paraboloid inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test paraboloid intersection against clipping parameters
	if phit.Z < p.Zmin || phit.Z > p.Zmax || phi > p.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.maxt {
			return false
		}
		// Compute paraboloid inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.Z < p.Zmin || phit.Z > p.Zmax || phi > p.phiMax {
			return false
		}
	}
	return true
}

func (p *Paraboloid) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (p *Paraboloid) Area() float64 {
	radius2 := p.radius * p.radius
	k := 4 * p.Zmax / radius2
	return (radius2 * radius2 * p.phiMax / (12.0 * p.Zmax * p.Zmax)) *
		(math.Pow(k*p.Zmax+1, 1.5) - math.Pow(k*p.Zmin+1, 1.5))
}

func (p *Paraboloid) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}

func (p *Paraboloid) Pdf(pshape *Point) float64 {
	return 1.0 / p.Area()
}

func (*Paraboloid) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}

func (p *Paraboloid) Pdf2(pp *Point, wi *Vector) float64 {
	return ShapePdf(p, pp, wi)
}

func (p *Paraboloid) ObjectToWorld() *Transform {
	return p.objectToWorld
}

func (p *Paraboloid) WorldToObject() *Transform {
	return p.worldToObject
}

func (p *Paraboloid) ReverseOrientation() bool {
	return p.reverseOrientation
}

func (p *Paraboloid) TransformSwapsHandedness() bool {
	return p.transformSwapsHandedness
}

func (p *Paraboloid) ShapeId() uint32 {
	return p.shapeId
}

func CreateParaboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Paraboloid {
	radius := params.FindFloatParam("radius", 1)
	zmin := params.FindFloatParam("zmin", 0)
	zmax := params.FindFloatParam("zmax", 1)
	phimax := params.FindFloatParam("phimax", 360)
	return NewParaboloid(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}
