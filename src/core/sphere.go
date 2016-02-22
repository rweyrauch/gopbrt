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

type Sphere struct {
	ShapeData
	radius, phiMax, Zmin, Zmax, thetaMin, thetaMax float64
}

func CreateSphere(o2w, w2o *Transform, ro bool, rad, z0, z1, pm float64) *Sphere {
	s := new(Sphere)
	s.objectToWorld = o2w
	s.worldToObject = w2o
	s.reverseOrientation = ro
	s.transformSwapsHandedness = SwapsHandednessTransform(s.objectToWorld)
	s.shapeId = GenerateShapeId()
	s.radius = rad
	s.Zmin = Clamp(math.Min(z0, z1), -s.radius, s.radius)
	s.Zmax = Clamp(math.Max(z0, z1), -s.radius, s.radius)
	s.thetaMin = math.Acos(Clamp(s.Zmin/s.radius, -1.0, 1.0))
	s.thetaMax = math.Acos(Clamp(s.Zmax/s.radius, -1.0, 1.0))
	s.phiMax = Radians(Clamp(pm, 0.0, 360.0))

	return s
}

func (s *Sphere) String() string {
	return fmt.Sprintf("sphere[r: %f zmin: %f zmax: %f phimax: %f obj2world: %v]", s.radius, s.Zmin, s.Zmax, Degrees(s.phiMax), s.objectToWorld)
}

func (s *Sphere) ObjectBound() *BBox {
	return &BBox{Point{-s.radius, -s.radius, s.Zmin}, Point{s.radius, s.radius, s.Zmax}}
}

func (s *Sphere) WorldBound() *BBox {
	return BBoxTransform(s.objectToWorld, s.ObjectBound())
}

func (*Sphere) CanIntersect() bool {
	return true
}

func (*Sphere) Refine(refined []Shape) []Shape {
	return refined
}

func (s *Sphere) Intersect(r RayBase) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := r.Transform(s.worldToObject)

	// Compute quadratic sphere coefficients
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y + ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y + ray.Dir().Z*ray.Origin().Z)
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y + ray.Origin().Z*ray.Origin().Z - s.radius*s.radius

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

	// Compute sphere hit position and $\phi$
	phit := ray.PointAt(thit)
	if phit.X == 0.0 && phit.Y == 0.0 {
		phit.X = 1.0e-5 * s.radius
	}
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.Zmin > -s.radius && phit.Z < s.Zmin) ||
		(s.Zmax < s.radius && phit.Z > s.Zmax) || phi > s.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		if t1 > ray.Maxt() {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		// Compute sphere hit position and $\phi$
		phit = ray.PointAt(thit)
		if phit.X == 0.0 && phit.Y == 0.0 {
			phit.X = 1.0e-5 * s.radius
		}
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if (s.Zmin > -s.radius && phit.Z < s.Zmin) ||
			(s.Zmax < s.radius && phit.Z > s.Zmax) || phi > s.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of sphere hit
	u := phi / s.phiMax
	theta := math.Acos(Clamp(phit.Z/s.radius, -1.0, 1.0))
	v := (theta - s.thetaMin) / (s.thetaMax - s.thetaMin)

	// Compute sphere $\dpdu$ and $\dpdv$
	zradius := math.Sqrt(phit.X*phit.X + phit.Y*phit.Y)
	invzradius := 1.0 / zradius
	cosphi := phit.X * invzradius
	sinphi := phit.Y * invzradius
	dpdu := CreateVector(-s.phiMax*phit.Y, s.phiMax*phit.X, 0.0)
	dpdv := CreateVector(phit.Z*cosphi, phit.Z*sinphi, -s.radius*math.Sin(theta)).Scale(s.thetaMax - s.thetaMin)

	// Compute sphere $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.X, phit.Y, 0).Scale(-s.phiMax * s.phiMax)
	d2Pduv := CreateVector(-sinphi, cosphi, 0.0).Scale((s.thetaMax - s.thetaMin) * phit.Z * s.phiMax)
	d2Pdvv := CreateVector(phit.X, phit.Y, phit.Z).Scale(-(s.thetaMax - s.thetaMin) * (s.thetaMax - s.thetaMin))

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
	dg = CreateDiffGeometry(PointTransform(s.objectToWorld, phit), VectorTransform(s.objectToWorld, dpdu), VectorTransform(s.objectToWorld, dpdv),
		NormalTransform(s.objectToWorld, dndu), NormalTransform(s.objectToWorld, dndv), u, v, s)

	// Update _tHit_ for quadric intersection
	tHit = thit

	// Compute _rayEpsilon_ for quadric intersection
	rayEpsilon = 5.0e-4 * tHit

	return true, tHit, rayEpsilon, dg
}

func (s *Sphere) IntersectP(r RayBase) bool {
	// Transform _Ray_ to object space
	ray := r.Transform(s.worldToObject)

	// Compute quadratic sphere coefficients
	A := ray.Dir().X*ray.Dir().X + ray.Dir().Y*ray.Dir().Y + ray.Dir().Z*ray.Dir().Z
	B := 2.0 * (ray.Dir().X*ray.Origin().X + ray.Dir().Y*ray.Origin().Y + ray.Dir().Z*ray.Origin().Z)
	C := ray.Origin().X*ray.Origin().X + ray.Origin().Y*ray.Origin().Y + ray.Origin().Z*ray.Origin().Z - s.radius*s.radius

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

	// Compute sphere hit position and $\phi$
	phit := ray.PointAt(thit)
	if phit.X == 0.0 && phit.Y == 0.0 {
		phit.X = 1.0e-5 * s.radius
	}
	phi := math.Atan2(phit.Y, phit.X)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.Zmin > -s.radius && phit.Z < s.Zmin) ||
		(s.Zmax < s.radius && phit.Z > s.Zmax) || phi > s.phiMax {
		if thit == t1 {
			return false
		}
		if t1 > ray.Maxt() {
			return false
		}
		thit = t1
		// Compute sphere hit position and $\phi$
		phit = ray.PointAt(thit)
		if phit.X == 0.0 && phit.Y == 0.0 {
			phit.X = 1.0e-5 * s.radius
		}
		phi = math.Atan2(phit.Y, phit.X)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if (s.Zmin > -s.radius && phit.Z < s.Zmin) ||
			(s.Zmax < s.radius && phit.Z > s.Zmax) || phi > s.phiMax {
			return false
		}
	}

	return true
}

func (s *Sphere) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (s *Sphere) Area() float64 {
	return s.phiMax * s.radius * (s.Zmax - s.Zmin)
}

func (s *Sphere) Sample(u1, u2 float64) (*Point, *Normal) {
	p := CreatePoint(0, 0, 0).Add(UniformSampleSphere(u1, u2).Scale(s.radius))
	ns := NormalizeNormal(NormalTransform(s.objectToWorld, &Normal{p.X, p.Y, p.Z}))
	if s.ReverseOrientation() {
		ns = ns.Negate()
	}
	return PointTransform(s.objectToWorld, p), ns
}

func (s *Sphere) Pdf(pshape *Point) float64 {
	return 1.0 / s.Area()
}

func (s *Sphere) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	// Compute coordinate system for sphere sampling
	Pcenter := PointTransform(s.objectToWorld, &Point{0, 0, 0})
	wc := NormalizeVector(Pcenter.Sub(p))
	wcX, wcY := CoordinateSystem(wc)

	// Sample uniformly on sphere if $\pt{}$ is inside it
	if DistanceSquaredPoint(p, Pcenter)-s.radius*s.radius < 1.0e-4 {
		return s.Sample(u1, u2)
	}

	// Sample sphere uniformly inside subtended cone
	sinThetaMax2 := s.radius * s.radius / DistanceSquaredPoint(p, Pcenter)
	cosThetaMax := math.Sqrt(math.Max(0.0, 1.0-sinThetaMax2))

	r := CreateRay(p, UniformSampleConeVector(u1, u2, cosThetaMax, wcX, wcY, wc), 1.0e-3, INFINITY, 0.0, 0)
	var thit float64
	var ok bool
	if ok, thit, _, _ = s.Intersect(r); !ok {
		thit = DotVector(Pcenter.Sub(p), NormalizeVector(r.Dir()))
	}
	ps := r.PointAt(thit)
	ns := CreateNormalFromVector(NormalizeVector(ps.Sub(Pcenter)))
	if s.ReverseOrientation() {
		ns = ns.Negate()
	}
	return ps, ns
}

func (s *Sphere) Pdf2(p *Point, wi *Vector) float64 {
	Pcenter := PointTransform(s.objectToWorld, &Point{0, 0, 0})
	// Return uniform weight if point inside sphere
	if DistanceSquaredPoint(p, Pcenter)-s.radius*s.radius < 1.0e-4 {
		return ShapePdf(s, p, wi)
	}
	// Compute general sphere weight
	sinThetaMax2 := s.radius * s.radius / DistanceSquaredPoint(p, Pcenter)
	cosThetaMax := math.Sqrt(math.Max(0.0, 1.0-sinThetaMax2))
	return UniformConePdf(cosThetaMax)
}

func (s *Sphere) ObjectToWorld() *Transform {
	return s.objectToWorld
}

func (s *Sphere) WorldToObject() *Transform {
	return s.worldToObject
}

func (s *Sphere) ReverseOrientation() bool {
	return s.reverseOrientation
}

func (s *Sphere) TransformSwapsHandedness() bool {
	return s.transformSwapsHandedness
}

func (s *Sphere) ShapeId() uint32 {
	return s.shapeId
}
