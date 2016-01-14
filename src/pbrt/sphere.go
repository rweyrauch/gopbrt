package pbrt

import (
	"math"
	"strings"
)

type Sphere struct {
	ShapeData
	radius, phiMax, zmin, zmax, thetaMin, thetaMax float64
}

func CreateSphere(o2w, w2o *Transform, ro bool, rad, z0, z1, pm float64) *Sphere {
	s := new(Sphere)
	s.objectToWorld = o2w
	s.worldToObject = w2o
	s.reverseOrientation = ro
	s.transformSwapsHandedness = false
	s.shapeId = GenerateShapeId()
	s.radius = rad
	s.zmin = Clamp(math.Min(z0, z1), -s.radius, s.radius)
	s.zmax = Clamp(math.Max(z0, z1), -s.radius, s.radius)
	s.thetaMin = math.Acos(Clamp(s.zmin/s.radius, -1.0, 1.0))
	s.thetaMax = math.Acos(Clamp(s.zmax/s.radius, -1.0, 1.0))
	s.phiMax = Radians(Clamp(pm, 0.0, 360.0))

	return s
}

func (s *Sphere) ObjectBound() *BBox {
	return &BBox{Point{-s.radius, -s.radius, s.zmin}, Point{s.radius, s.radius, s.zmax}}
}

func (s *Sphere) WorldBound() *BBox {
	return BBoxTransform(s.objectToWorld, s.ObjectBound())
}
func (*Sphere) CanIntersect() bool {
	return true
}
func (*Sphere) Refine() (refined []*Shape) {
	return nil
}
func (s *Sphere) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(s.worldToObject, r)

	// Compute quadratic sphere coefficients
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y + ray.dir.z*ray.dir.z
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y + ray.dir.z*ray.origin.z)
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y + ray.origin.z*ray.origin.z - s.radius*s.radius

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); ok {
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

	// Compute sphere hit position and $\phi$
	phit := ray.PointAt(thit)
	if phit.x == 0.0 && phit.y == 0.0 {
		phit.x = 1.0e-5 * s.radius
	}
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zmin > -s.radius && phit.z < s.zmin) ||
		(s.zmax < s.radius && phit.z > s.zmax) || phi > s.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		if t1 > ray.maxt {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		// Compute sphere hit position and $\phi$
		phit = ray.PointAt(thit)
		if phit.x == 0.0 && phit.y == 0.0 {
			phit.x = 1.0e-5 * s.radius
		}
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if (s.zmin > -s.radius && phit.z < s.zmin) ||
			(s.zmax < s.radius && phit.z > s.zmax) || phi > s.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of sphere hit
	u := phi / s.phiMax
	theta := math.Acos(Clamp(phit.z/s.radius, -1.0, 1.0))
	v := (theta - s.thetaMin) / (s.thetaMax - s.thetaMin)

	// Compute sphere $\dpdu$ and $\dpdv$
	zradius := math.Sqrt(phit.x*phit.x + phit.y*phit.y)
	invzradius := 1.0 / zradius
	cosphi := phit.x * invzradius
	sinphi := phit.y * invzradius
	dpdu := CreateVector(-s.phiMax*phit.y, s.phiMax*phit.x, 0.0)
	dpdv := CreateVector(phit.z*cosphi, phit.z*sinphi, -s.radius*math.Sin(theta)).Scale(s.thetaMax - s.thetaMin)

	// Compute sphere $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.x, phit.y, 0).Scale(-s.phiMax * s.phiMax)
	d2Pduv := CreateVector(-sinphi, cosphi, 0.0).Scale((s.thetaMax - s.thetaMin) * phit.z * s.phiMax)
	d2Pdvv := CreateVector(phit.x, phit.y, phit.z).Scale(-(s.thetaMax - s.thetaMin) * (s.thetaMax - s.thetaMin))

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

func (s *Sphere) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(s.worldToObject, r)

	// Compute quadratic sphere coefficients
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y + ray.dir.z*ray.dir.z
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y + ray.dir.z*ray.origin.z)
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y + ray.origin.z*ray.origin.z - s.radius*s.radius

	// Solve quadratic equation for _t_ values
	var t0, t1 float64
	var ok bool
	if ok, t0, t1 = Quadratic(A, B, C); ok {
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

	// Compute sphere hit position and $\phi$
	phit := ray.PointAt(thit)
	if phit.x == 0.0 && phit.y == 0.0 {
		phit.x = 1.0e-5 * s.radius
	}
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zmin > -s.radius && phit.z < s.zmin) ||
		(s.zmax < s.radius && phit.z > s.zmax) || phi > s.phiMax {
		if thit == t1 {
			return false
		}
		if t1 > ray.maxt {
			return false
		}
		thit = t1
		// Compute sphere hit position and $\phi$
		phit = ray.PointAt(thit)
		if phit.x == 0.0 && phit.y == 0.0 {
			phit.x = 1.0e-5 * s.radius
		}
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if (s.zmin > -s.radius && phit.z < s.zmin) ||
			(s.zmax < s.radius && phit.z > s.zmax) || phi > s.phiMax {
			return false
		}
	}
	return true
}

func (s *Sphere) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (s *Sphere) Area() float64 {
	return s.phiMax * s.radius * (s.zmax - s.zmin)
}
func (s *Sphere) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (s *Sphere) Pdf(pshape *Point) float64 {
	return 0.0
}
func (s *Sphere) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (s *Sphere) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
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

func CreateSphereShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Sphere {
	radius := 1.0
	zmin := -radius
	zmax := radius
	phimax := 360.0

	for i, t := range params.tokens {
		ps, ok := params.params[i].([]Object)
		if !ok {
			continue
		}
		if strings.Compare("radius", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				radius = v
			}
		} else if strings.Compare("zmin", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				zmin = v
			}
		} else if strings.Compare("zmax", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				zmax = v
			}
		} else if strings.Compare("phimax", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				phimax = v
			}
		}
	}
	return CreateSphere(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}
