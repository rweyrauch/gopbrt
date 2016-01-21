package core

import (
	"math"
)

type Paraboloid struct {
	ShapeData
	radius, zmin, zmax, phiMax float64
}

func NewParaboloid(o2w, w2o *Transform, ro bool, rad, zmin, zmax, pmax float64) *Paraboloid {
	p := new(Paraboloid)
	p.objectToWorld = o2w
	p.worldToObject = w2o
	p.reverseOrientation = ro
	p.transformSwapsHandedness = SwapsHandednessTransform(p.objectToWorld)
	p.shapeId = GenerateShapeId()
	p.radius = rad
	p.zmin = zmin
	p.zmax = zmax
	p.phiMax = Radians(Clamp(pmax, 0.0, 360.0))

	return p
}

func (p *Paraboloid) ObjectBound() *BBox {
	return &BBox{Point{-p.radius, -p.radius, p.zmin}, Point{p.radius, p.radius, p.zmax}}
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
	k := p.zmax / (p.radius * p.radius)
	A := k * (ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y)
	B := 2*k*(ray.dir.x*ray.origin.x+ray.dir.y*ray.origin.y) - ray.dir.z
	C := k*(ray.origin.x*ray.origin.x+ray.origin.y*ray.origin.y) - ray.origin.z

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
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test paraboloid intersection against clipping parameters
	if phit.z < p.zmin || phit.z > p.zmax || phi > p.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.maxt {
			return false, 0.0, 0.0, nil
		}
		// Compute paraboloid inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < p.zmin || phit.z > p.zmax || phi > p.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of paraboloid hit
	u := phi / p.phiMax
	v := (phit.z - p.zmin) / (p.zmax - p.zmin)

	// Compute parabaloid $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-p.phiMax*phit.y, p.phiMax*phit.x, 0.0)
	dpdv := CreateVector(phit.x/(2.0*phit.z), phit.y/(2.0*phit.z), 1.0).Scale(p.zmax - p.zmin)

	// Compute parabaloid $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.x, phit.y, 0).Scale(-p.phiMax * p.phiMax)
	d2Pduv := CreateVector(-phit.y/(2.0*phit.z), phit.x/(2.0*phit.z), 0).Scale((p.zmax - p.zmin) * p.phiMax)
	d2Pdvv := CreateVector(phit.x/(4.0*phit.z*phit.z), phit.y/(4.0*phit.z*phit.z), 0.0).Scale(-(p.zmax - p.zmin) * (p.zmax - p.zmin))

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
	k := p.zmax / (p.radius * p.radius)
	A := k * (ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y)
	B := 2*k*(ray.dir.x*ray.origin.x+ray.dir.y*ray.origin.y) - ray.dir.z
	C := k*(ray.origin.x*ray.origin.x+ray.origin.y*ray.origin.y) - ray.origin.z

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
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test paraboloid intersection against clipping parameters
	if phit.z < p.zmin || phit.z > p.zmax || phi > p.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.maxt {
			return false
		}
		// Compute paraboloid inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < p.zmin || phit.z > p.zmax || phi > p.phiMax {
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
	k := 4 * p.zmax / radius2
	return (radius2 * radius2 * p.phiMax / (12.0 * p.zmax * p.zmax)) *
		(math.Pow(k*p.zmax+1, 1.5) - math.Pow(k*p.zmin+1, 1.5))
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
