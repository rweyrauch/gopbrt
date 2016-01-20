package core

import (
	"fmt"
	"math"
)

type Cylinder struct {
	ShapeData
	radius, zmin, zmax, phiMax float64
}

func NewCylinder(o2w, w2o *Transform, ro bool, rad, zmin, zmax, pmax float64) *Cylinder {
	c := new(Cylinder)
	c.objectToWorld = o2w
	c.worldToObject = w2o
	c.reverseOrientation = ro
	c.transformSwapsHandedness = SwapsHandednessTransform(c.objectToWorld)
	c.shapeId = GenerateShapeId()
	c.radius = rad
	c.zmin = zmin
	c.zmax = zmax
	c.phiMax = Radians(Clamp(pmax, 0.0, 360.0))

	return c
}

func (c *Cylinder) String() string {
	return fmt.Sprintf("cylinder[r: %f zmin: %f zmax: %f phimax: %f obj2world: %v]",
		c.radius, c.zmin, c.zmax, Degrees(c.phiMax), c.objectToWorld)
}

func (c *Cylinder) ObjectBound() *BBox {
	return &BBox{Point{-c.radius, -c.radius, c.zmin}, Point{c.radius, c.radius, c.zmax}}
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

func (c *Cylinder) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(c.worldToObject, r)

	// Compute quadratic cylinder coefficients
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y)
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y - c.radius*c.radius

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

	// Compute cylinder hit point and $\phi$
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cylinder intersection against clipping parameters
	if phit.z < c.zmin || phit.z > c.zmax || phi > c.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.maxt {
			return false, 0.0, 0.0, nil
		}
		// Compute cylinder hit point and $\phi$
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < c.zmin || phit.z > c.zmax || phi > c.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of cylinder hit
	u := phi / c.phiMax
	v := (phit.z - c.zmin) / (c.zmax - c.zmin)

	// Compute cylinder $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-c.phiMax*phit.y, c.phiMax*phit.x, 0.0)
	dpdv := CreateVector(0.0, 0.0, c.zmax-c.zmin)

	// Compute cylinder $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.x, phit.y, 0).Scale(-c.phiMax * c.phiMax)
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

func (c *Cylinder) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(c.worldToObject, r)

	// Compute quadratic cylinder coefficients
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y)
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y - c.radius*c.radius

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

	// Compute cylinder hit point and $\phi$
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cylinder intersection against clipping parameters
	if phit.z < c.zmin || phit.z > c.zmax || phi > c.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.maxt {
			return false
		}
		// Compute cylinder hit point and $\phi$
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < c.zmin || phit.z > c.zmax || phi > c.phiMax {
			return false
		}
	}
	return true
}

func (c *Cylinder) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (c *Cylinder) Area() float64 {
	return (c.zmax - c.zmin) * c.phiMax * c.radius
}

func (c *Cylinder) Sample(u1, u2 float64) (*Point, *Normal) {
	z := Lerp(u1, c.zmin, c.zmax)
	t := u2 * c.phiMax
	p := CreatePoint(c.radius*math.Cos(t), c.radius*math.Sin(t), z)
	Ns := NormalizeNormal(NormalTransform(c.objectToWorld, CreateNormal(p.x, p.y, 0.0)))
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

func CreateCylinderShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cylinder {
	radius := params.FindFloatParam("radius", 1.0)
	zmin := params.FindFloatParam("zmin", -1.0)
	zmax := params.FindFloatParam("zmax", 1.0)
	phimax := params.FindFloatParam("phimax", 360)
	return NewCylinder(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}
