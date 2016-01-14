package pbrt

import (
	"fmt"
	"math"
	"strings"
)

type Cone struct {
	ShapeData
	radius, height, phiMax float64
}

func CreateCone(o2w, w2o *Transform, ro bool, rad, h, pm float64) *Cone {
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

func (c *Cone) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(c.worldToObject, r)

	// Compute quadratic cone coefficients
	k := c.radius / c.height
	k = k * k
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y - k*ray.dir.z*ray.dir.z
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y - k*ray.dir.z*(ray.origin.z-c.height))
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y - k*(ray.origin.z-c.height)*(ray.origin.z-c.height)

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

	// Compute cone inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cone intersection against clipping parameters
	if phit.z < 0 || phit.z > c.height || phi > c.phiMax {
		if thit == t1 {
			return false, 0.0, 0.0, nil
		}
		thit = t1
		if t1 > ray.maxt {
			return false, 0.0, 0.0, nil
		}
		// Compute cone inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < 0 || phit.z > c.height || phi > c.phiMax {
			return false, 0.0, 0.0, nil
		}
	}

	// Find parametric representation of cone hit
	u := phi / c.phiMax
	v := phit.z / c.height

	// Compute cone $\dpdu$ and $\dpdv$
	dpdu := CreateVector(-c.phiMax*phit.y, c.phiMax*phit.x, 0.0)
	dpdv := CreateVector(-phit.x/(1.0-v), -phit.y/(1.0-v), c.height)

	// Compute cone $\dndu$ and $\dndv$
	d2Pduu := CreateVector(phit.x, phit.y, 0.0).Scale(-c.phiMax * c.phiMax)
	d2Pduv := CreateVector(phit.y, -phit.x, 0.0).Scale(c.phiMax / (1.0 - v))

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

func (c *Cone) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(c.worldToObject, r)

	// Compute quadratic cone coefficients
	k := c.radius / c.height
	k = k * k
	A := ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y - k*ray.dir.z*ray.dir.z
	B := 2.0 * (ray.dir.x*ray.origin.x + ray.dir.y*ray.origin.y - k*ray.dir.z*(ray.origin.z-c.height))
	C := ray.origin.x*ray.origin.x + ray.origin.y*ray.origin.y - k*(ray.origin.z-c.height)*(ray.origin.z-c.height)

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

	// Compute cone inverse mapping
	phit := ray.PointAt(thit)
	phi := math.Atan2(phit.y, phit.x)
	if phi < 0.0 {
		phi += 2.0 * math.Pi
	}

	// Test cone intersection against clipping parameters
	if phit.z < 0 || phit.z > c.height || phi > c.phiMax {
		if thit == t1 {
			return false
		}
		thit = t1
		if t1 > ray.maxt {
			return false
		}
		// Compute cone inverse mapping
		phit = ray.PointAt(thit)
		phi = math.Atan2(phit.y, phit.x)
		if phi < 0.0 {
			phi += 2.0 * math.Pi
		}
		if phit.z < 0 || phit.z > c.height || phi > c.phiMax {
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

func CreateConeShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cone {
	radius := 1.0
	height := radius
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
		} else if strings.Compare("height", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				height = v
			}
		} else if strings.Compare("phimax", t) == 0 {
			v, ok := extractFloatParam(ps)
			if ok {
				phimax = v
			}
		}
	}
	return CreateCone(o2w, w2o, reverseOrientation, radius, height, phimax)
}
