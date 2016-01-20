package core

import (
	"math"
	)

type Hyperboloid struct {
	ShapeData
	p1, p2 Point
	zmin, zmax, phiMax, rmax, a, c float64
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
	
     radius1 := math.Sqrt(h.p1.x*h.p1.x + h.p1.y*h.p1.y)
     radius2 := math.Sqrt(h.p2.x*h.p2.x + h.p2.y*h.p2.y)
    h.rmax = math.Max(radius1, radius2)
    h.zmin = math.Min(h.p1.z, h.p2.z)
    h.zmax = math.Max(h.p1.z, h.p2.z)
	
	return h
}

func (h *Hyperboloid) ObjectBound() *BBox {
    return &BBox{Point{-h.rmax, -h.rmax, h.zmin}, Point{h.rmax,  h.rmax, h.zmax}}
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

func (h *Hyperboloid) Intersect(r *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	// Transform _Ray_ to object space
	ray := RayTransform(h.worldToObject, r)

    // Compute quadratic hyperboloid coefficients
     A := h.a*ray.dir.x*ray.dir.x + h.a*ray.dir.y*ray.dir.y - h.c*ray.dir.z*ray.dir.z
     B := 2.0 * (h.a*ray.dir.x*ray.origin.x + h.a*ray.dir.y*ray.origin.y - h.c*ray.dir.z*ray.origin.z)
     C := h.a*ray.origin.x*ray.origin.x + h.a*ray.origin.y*ray.origin.y - h.c*ray.origin.z*ray.origin.z - 1

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

    // Compute hyperboloid inverse mapping
    phit := ray.PointAt(thit)
    v := (phit.z - h.p1.z)/(h.p2.z - h.p1.z)
    pr := h.p1.Scale(1.0-v).Sub(h.p2.Scale(-v)) // using Sub(-v) rather than Add(v)
    phi := math.Atan2(pr.x*phit.y - phit.x*pr.y, phit.x*pr.x + phit.y*pr.y)
    if phi < 0.0 {
        phi += 2.0*math.Pi
	}
    
    // Test hyperboloid intersection against clipping parameters
    if phit.z < h.zmin || phit.z > h.zmax || phi > h.phiMax {
        if thit == t1 { return false, 0.0, 0.0, nil }
        thit = t1
        if t1 > ray.maxt { return false, 0.0, 0.0, nil }
        // Compute hyperboloid inverse mapping
        phit = ray.PointAt(thit)
        v = (phit.z - h.p1.z)/(h.p2.z - h.p1.z)
        pr := h.p1.Scale(1.0-v).Sub(h.p2.Scale(-v))
        phi = math.Atan2(pr.x*phit.y - phit.x*pr.y, phit.x*pr.x + phit.y*pr.y)
        if phi < 0 {
            phi += 2*math.Pi
            }
        if phit.z < h.zmin || phit.z > h.zmax || phi > h.phiMax {
            return false, 0.0, 0.0, nil
       }
    }

    // Compute parametric representation of hyperboloid hit
    u := phi / h.phiMax

    // Compute hyperboloid $\dpdu$ and $\dpdv$
    cosphi, sinphi := math.Cos(phi), math.Sin(phi)
    dpdu := CreateVector(-h.phiMax * phit.y, h.phiMax * phit.x, 0.0)
    dpdv := CreateVector((h.p2.x-h.p1.x) * cosphi - (h.p2.y-h.p1.y) * sinphi,
        (h.p2.x-h.p1.x) * sinphi + (h.p2.y-h.p1.y) * cosphi,
        h.p2.z-h.p1.z)

    // Compute hyperboloid $\dndu$ and $\dndv$
     d2Pduu := CreateVector(phit.x, phit.y, 0).Scale(-h.phiMax * h.phiMax)
     d2Pduv := CreateVector(-dpdv.y, dpdv.x, 0).Scale(h.phiMax)
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

func (h *Hyperboloid) IntersectP(r *Ray) bool {
	// Transform _Ray_ to object space
	ray := RayTransform(h.worldToObject, r)

    // Compute quadratic hyperboloid coefficients
     A := h.a*ray.dir.x*ray.dir.x + h.a*ray.dir.y*ray.dir.y - h.c*ray.dir.z*ray.dir.z
     B := 2.0 * (h.a*ray.dir.x*ray.origin.x + h.a*ray.dir.y*ray.origin.y - h.c*ray.dir.z*ray.origin.z)
     C := h.a*ray.origin.x*ray.origin.x + h.a*ray.origin.y*ray.origin.y - h.c*ray.origin.z*ray.origin.z - 1

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

    // Compute hyperboloid inverse mapping
    phit := ray.PointAt(thit)
    v := (phit.z - h.p1.z)/(h.p2.z - h.p1.z)
    pr := h.p1.Scale(1.0-v).Sub(h.p2.Scale(-v)) // using Sub(-v) rather than Add(v)
    phi := math.Atan2(pr.x*phit.y - phit.x*pr.y, phit.x*pr.x + phit.y*pr.y)
    if phi < 0.0 {
        phi += 2.0*math.Pi
	}
    
    // Test hyperboloid intersection against clipping parameters
    if phit.z < h.zmin || phit.z > h.zmax || phi > h.phiMax {
        if thit == t1 { return false }
        thit = t1
        if t1 > ray.maxt { return false }
        // Compute hyperboloid inverse mapping
        phit = ray.PointAt(thit)
        v = (phit.z - h.p1.z)/(h.p2.z - h.p1.z)
        pr := h.p1.Scale(1.0-v).Sub(h.p2.Scale(-v))
        phi = math.Atan2(pr.x*phit.y - phit.x*pr.y, phit.x*pr.x + phit.y*pr.y)
        if phi < 0 {
            phi += 2*math.Pi
            }
        if phit.z < h.zmin || phit.z > h.zmax || phi > h.phiMax {
            return false
       }
    }
    return true
}

func (h *Hyperboloid) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return dg
}

func (h *Hyperboloid) Area() float64 {
	SQR := func(a float64) float64 { return a*a }
	QUAD := func(a float64) float64 { return SQR(a)*SQR(a) }

    return h.phiMax/6.0 *
       (2.0*QUAD(h.p1.x) - 2.0*h.p1.x*h.p1.x*h.p1.x*h.p2.x +
         2.0*QUAD(h.p2.x) +
            2.0*(h.p1.y*h.p1.y + h.p1.y*h.p2.y + h.p2.y*h.p2.y)*
            (SQR(h.p1.y - h.p2.y) + SQR(h.p1.z - h.p2.z)) +
           h.p2.x*h.p2.x*(5.0*h.p1.y*h.p1.y + 2.0*h.p1.y*h.p2.y -
              4.0*h.p2.y*h.p2.y + 2.0*SQR(h.p1.z - h.p2.z)) +
           h.p1.x*h.p1.x*(-4.0*h.p1.y*h.p1.y + 2.0*h.p1.y*h.p2.y +
              5.0*h.p2.y*h.p2.y + 2.0*SQR(h.p1.z - h.p2.z)) -
           2.0*h.p1.x*h.p2.x*(h.p2.x*h.p2.x - h.p1.y*h.p1.y +
              5.0*h.p1.y*h.p2.y - h.p2.y*h.p2.y - h.p1.z*h.p1.z +
              2.0*h.p1.z*h.p2.z - h.p2.z*h.p2.z))
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

func CreateHyperboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Hyperboloid {
    p1 := params.FindPointParam("p1", Point{0,0,0})
    p2 := params.FindPointParam("p2", Point{1,1,1})
    phimax := params.FindFloatParam("phimax", 360)
    return NewHyperboloid(o2w, w2o, reverseOrientation, p1, p2, phimax)
}
