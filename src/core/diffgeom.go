package pbrt

import (
	"math"
)

type DifferentialGeometry struct {
    p Point
    nn Normal
    u, v float64
    shape Shape
    dpdu, dpdv Vector
    dndu, dndv Normal
    dpdx, dpdy Vector
    dudx, dvdx, dudy, dvdy float64
}

func CreateDiffGeometry(p Point, dpdu, dpdv Vector, dndu, dndv Normal, uu, vv float64, sh Shape) *DifferentialGeometry {
	dg := &DifferentialGeometry{p: p, u: uu, v: vv, shape: sh, dpdu: dpdu, dpdv: dpdv, dndu: dndu, dndv: dndv}
	dg.nn = *CreateNormalFromVector(NormalizeVector(CrossVector(&dpdu, &dpdv)))
	dg.dudx = 0.0
	dg.dudy = 0.0
	dg.dvdx = 0.0
	dg.dvdy = 0.0
	
    // Adjust normal based on orientation and handedness
    if dg.shape != nil && ((dg.shape.ReverseOrientation() && !dg.shape.TransformSwapsHandedness()) || (!dg.shape.ReverseOrientation() && dg.shape.TransformSwapsHandedness())) {
        dg.nn = *dg.nn.Negate()
	}
	return dg
}

func (dg *DifferentialGeometry) ComputeDifferentials(ray *RayDifferential) {
    if ray.hasDifferentials {
        // Estimate screen space change in $\pt{}$ and $(u,v)$

        // Compute auxiliary intersection points with plane
        d := -DotNormalVector(&dg.nn, &Vector{dg.p.x, dg.p.y, dg.p.z})
        rxv := Vector{ray.rxOrigin.x, ray.rxOrigin.y, ray.rxOrigin.z}
        tx := -(DotNormalVector(&dg.nn, &rxv) + d) / DotNormalVector(&dg.nn, &ray.rxDirection)
        if math.IsNaN(tx) { goto fail }
        px := ray.rxOrigin.Add(ray.rxDirection.Scale(tx))
        ryv := Vector{ray.ryOrigin.x, ray.ryOrigin.y, ray.ryOrigin.z}
        ty := -(DotNormalVector(&dg.nn, &ryv) + d) / DotNormalVector(&dg.nn, &ray.ryDirection)
        if math.IsNaN(ty) { goto fail }
        py := ray.ryOrigin.Add(ray.ryDirection.Scale(ty))
        dg.dpdx = *px.Sub(&dg.p)
        dg.dpdy = *py.Sub(&dg.p)

        // Compute $(u,v)$ offsets at auxiliary points

        // Initialize _A_, _Bx_, and _By_ matrices for offset computation
        var A [2][2]float64
        var Bx, By [2]float64
        var axes [2]int
        if math.Abs(dg.nn.x) > math.Abs(dg.nn.y) && math.Abs(dg.nn.x) > math.Abs(dg.nn.z) {
            axes[0] = 1 
            axes[1] = 2
        } else if math.Abs(dg.nn.y) > math.Abs(dg.nn.z) {
            axes[0] = 0
            axes[1] = 2
        } else {
            axes[0] = 0
            axes[1] = 1
        }

        // Initialize matrices for chosen projection plane
        A[0][0] = dg.dpdu.At(axes[0])
        A[0][1] = dg.dpdv.At(axes[0])
        A[1][0] = dg.dpdu.At(axes[1])
        A[1][1] = dg.dpdv.At(axes[1])
        Bx[0] = px.At(axes[0]) - dg.p.At(axes[0])
        Bx[1] = px.At(axes[1]) - dg.p.At(axes[1])
        By[0] = py.At(axes[0]) - dg.p.At(axes[0])
        By[1] = py.At(axes[1]) - dg.p.At(axes[1])
        var ok bool
        if ok, dg.dudx, dg.dvdx = SolveLinearSystem2x2(A, Bx); !ok {
            dg.dudx = 0.0
            dg.dvdx = 0.0
        }
        if ok, dg.dudy, dg.dvdy = SolveLinearSystem2x2(A, By); !ok {
            dg.dudy = 0.0 
            dg.dvdy = 0.0
        }
    } else {
fail:
        dg.dudx, dg.dvdx = 0.0, 0.0
        dg.dudy, dg.dvdy = 0.0, 0.0
        dg.dpdx, dg.dpdy = Vector{0,0,0}, Vector{0,0,0}
    }
	
}