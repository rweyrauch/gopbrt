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

type DifferentialGeometry struct {
	p                      *Point
	nn                     *Normal
	u, v                   float64
	shape                  Shape
	dpdu, dpdv             *Vector
	dndu, dndv             *Normal
	dpdx, dpdy             *Vector
	dudx, dvdx, dudy, dvdy float64
}

func CreateDiffGeometry(p *Point, dpdu, dpdv *Vector, dndu, dndv *Normal, uu, vv float64, sh Shape) *DifferentialGeometry {
	dg := &DifferentialGeometry{p: p, u: uu, v: vv, shape: sh, dpdu: dpdu, dpdv: dpdv, dndu: dndu, dndv: dndv}
	dg.nn = CreateNormalFromVector(NormalizeVector(CrossVector(dpdu, dpdv)))
	dg.dudx = 0.0
	dg.dudy = 0.0
	dg.dvdx = 0.0
	dg.dvdy = 0.0

	// Adjust normal based on orientation and handedness
	if dg.shape != nil && Xor(dg.shape.ReverseOrientation(), dg.shape.TransformSwapsHandedness()) {
		dg.nn = dg.nn.Negate()
	}
	return dg
}

func (dg *DifferentialGeometry) ComputeDifferentials(ray *RayDifferential) {
	dg.dudx, dg.dvdx = 0.0, 0.0
	dg.dudy, dg.dvdy = 0.0, 0.0
	dg.dpdx, dg.dpdy = &Vector{0, 0, 0}, &Vector{0, 0, 0}

	if ray.HasDifferentials {
		// Estimate screen space change in $\pt{}$ and $(u,v)$

		// Compute auxiliary intersection points with plane
		d := -DotNormalVector(dg.nn, &Vector{dg.p.X, dg.p.Y, dg.p.Z})
		rxv := &Vector{ray.RxOrigin.X, ray.RxOrigin.Y, ray.RxOrigin.Z}
		tx := -(DotNormalVector(dg.nn, rxv) + d) / DotNormalVector(dg.nn, &ray.RxDirection)
		if math.IsNaN(tx) {
			goto fail
		}
		px := ray.RxOrigin.Add(ray.RxDirection.Scale(tx))
		ryv := &Vector{ray.RyOrigin.X, ray.RyOrigin.Y, ray.RyOrigin.Z}
		ty := -(DotNormalVector(dg.nn, ryv) + d) / DotNormalVector(dg.nn, &ray.RyDirection)
		if math.IsNaN(ty) {
			goto fail
		}
		py := ray.RyOrigin.Add(ray.RyDirection.Scale(ty))
		dg.dpdx = px.Sub(dg.p)
		dg.dpdy = py.Sub(dg.p)

		// Compute $(u,v)$ offsets at auxiliary points

		// Initialize _A_, _Bx_, and _By_ matrices for offset computation
		var A [2][2]float64
		var Bx, By [2]float64
		var axes [2]int
		if math.Abs(dg.nn.X) > math.Abs(dg.nn.Y) && math.Abs(dg.nn.X) > math.Abs(dg.nn.Z) {
			axes[0] = 1
			axes[1] = 2
		} else if math.Abs(dg.nn.Y) > math.Abs(dg.nn.Z) {
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
	}
fail:
}

func (dg *DifferentialGeometry) String() string {
	return fmt.Sprintf("dg[p: %v n: %v u,v: %f,%f]", dg.p, dg.nn, dg.u, dg.v)
}
