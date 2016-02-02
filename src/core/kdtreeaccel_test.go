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
	"testing"
)

func TestKdTreeAceel(t *testing.T) {
	var options Options
	options.Verbose = true
	PbrtInit(&options)

	ident := NewIdentityMatrix4x4()

	o2w, _ := NewTransform(ident)
	w2o := InverseTransform(o2w)
	sphere := CreateSphere(o2w, w2o, false, 4.0, -4.0, 4.0, 360.0)

	vi := []int{0, 1, 2, 2, 3, 0}
	P := []Point{{-5, -5, -5}, {5, -5, -5}, {5, -5, 5}, {-5, -5, 5}}
	mesh := NewTriangleMesh(o2w, w2o, false, vi, P, nil, nil, nil, nil)

	primS := NewGeometricPrimitive(sphere, nil, nil)
	primM := NewGeometricPrimitive(mesh, nil, nil)

	prims := make([]Primitive, 0, 2)
	prims = append(prims, primS)
	prims = append(prims, primM)

	accel := NewKdTreeAccel(prims, 80, 1, 0.5, 1, -1)
	//accel := NewGridAccel(prims, true)

	ray := CreateRay(&Point{10.0, 0.0, 0.0}, &Vector{-1.0, 0.0, 0.0}, 0.0, 100.0, 0.0, 0)
	raydiff := CreateRayDifferentialFromRay(ray)
	fmt.Printf("Ray: %v\n", ray)

	if accel.IntersectP(ray) {
		fmt.Printf("Ray hit sphere.\n")
	} else {
		t.Fail()
	}

	if ok, _ := accel.Intersect(raydiff); ok {
		fmt.Printf("Ray hit sphere at: %f\n", raydiff.maxt)
	} else {
		t.Fail()
	}

	ray = CreateRay(&Point{10.0, 2.0, 0.0}, &Vector{-1.0, 0.0, 0.0}, 0.0, 100.0, 0.0, 0)
	if accel.IntersectP(ray) {
		fmt.Printf("Ray hit sphere.\n")
	} else {
		t.Fail()
	}

	raydiff = CreateRayDifferentialFromRay(ray)
	if ok, _ := accel.Intersect(raydiff); ok {
		fmt.Printf("Ray hit sphere at: %f\n", raydiff.maxt)
	} else {
		t.Fail()
	}

	ray = CreateRay(&Point{1.0, 10.0, -2.0}, &Vector{0.0, -1.0, 0.0}, 0.0, 100.0, 0.0, 0)
	fmt.Printf("Ray: %v\n", ray)
	raydiff = CreateRayDifferentialFromRay(ray)
	if accel.IntersectP(ray) {
		fmt.Printf("Ray hit plane.\n")
	} else {
		t.Fail()
	}

	if ok, _ := accel.Intersect(raydiff); ok {
		fmt.Printf("Ray hit plane at: %f\n", raydiff.maxt)
	} else {
		t.Fail()
	}

	PbrtCleanup()
}
