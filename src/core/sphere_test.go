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

func TestSphereShape(t *testing.T) {

	ident := NewIdentityMatrix4x4()
	fmt.Printf("Ident: %v\n", ident)

	o2w, _ := NewTransform(ident)
	w2o := InverseTransform(o2w)
	sphere := CreateSphere(o2w, w2o, false, 5.0, -5.0, 5.0, 360.0)

	fmt.Printf("Sphere: %v\n", sphere)

	ray := CreateRay(&Point{10.0, 0.0, 0.0}, &Vector{-1.0, 0.0, 0.0}, 0.0, 100.0, 0.0, 0)
	fmt.Printf("Ray: %v\n", ray)

	if sphere.IntersectP(ray) {
		fmt.Printf("Ray hit sphere.\n")
	} else {
		t.Fail()
	}

	if ok, thit, _, _ := sphere.Intersect(ray); ok {
		fmt.Printf("Ray hit sphere at: %f\n", thit)
	} else {
		t.Fail()
	}

	ray = CreateRay(&Point{10.0, 2.0, 0.0}, &Vector{-1.0, 0.0, 0.0}, 0.0, 100.0, 0.0, 0)
	if sphere.IntersectP(ray) {
		fmt.Printf("Ray hit sphere.\n")
	} else {
		t.Fail()
	}

	if ok, thit, _, _ := sphere.Intersect(ray); ok {
		fmt.Printf("Ray hit sphere at: %f\n", thit)
	} else {
		t.Fail()
	}
}
