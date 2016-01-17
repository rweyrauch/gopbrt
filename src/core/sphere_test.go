package pbrt

import (
	"fmt"
	"testing"
)

func TestSphereShape(t *testing.T) {

	ident := CreateIdentityMatrix4x4()
	fmt.Printf("Ident: %v\n", ident)

	o2w, _ := CreateTransform(ident)
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
