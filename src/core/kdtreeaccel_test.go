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
	P := []Point{{-5, -5, -5},  {5, -5, -5},  {5, -5, 5},  {-5, -5, 5}, }	
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
