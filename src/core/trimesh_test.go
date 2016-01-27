package core

import (
	"fmt"
	"testing"
)

func TestTriangleMeshShape(t *testing.T) {

	ident := NewIdentityMatrix4x4()
	fmt.Printf("Ident: %v\n", ident)

	o2w, _ := NewTransform(ident)
	w2o := InverseTransform(o2w)
	
	vi := []int{0, 1, 2, 2, 3, 0}
	P := []Point{{-5, 0, -5},  {5, 0, -5},  {5, 0, 5},  {-5, 0, 5}, }
	
	mesh := NewTriangleMesh(o2w, w2o, false, vi, P, nil, nil, nil, nil)

	fmt.Printf("Mesh: %v\n", mesh)

	ray := CreateRay(&Point{0.0, 10.0, 0.0}, &Vector{0.0, -1.0, 0.0}, 0.0, 100.0, 0.0, 0)
	fmt.Printf("Ray: %v\n", ray)

	tris := make([]Shape, 0, 2)
	tris = mesh.Refine(tris)
	if tris == nil || len(tris) != 2 {
		t.Fail()
	}
	
	var hit bool = false
	for _, tri := range tris {
		if tri.IntersectP(ray) { 
			hit = true
		} 
		hit2,_,_,_ := tri.Intersect(ray)
		if hit != hit2 { t.Fail() }
	}
	// should hit one of the two triangle
	if !hit { t.Fail() }
	
	prim := NewGeometricPrimitive(mesh, nil, nil)
	prims := make([]Primitive, 0, 2)
	prims = prim.Refine(prims)

	hit = false
	for _, pm := range prims {
		if pm.IntersectP(ray) { 
			hit = true
		} 
		hit2,_ := pm.Intersect(CreateRayDifferentialFromRay(ray))
		if hit != hit2 { t.Fail() }
	}
	// should hit one of the two triangle
	if !hit { t.Fail() }
	
}
