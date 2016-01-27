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

	ray := CreateRay(&Point{1.0, 10.0, -2.0}, &Vector{0.0, -1.0, 0.0}, 0.0, 100.0, 0.0, 0)
	fmt.Printf("Ray: %v\n", ray)
	raydiff := CreateRayDifferentialFromRay(ray)

	tris := make([]Shape, 0, 2)
	tris = mesh.Refine(tris)
	if tris == nil || len(tris) != 2 {
		t.Error("Triangle mesh did not refined into 2 triangles.")
	}
	
	var gotahit bool = false
	for _, s := range tris {
		hit := s.IntersectP(ray)
		hit2,_,_,dg := s.Intersect(ray)
		if hit != hit2 { 
			t.Errorf("Inconsistant triangle Intersect and IntersectP tests.") 
			}
		if hit2 {
			fmt.Printf("DiffGeom: %v\n", dg)
		}
		gotahit = gotahit || hit
		
		tri, ok := s.(*Triangle)
		if ok {
			var uvs [3][2]float64
			tri.GetUVs(&uvs)
			fmt.Printf("Tri UVS: %v\n", uvs)
		}
	}
	// should hit one of the two triangles
	if !gotahit { t.Errorf("Ray did not hit the mesh.") }
	
	mesh = NewTriangleMesh(o2w, w2o, false, vi, P, nil, nil, nil, nil)
	prim := NewGeometricPrimitive(mesh, nil, nil)
	prims := make([]Primitive, 0, 2)
	prims = prim.Refine(prims)

	gotahit = false
	for _, pm := range prims {
		hit := pm.IntersectP(ray)
		hit2,isect := pm.Intersect(raydiff)
		if hit != hit2 { t.Errorf("Inconsistant triangle Intersect and IntersectP tests.")  }
		if hit2 {
			fmt.Printf("DiffGeom: %v\n", isect.dg)
		}
		gotahit = gotahit || hit
	}
	// should hit one of the two triangles
	if !gotahit { t.Errorf("Ray did not hit the mesh.") }
	
}
