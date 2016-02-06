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

import ()

type Heightfield struct {
	ShapeData
	z []float64
	nx, ny int
}

func NewHeightfield(o2w, w2o *Transform, ro bool, nx, ny int, zs []float64) *Heightfield {
	field := new(Heightfield)
	field.objectToWorld = o2w
	field.worldToObject = w2o
	field.reverseOrientation = ro
	field.transformSwapsHandedness = SwapsHandednessTransform(field.objectToWorld)
	field.shapeId = GenerateShapeId()
	field.nx = nx
	field.ny = ny
	field.z = make([]float64, len(zs), len(zs))
	copy(field.z, zs)
	
	return field
}

func (field *Heightfield) ObjectBound() *BBox {
    minz, maxz := field.z[0], field.z[0]
    for _, z := range field.z {
        if z < minz { minz = z }
        if z > maxz { maxz = z }
    }
    return CreateBBoxFromPoints(CreatePoint(0,0,minz), CreatePoint(1,1,maxz))
}

func (field *Heightfield) WorldBound() *BBox {
	return BBoxTransform(field.objectToWorld, field.ObjectBound())
}
func (field *Heightfield) CanIntersect() bool {
	return false
}
func (field *Heightfield) Refine(refined []Shape) []Shape {
    ntris := 2*(field.nx-1)*(field.ny-1)
    nverts := field.nx*field.ny
   
    verts := make([]int, 3*ntris, 3*ntris)
    P := make([]Point, nverts, nverts)
    uvs := make([]float64, 2*nverts, 2*nverts)
    // Compute heightfield vertex positions
    pos := 0
    for y := 0; y < field.ny; y++ {
        for x := 0; x < field.nx; x++ {
            P[pos].X = float64(x) / float64(field.nx-1)
            P[pos].Y = float64(y) / float64(field.ny-1)
            P[pos].Z = field.z[pos]
            uvs[2*pos] = P[pos].X
            uvs[2*pos+1] = P[pos].Y
            pos++
        }
    }

    // Fill in heightfield vertex offset array
    VERT := func(x,y int) int { return x+y*field.nx }
    vi := 0
    for y := 0; y < field.ny-1; y++ {
        for x := 0; x < field.nx-1; x++ {
            verts[vi] = VERT(x, y); vi++
            verts[vi] = VERT(x+1, y); vi++
            verts[vi] = VERT(x+1, y+1); vi++
    
            verts[vi] = VERT(x, y); vi++
            verts[vi] = VERT(x+1, y+1); vi++
            verts[vi] = VERT(x, y+1); vi++
        }
    }
	refined = append(refined, NewTriangleMesh(field.objectToWorld,
		field.worldToObject, field.reverseOrientation, verts, P, nil, nil, uvs, nil))

	return refined
}
func (field *Heightfield) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (field *Heightfield) IntersectP(ray *Ray) bool {
	return false
}
func (field *Heightfield) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (field *Heightfield) Area() float64 {
	return 0.0
}
func (field *Heightfield) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (field *Heightfield) Pdf(pshape *Point) float64 {
	return 0.0
}
func (field *Heightfield) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (field *Heightfield) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (field *Heightfield) ObjectToWorld() *Transform {
	return field.objectToWorld
}
func (field *Heightfield) WorldToObject() *Transform {
	return field.worldToObject
}
func (field *Heightfield) ReverseOrientation() bool {
	return field.reverseOrientation
}
func (field *Heightfield) TransformSwapsHandedness() bool {
	return field.transformSwapsHandedness
}
func (field *Heightfield) ShapeId() uint32 {
	return field.shapeId
}

func CreateHeightfieldShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Heightfield {
    nu := params.FindIntParam("nu", -1)
    nv := params.FindIntParam("nv", -1)
    Pz := params.FindFloatArrayParam("Pz")
    Assert(nu != -1 && nv != -1 && Pz != nil)
    Assert(len(Pz) == nu*nv)
    return NewHeightfield(o2w, w2o, reverseOrientation, nu, nv, Pz)
}
