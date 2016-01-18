package core

import ()

type TriangleMesh struct {
	ShapeData
	ntris, nverts int
	vertexIndex []int
	p []Point
	n []Normal
	s []Vector
	uvs []float64
	alphaTexture *TextureFloat
}

type Triangle struct {
	ShapeData
	mesh *TriangleMesh
	v []int
}

func CreateTriangleMeshShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet, textures map[string]TextureFloat) *TriangleMesh {
	return nil
}

func NewTriangleMesh(o2w, w2o *Transform, reverseOrientation bool, nt, nv int, vi []int, P []Point, N []Normal, S []Vector, uv []float64, atex *TextureFloat) *TriangleMesh {
	t := new(TriangleMesh)
	t.objectToWorld = o2w
	t.worldToObject = w2o
	t.reverseOrientation = reverseOrientation
	t.transformSwapsHandedness = SwapsHandednessTransform(t.objectToWorld)
	t.shapeId = GenerateShapeId()
	t.alphaTexture = atex
	t.ntris = nt
	t.nverts = nv
	t.vertexIndex = make([]int, 3*t.ntris, 3*t.ntris)
	
   	if uv != nil && len(uv) == 2*t.nverts {
        t.uvs = make([]float64, 0, 2*t.nverts)
        copy(t.uvs, uv)
    } else {
    	t.uvs = nil
   	}
    t.p = make([]Point, 0, t.nverts)
    if N != nil {
        t.n = make([]Normal, 0, t.nverts)
        copy(t.n, N)
    } else { 
    	t.n = nil
   	}
    if S != nil {
        t.s = make([]Vector, 0, t.nverts)
        copy(t.s, S)
    } else {
    	t.s = nil
	}
    // Transform mesh vertices to world space
    for i := 0; i < t.nverts; i++ {
        t.p[i] = *PointTransform(t.objectToWorld, &P[i])
	}
	return t	
}

func (t *TriangleMesh) ObjectBound() *BBox {
    objectBounds := CreateEmptyBBox()
    for i := 0; i < t.nverts; i++ {
        objectBounds = UnionBBoxPoint(objectBounds, PointTransform(t.worldToObject, &t.p[i]))
    }    
    return objectBounds
}

func (t *TriangleMesh) WorldBound() *BBox {
    worldBounds := CreateEmptyBBox()
    for i := 0; i < t.nverts; i++ {
        worldBounds = UnionBBoxPoint(worldBounds, PointTransform(t.worldToObject, &t.p[i]))
    }    
    return worldBounds
}

func (t *TriangleMesh) CanIntersect() bool {
	return true
}
func (t *TriangleMesh) Refine(refined []Shape) []Shape {
    for i := 0; i < t.ntris; i++ {
        refined = append(refined, NewTriangle(t.objectToWorld, t.worldToObject, t.reverseOrientation, t, i))
	}
	return refined
}

func (t *TriangleMesh) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (t *TriangleMesh) IntersectP(ray *Ray) bool {
	return false
}
func (t *TriangleMesh) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (t *TriangleMesh) Area() float64 {
	return 0.0
}
func (t *TriangleMesh) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *TriangleMesh) Pdf(pshape *Point) float64 {
	return 0.0
}
func (t *TriangleMesh) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *TriangleMesh) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (t *TriangleMesh) ObjectToWorld() *Transform {
	return t.objectToWorld
}
func (t *TriangleMesh) WorldToObject() *Transform {
	return t.worldToObject
}
func (t *TriangleMesh) ReverseOrientation() bool {
	return t.reverseOrientation
}
func (t *TriangleMesh) TransformSwapsHandedness() bool {
	return t.transformSwapsHandedness
}
func (t *TriangleMesh) ShapeId() uint32 {
	return t.shapeId
}

func NewTriangle(o2w, w2o *Transform, ro bool, m *TriangleMesh, n int) *Triangle {
	return nil
}
func (t *Triangle) ObjectBound() *BBox {
	return nil
}

func (t *Triangle) WorldBound() *BBox {
	return nil
}
func (t *Triangle) CanIntersect() bool {
	return false
}
func (t *Triangle) Refine(refined []Shape) []Shape {
	return refined
}
func (t *Triangle) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (t *Triangle) IntersectP(ray *Ray) bool {
	return false
}
func (t *Triangle) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (t *Triangle) Area() float64 {
	return 0.0
}
func (t *Triangle) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *Triangle) Pdf(pshape *Point) float64 {
	return 0.0
}
func (t *Triangle) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (t *Triangle) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (t *Triangle) ObjectToWorld() *Transform {
	return t.objectToWorld
}
func (t *Triangle) WorldToObject() *Transform {
	return t.worldToObject
}
func (t *Triangle) ReverseOrientation() bool {
	return t.reverseOrientation
}
func (t *Triangle) TransformSwapsHandedness() bool {
	return t.transformSwapsHandedness
}
func (t *Triangle) ShapeId() uint32 {
	return t.mesh.shapeId
}

func (t *Triangle) GetUVs(uv [3][2]float64) {
	
}