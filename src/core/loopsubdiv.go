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
	"container/heap"
	"math"
	"sort"
)

type (
	LoopSubdiv struct {
		ShapeData
		nLevels  int
		vertices []*SDVertex
		faces    []*SDFace
	}

	SDVertex struct {
		P                 Point
		startFace         *SDFace
		child             *SDVertex
		regular, boundary bool
		id                int
	}

	SDFace struct {
		v        [3]*SDVertex
		f        [3]*SDFace
		children [4]*SDFace
	}

	SDEdge struct {
		v         [2]*SDVertex
		f         [2]*SDFace
		f0edgeNum int
	}

	SDEdgeHeap []SDEdge
)

func CreateLoopSubdivShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *LoopSubdiv {
	nlevels := params.FindIntParam("nlevels", 3)
	vi := params.FindIntArrayParam("indices")
	P := params.FindPointArrayParam("P")
	if vi == nil || P == nil {
		return nil
	}

	// don't actually use this for now...
	//scheme := params.FindStringParam("scheme", "loop")

	return NewLoopSubdiv(o2w, w2o, reverseOrientation, len(vi)/3, len(P),
		vi, P, nlevels)
}

var nextVertexId int = 0

func NewLoopSubdiv(o2w, w2o *Transform, ro bool, nfaces, nvertices int, vertexIndices []int, P []Point, nLevels int) *LoopSubdiv {
	subdiv := new(LoopSubdiv)
	subdiv.objectToWorld = o2w
	subdiv.worldToObject = w2o
	subdiv.reverseOrientation = ro
	subdiv.transformSwapsHandedness = SwapsHandednessTransform(subdiv.objectToWorld)
	subdiv.shapeId = GenerateShapeId()

	subdiv.nLevels = nLevels
	// Allocate _LoopSubdiv_ vertices and faces
	subdiv.vertices = make([]*SDVertex, nvertices, nvertices)
	for i := 0; i < nvertices; i++ {
		subdiv.vertices[i] = &SDVertex{P[i], nil, nil, false, false, nextVertexId}
		nextVertexId++
	}
	subdiv.faces = make([]*SDFace, nfaces, nfaces)
	for i := 0; i < nfaces; i++ {
		subdiv.faces[i] = new(SDFace)
	}

	// Set face to vertex pointers
	vp := 0
	for i := 0; i < nfaces; i++ {
		f := subdiv.faces[i]
		for j := 0; j < 3; j++ {
			v := subdiv.vertices[vertexIndices[vp+j]]
			f.v[j] = v
			v.startFace = f
		}
		vp += 3
	}

	// Set neighbor pointers in _faces_
	edges := new(SDEdgeHeap)
	heap.Init(edges)
	for i := 0; i < nfaces; i++ {
		f := subdiv.faces[i]
		for edgeNum := 0; edgeNum < 3; edgeNum++ {
			// Update neighbor pointer for _edgeNum_
			v0, v1 := edgeNum, (edgeNum+1)%3
			e := *NewEdge(f.v[v0], f.v[v1])
			found, ei := edges.Find(e)
			if !found {
				// Handle new edge
				e.f[0] = f
				e.f0edgeNum = edgeNum
				heap.Push(edges, e)
			} else {
				// Handle previously seen edge
				e = (*edges)[ei]
				e.f[0].f[e.f0edgeNum] = f
				f.f[edgeNum] = e.f[0]
				heap.Remove(edges, ei)
			}
		}
	}

	// Finish vertex initialization
	for i := 0; i < nvertices; i++ {
		v := subdiv.vertices[i]
		f := v.startFace
		for {
			f = f.nextFace(v)
			if f == nil || f == v.startFace {
				break
			}
		}
		v.boundary = (f == nil)
		if !v.boundary && v.valence() == 6 {
			v.regular = true
		} else if v.boundary && v.valence() == 4 {
			v.regular = true
		} else {
			v.regular = false
		}
	}

	return subdiv
}

func (subdiv *LoopSubdiv) ObjectBound() *BBox {
	b := CreateEmptyBBox()
	for _, vtx := range subdiv.vertices {
		b = UnionBBoxPoint(b, &vtx.P)
	}
	return b
}

func (subdiv *LoopSubdiv) WorldBound() *BBox {
	b := CreateEmptyBBox()
	for _, vtx := range subdiv.vertices {
		b = UnionBBoxPoint(b, PointTransform(subdiv.objectToWorld, &vtx.P))
	}
	return b
}

func (subdiv *LoopSubdiv) CanIntersect() bool {
	return false
}

func (subdiv *LoopSubdiv) Refine(refined []Shape) []Shape {
	f := subdiv.faces
	v := subdiv.vertices
	for i := 0; i < subdiv.nLevels; i++ {
		// Update _f_ and _v_ for next level of subdivision
		newFaces := make([]*SDFace, 0, len(f))
		newVertices := make([]*SDVertex, 0, len(v))

		// Allocate next level of children in mesh tree
		for j := 0; j < len(v); j++ {
			v[j].child = &SDVertex{Point{0, 0, 0}, nil, nil, v[j].regular, v[j].boundary, nextVertexId}
			nextVertexId++
			newVertices = append(newVertices, v[j].child)
		}
		for j := 0; j < len(f); j++ {
			for k := 0; k < 4; k++ {
				f[j].children[k] = new(SDFace)
				newFaces = append(newFaces, f[j].children[k])
			}
		}
		// Update vertex positions and create new edge vertices

		// Update vertex positions for even vertices
		for j := 0; j < len(v); j++ {
			if !v[j].boundary {
				// Apply one-ring rule for even vertex
				if v[j].regular {
					v[j].child.P = *subdiv.weightOneRing(v[j], 1.0/16.0)
				} else {
					v[j].child.P = *subdiv.weightOneRing(v[j], beta(v[j].valence()))
				}
			} else {
				// Apply boundary rule for even vertex
				v[j].child.P = *subdiv.weightBoundary(v[j], 1.0/8.0)
			}
		}

		// Compute new odd edge vertices
		edgeVerts := make(map[SDEdge]*SDVertex)
		for j := 0; j < len(f); j++ {
			face := f[j]
			for k := 0; k < 3; k++ {
				// Compute odd vertex on _k_th edge
				edge := *NewEdge(face.v[k], face.v[(k+1)%3])
				vert := edgeVerts[edge]
				if vert == nil {
					// Create and initialize new odd vertex
					vert = new(SDVertex)
					newVertices = append(newVertices, vert)
					vert.regular = true
					vert.boundary = (face.f[k] == nil)
					vert.startFace = face.children[3]

					// Apply edge rules to compute new vertex position
					if vert.boundary {
						vert.P = *edge.v[0].P.Scale(0.5)
						vert.P = *vert.P.AddPoint(edge.v[1].P.Scale(0.5))
					} else {
						vert.P = *edge.v[0].P.Scale(3.0 / 8.0)
						vert.P = *vert.P.AddPoint(edge.v[1].P.Scale(3.0 / 8.0))
						vert.P = *vert.P.AddPoint(face.otherVert(edge.v[0], edge.v[1]).P.Scale(1.0 / 8.0))
						vert.P = *vert.P.AddPoint(face.f[k].otherVert(edge.v[0], edge.v[1]).P.Scale(1.0 / 8.0))
					}
					edgeVerts[edge] = vert
				}
			}
		}

		// Update new mesh topology

		// Update even vertex face pointers
		for j := 0; j < len(v); j++ {
			vert := v[j]
			vertNum := vert.startFace.vnum(vert)
			vert.child.startFace =
				vert.startFace.children[vertNum]
		}

		// Update face neighbor pointers
		for j := 0; j < len(f); j++ {
			face := f[j]
			for k := 0; k < 3; k++ {
				// Update children _f_ pointers for siblings
				face.children[3].f[k] = face.children[(k+1)%3]
				face.children[k].f[(k+1)%3] = face.children[3]

				// Update children _f_ pointers for neighbor children
				f2 := face.f[k]
				if f2 != nil {
					face.children[k].f[k] = f2.children[f2.vnum(face.v[k])]
				} else {
					face.children[k].f[k] = nil
				}
				f2 = face.f[(k+2)%3]
				if f2 != nil {
					face.children[k].f[(k+2)%3] = f2.children[f2.vnum(face.v[k])]
				} else {
					face.children[k].f[(k+2)%3] = nil
				}
			}
		}

		// Update face vertex pointers
		for j := 0; j < len(f); j++ {
			face := f[j]
			for k := 0; k < 3; k++ {
				// Update child vertex pointer to new even vertex
				face.children[k].v[k] = face.v[k].child

				// Update child vertex pointer to new odd vertex
				vert := edgeVerts[*NewEdge(face.v[k], face.v[(k+1)%3])]
				face.children[k].v[(k+1)%3] = vert
				face.children[(k+1)%3].v[k] = vert
				face.children[3].v[k] = vert
			}
		}

		// Prepare for next level of subdivision
		f = newFaces
		v = newVertices
	}
	// Push vertices to limit surface
	Plimit := make([]Point, len(v), len(v))
	for i := 0; i < len(v); i++ {
		if v[i].boundary {
			Plimit[i] = *subdiv.weightBoundary(v[i], 1.0/5.0)
		} else {
			Plimit[i] = *subdiv.weightOneRing(v[i], gamma(v[i].valence()))
		}
	}
	for i := 0; i < len(v); i++ {
		v[i].P = Plimit[i]
	}

	// Compute vertex tangents on limit surface
	Ns := make([]Normal, 0, len(v))
	Pring := make([]Point, 16, 16)
	for i := 0; i < len(v); i++ {
		vert := v[i]
		S := CreateVector(0, 0, 0)
		T := CreateVector(0, 0, 0)
		valence := vert.valence()
		if valence > len(Pring) {
			Pring = make([]Point, valence, valence)
		}
		vert.oneRing(Pring)
		if !vert.boundary {
			// Compute tangents of interior face
			for k := 0; k < valence; k++ {
				S = S.Add(CreateVectorFromPoint(&Pring[k]).Scale(math.Cos(2.0 * math.Pi * float64(k) / float64(valence))))
				T = T.Add(CreateVectorFromPoint(&Pring[k]).Scale(math.Sin(2.0 * math.Pi * float64(k) / float64(valence))))
			}
		} else {
			// Compute tangents of boundary face
			S = Pring[valence-1].Sub(&Pring[0])
			if valence == 2 {
				T = Pring[0].AddPoint(&Pring[1]).Sub(vert.P.Scale(2))
			} else if valence == 3 {
				T = Pring[1].Sub(&vert.P)
			} else if valence == 4 { // regular
				T = CreateVectorFromPoint(Pring[0].Scale(-1).AddPoint(Pring[1].Scale(2)).AddPoint(Pring[2].Scale(2)).AddPoint(Pring[3].Scale(-1)).AddPoint(vert.P.Scale(-2)))
			} else {
				theta := math.Pi / float64(valence-1)
				T = CreateVectorFromPoint(Pring[0].AddPoint(&Pring[valence-1]).Scale(math.Sin(theta)))
				for k := 1; k < valence-1; k++ {
					wt := (2.0*math.Cos(theta) - 2.0) * math.Sin(float64(k)*theta)
					T = T.Add(CreateVectorFromPoint(Pring[k].Scale(wt)))
				}
				T = T.Negate()
			}
		}
		Ns = append(Ns, *CreateNormalFromVector(CrossVector(S, T)))
	}

	// Create _TriangleMesh_ from subdivision mesh
	ntris := len(f)
	verts := make([]int, 3*ntris, 3*ntris)
	vp := 0
	totVerts := len(v)
	usedVerts := make(map[*SDVertex]int)
	for i := 0; i < totVerts; i++ {
		usedVerts[v[i]] = i
	}
	for i := 0; i < ntris; i++ {
		for j := 0; j < 3; j++ {
			verts[vp] = usedVerts[f[i].v[j]]
			vp++
		}
	}
	refined = append(refined, NewTriangleMesh(subdiv.objectToWorld,
		subdiv.worldToObject, subdiv.reverseOrientation, verts, Plimit, Ns, nil, nil, nil))

	return refined
}

func (subdiv *LoopSubdiv) Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
	return false, 0.0, 0.0, nil
}
func (subdiv *LoopSubdiv) IntersectP(ray *Ray) bool {
	return false
}
func (subdiv *LoopSubdiv) GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
	return nil
}
func (subdiv *LoopSubdiv) Area() float64 {
	return 0.0
}
func (subdiv *LoopSubdiv) Sample(u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (subdiv *LoopSubdiv) Pdf(pshape *Point) float64 {
	return 0.0
}
func (subdiv *LoopSubdiv) SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
	return nil, nil
}
func (subdiv *LoopSubdiv) Pdf2(p *Point, wi *Vector) float64 {
	return 0.0
}
func (subdiv *LoopSubdiv) ObjectToWorld() *Transform {
	return subdiv.objectToWorld
}
func (subdiv *LoopSubdiv) WorldToObject() *Transform {
	return subdiv.worldToObject
}
func (subdiv *LoopSubdiv) ReverseOrientation() bool {
	return subdiv.reverseOrientation
}
func (subdiv *LoopSubdiv) TransformSwapsHandedness() bool {
	return subdiv.transformSwapsHandedness
}
func (subdiv *LoopSubdiv) ShapeId() uint32 {
	return subdiv.shapeId
}

func (subdiv *LoopSubdiv) weightOneRing(vert *SDVertex, beta float64) *Point {
	// Put _vert_ one-ring in _Pring_
	valence := vert.valence()
	Pring := make([]Point, valence, valence)
	vert.oneRing(Pring)
	P := vert.P.Scale(1.0 - float64(valence)*beta)
	for i := 0; i < valence; i++ {
		P = P.AddPoint(Pring[i].Scale(beta))
	}
	return P
}

func (subdiv *LoopSubdiv) weightBoundary(vert *SDVertex, beta float64) *Point {
	// Put _vert_ one-ring in _Pring_
	valence := vert.valence()
	Pring := make([]Point, valence, valence)
	vert.oneRing(Pring)
	P := vert.P.Scale(1 - 2*beta)
	P = P.AddPoint(Pring[0].Scale(beta))
	P = P.AddPoint(Pring[valence-1].Scale(beta))
	return P
}

func beta(valence int) float64 {
	if valence == 0 {
		return 3.0 / 16.0
	} else {
		return 3.0 / 8.0 * float64(valence)
	}
}

func gamma(valence int) float64 {
	return 1.0 / (float64(valence) + 3.0/(8.0*beta(valence)))
}

func (vtx *SDVertex) valence() int {
	f := vtx.startFace
	if !vtx.boundary {
		// Compute valence of interior vertex
		nf := 1
		f = f.nextFace(vtx)
		for f != vtx.startFace {
			nf++
			f = f.nextFace(vtx)
		}
		return nf
	} else {
		// Compute valence of boundary vertex
		nf := 1
		f = f.nextFace(vtx)
		for f != nil {
			nf++
			f = f.nextFace(vtx)
		}
		f = vtx.startFace
		f = f.prevFace(vtx)
		for f != nil {
			nf++
			f = f.prevFace(vtx)
		}
		return nf + 1
	}
}

func (vtx *SDVertex) oneRing(P []Point) {
	pi := 0
	if !vtx.boundary {
		// Get one-ring vertices for interior vertex
		face := vtx.startFace

		for {
			P[pi] = face.nextVert(vtx).P
			pi++
			face = face.nextFace(vtx)
			if face == vtx.startFace {
				break
			}
		}
	} else {
		// Get one-ring vertices for boundary vertex
		face := vtx.startFace
		f2 := face.nextFace(vtx)
		for f2 != nil {
			face = f2
			f2 = face.nextFace(vtx)
		}
		P[pi] = face.nextVert(vtx).P
		pi++
		for {
			P[pi] = face.prevVert(vtx).P
			pi++
			face = face.prevFace(vtx)
			if face == nil {
				break
			}
		}
	}
}

func (face *SDFace) vnum(vert *SDVertex) int {
	for i := 0; i < 3; i++ {
		if face.v[i] == vert {
			return i
		}
	}
	Severe("Basic logic error in SDFace::vnum()")
	return -1
}
func (face *SDFace) nextFace(vert *SDVertex) *SDFace {
	return face.f[face.vnum(vert)]
}
func (face *SDFace) prevFace(vert *SDVertex) *SDFace {
	return face.f[(face.vnum(vert)+2)%3]
}
func (face *SDFace) nextVert(vert *SDVertex) *SDVertex {
	return face.v[(face.vnum(vert)+1)%3]
}
func (face *SDFace) prevVert(vert *SDVertex) *SDVertex {
	return face.v[(face.vnum(vert)+2)%3]
}
func (face *SDFace) otherVert(v0, v1 *SDVertex) *SDVertex {
	for i := 0; i < 3; i++ {
		if face.v[i] != v0 && face.v[i] != v1 {
			return face.v[i]
		}
	}
	Severe("Basic logic error in SDVertex::otherVert()")
	return nil
}

func NewEdge(v0, v1 *SDVertex) *SDEdge {
	edge := new(SDEdge)
	if v0.id < v1.id {
		edge.v[0] = v0
		edge.v[1] = v1
	} else {
		edge.v[0] = v1
		edge.v[1] = v0
	}
	edge.f[0] = nil
	edge.f[1] = nil
	edge.f0edgeNum = -1

	return edge
}

func (h SDEdgeHeap) Len() int { return len(h) }
func (h SDEdgeHeap) Less(i, j int) bool {
	if h[i].v[0] == h[j].v[0] {
		return h[i].v[1].id < h[j].v[1].id
	}
	return h[i].v[0].id < h[j].v[0].id
}
func (h SDEdgeHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] }

func (h *SDEdgeHeap) Push(x interface{}) {
	*h = append(*h, x.(SDEdge))
}

func (h *SDEdgeHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func (h SDEdgeHeap) Find(edge SDEdge) (found bool, index int) {
	index = sort.Search(len(h), func(i int) bool { return (edge.v[0] == h[i].v[0] && edge.v[1] == h[i].v[1]) })
	return index < len(h), index
}
