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

// Octree Declarations
type (
	OctNode struct {
		children [8]*OctNode
		data     []Object
	}

	Octree struct {
		maxDepth int
		bound    BBox
		root     OctNode
	}

	OctreeLookupProc func(data Object) bool
)

func NewOctree(bbox *BBox, maxDepth int) *Octree {
	octree := new(Octree)
	octree.maxDepth = maxDepth
	octree.bound = *bbox

	return octree
}

func (octree *Octree) Add(dataItem Object, dataBound *BBox) {
	octree.addPrivate(&octree.root, &octree.bound, dataItem, dataBound, DistanceSquaredPoint(&dataBound.pMin, &dataBound.pMax), 0)
}

func (octree *Octree) Lookup(p *Point, process OctreeLookupProc) bool {
	if !octree.bound.Inside(p) {
		return false
	}
	return octree.lookupPrivate(&octree.root, &octree.bound, p, process)
}

func octreeChildBound(child int, nodeBound *BBox, pMid *Point) BBox {
	var childBound BBox
	if child&4 != 0 {
		childBound.pMin.x = pMid.x
		childBound.pMax.x = nodeBound.pMax.x
	} else {
		childBound.pMin.x = nodeBound.pMin.x
		childBound.pMax.x = pMid.x
	}
	if child&2 != 0 {
		childBound.pMin.y = pMid.y
		childBound.pMax.y = nodeBound.pMax.y
	} else {
		childBound.pMin.y = nodeBound.pMin.y
		childBound.pMax.y = pMid.y
	}
	if child&1 != 0 {
		childBound.pMin.z = pMid.z
		childBound.pMax.z = nodeBound.pMax.z
	} else {
		childBound.pMin.z = nodeBound.pMin.z
		childBound.pMax.z = pMid.z
	}
	return childBound
}

// Octree Method Definitions
func (octree *Octree) addPrivate(node *OctNode, nodeBound *BBox, dataItem Object, dataBound *BBox, diag2 float64, depth int) {
	// Possibly add data item to current octree node
	if depth == octree.maxDepth ||
		DistanceSquaredPoint(&nodeBound.pMin, &nodeBound.pMax) < diag2 {
		node.data = append(node.data, dataItem)
		return
	}

	// Otherwise add data item to octree children
	pMid := nodeBound.pMin.AddPoint(&nodeBound.pMax).Scale(0.5)

	// Determine which children the item overlaps
	x := [2]bool{dataBound.pMin.x <= pMid.x, dataBound.pMax.x > pMid.x}
	y := [2]bool{dataBound.pMin.y <= pMid.y, dataBound.pMax.y > pMid.y}
	z := [2]bool{dataBound.pMin.z <= pMid.z, dataBound.pMax.z > pMid.z}
	over := [8]bool{bool(x[0] && y[0] && z[0]), bool(x[0] && y[0] && z[1]),
		bool(x[0] && y[1] && z[0]), bool(x[0] && y[1] && z[1]),
		bool(x[1] && y[0] && z[0]), bool(x[1] && y[0] && z[1]),
		bool(x[1] && y[1] && z[0]), bool(x[1] && y[1] && z[1])}
	for child := 0; child < 8; child++ {
		if !over[child] {
			continue
		}
		// Allocate octree node if needed and continue recursive traversal
		if node.children[child] == nil {
			node.children[child] = new(OctNode)
		}
		childBound := octreeChildBound(child, nodeBound, pMid)
		octree.addPrivate(node.children[child], &childBound,
			dataItem, dataBound, diag2, depth+1)
	}
}

func (octree *Octree) lookupPrivate(node *OctNode, nodeBound *BBox, p *Point, process OctreeLookupProc) bool {
	for i := 0; i < len(node.data); i++ {
		if !process(node.data[i]) {
			return false
		}
	}
	// Determine which octree child node _p_ is inside
	pMid := nodeBound.pMin.AddPoint(&nodeBound.pMax).Scale(0.5)
	child := 0
	if p.x > pMid.x {
		child += 4
	}
	if p.y > pMid.y {
		child += 2
	}
	if p.z > pMid.z {
		child += 1
	}

	if node.children[child] == nil {
		return true
	}
	childBound := octreeChildBound(child, nodeBound, pMid)
	return octree.lookupPrivate(node.children[child], &childBound, p, process)
}
