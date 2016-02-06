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
	"math"
	"sort"
)

type SplitAxis int

const (
	SPLIT_X    SplitAxis = 0
	SPLIT_Y              = 1
	SPLIT_Z              = 2
	SPLIT_LEAF           = 3
)

var SplitAxisNames map[SplitAxis]string = map[SplitAxis]string{
	SPLIT_X:    "X",
	SPLIT_Y:    "Y",
	SPLIT_Z:    "Z",
	SPLIT_LEAF: "leaf",
}

const (
	NO_RIGHT_CHILD = math.MaxInt64
)

type (
	NodeData interface {
		Location() *Point
	}

	KdNode struct {
		// KdNode Data
		splitPos     float64
		splitAxis    SplitAxis
		hasLeftChild bool
		rightChild   int
	}

	KdTree struct {
		nodes                []KdNode
		nodeData             []NodeData
		nNodes, nextFreeNode int
	}

	KdTreeLookupProc func(p *Point, nodeData NodeData, dist2 float64, maxDistSquared *float64)
)

func NewKdTree(data []NodeData) *KdTree {
	kdtree := new(KdTree)

	kdtree.nNodes = len(data)
	kdtree.nextFreeNode = 1
	kdtree.nodes = make([]KdNode, kdtree.nNodes, kdtree.nNodes)
	kdtree.nodeData = make([]NodeData, kdtree.nNodes, kdtree.nNodes)
	
	buildNodes := make([]NodeData, kdtree.nNodes, kdtree.nNodes)
	copy(buildNodes, data)

	kdtree.recursiveBuild(0, 0, kdtree.nNodes, buildNodes)

	return kdtree
}

func (kdtree *KdTree) Lookup(p *Point, proc KdTreeLookupProc, maxDistSquared *float64) {
	kdtree.privateLookup(0, p, proc, maxDistSquared)
}

func (kdtree *KdTree) recursiveBuild(nodeNum, start, end int, buildNodes []NodeData) {
	// Create leaf node of kd-tree if we've reached the bottom
	if start+1 == end {
		kdtree.nodes[nodeNum].initLeaf()
		kdtree.nodeData[nodeNum] = buildNodes[start]
		return
	}

	// Choose split direction and partition data

	// Compute bounds of data from _start_ to _end_
	bound := CreateEmptyBBox()
	for i := start; i < end; i++ {
		bound = UnionBBoxPoint(bound, buildNodes[i].Location())
	}
	splitAxis := bound.MaximumExtent()
	splitPos := (start + end) / 2

	kdtree_nth_element(buildNodes, start, splitPos, end, splitAxis)

	// Allocate kd-tree node and continue recursively
	kdtree.nodes[nodeNum].init(buildNodes[splitPos].Location().At(splitAxis), SplitAxis(splitAxis))
	kdtree.nodeData[nodeNum] = buildNodes[splitPos]
	if start < splitPos {
		kdtree.nodes[nodeNum].hasLeftChild = true
		childNum := kdtree.nextFreeNode
		kdtree.nextFreeNode++
		kdtree.recursiveBuild(childNum, start, splitPos, buildNodes)
	}
	if splitPos+1 < end {
		kdtree.nodes[nodeNum].rightChild = kdtree.nextFreeNode
		kdtree.nextFreeNode++
		kdtree.recursiveBuild(kdtree.nodes[nodeNum].rightChild, splitPos+1, end, buildNodes)
	}
}

func (kdtree *KdTree) privateLookup(nodeNum int, p *Point, process KdTreeLookupProc, maxDistSquared *float64) {
	node := &kdtree.nodes[nodeNum]
	// Process kd-tree node's children
	axis := node.splitAxis
	if axis != SPLIT_LEAF {
		dist2 := (p.At(int(axis)) - node.splitPos) * (p.At(int(axis)) - node.splitPos)
		if p.At(int(axis)) <= node.splitPos {
			if node.hasLeftChild {
				kdtree.privateLookup(nodeNum+1, p, process, maxDistSquared)
			}
			if dist2 < *maxDistSquared && node.rightChild < kdtree.nNodes {
				kdtree.privateLookup(node.rightChild, p, process, maxDistSquared)
			}
		} else {
			if node.rightChild < kdtree.nNodes {
				kdtree.privateLookup(node.rightChild, p, process, maxDistSquared)
			}
			if dist2 < *maxDistSquared && node.hasLeftChild {
				kdtree.privateLookup(nodeNum+1, p, process, maxDistSquared)
			}
		}
	}

	// Hand kd-tree node to processing function
	dist2 := DistanceSquaredPoint(kdtree.nodeData[nodeNum].Location(), p)
	if dist2 < *maxDistSquared {
		process(p, kdtree.nodeData[nodeNum], dist2, maxDistSquared)
	}
}

type kdNodeSorter struct {
	axis  int
	nodes []NodeData
}

func (s *kdNodeSorter) Len() int {
	return len(s.nodes)
}
func (s *kdNodeSorter) Swap(i, j int) {
	s.nodes[i], s.nodes[j] = s.nodes[j], s.nodes[i]
}
func (s *kdNodeSorter) Less(i, j int) bool {
	return s.nodes[i].Location().At(s.axis) < s.nodes[j].Location().At(s.axis)
}

func kdtree_nth_element(nodes []NodeData, first, nth, last int, axis int) {
	nodeSorter := &kdNodeSorter{axis, nodes[first:last]}
	sort.Sort(nodeSorter)
}

func (kdnode *KdNode) init(p float64, axis SplitAxis) {
	kdnode.splitPos = p
	kdnode.splitAxis = axis
	kdnode.rightChild = NO_RIGHT_CHILD
	kdnode.hasLeftChild = false
}

func (kdnode *KdNode) initLeaf() {
	kdnode.splitAxis = SPLIT_LEAF
	kdnode.rightChild = NO_RIGHT_CHILD
	kdnode.hasLeftChild = false
}
