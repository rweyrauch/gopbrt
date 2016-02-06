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
	"sort"
)

type EdgeType int

const (
	EDGE_START EdgeType = iota
	EDGE_END
)

type (
	KdTreeAccel struct {
		PrimitiveData
		isectCost, traversalCost, maxPrims, maxDepth int
		emptyBonus                                   float64
		primitives                                   []Primitive
		nodes                                        []KdAccelNode
		nAllocedNodes, nextFreeNode                  int
		bounds                                       BBox
		arena                                        *MemoryArena
	}

	KdAccelNode struct {
		axis       SplitAxis
		split      float64 // Interior
		aboveChild int     // Interior

		onePrimitive int   // Leaf
		primitives   []int // Leaf
		nPrims       int   // Leaf
	}

	BoundEdge struct {
		t        float64
		primNum  int
		edgeType EdgeType
	}

	KdToDo struct {
		nodeIdx    int
		tmin, tmax float64
	}
)

func NewKdTreeAccel(prims []Primitive, icost, tcost int, ebonus float64, maxp, md int) *KdTreeAccel {
	accel := new(KdTreeAccel)
	accel.primitiveId = GeneratePrimitiveId()
	accel.primitives = make([]Primitive, 0, 16)
	accel.isectCost = icost
	accel.traversalCost = tcost
	accel.maxPrims = maxp
	accel.maxDepth = md
	accel.emptyBonus = ebonus

	//PBRT_KDTREE_STARTED_CONSTRUCTION(this, p.size());
	for _, p := range prims {
		accel.primitives = p.FullyRefine(accel.primitives)
	}
	// Build kd-tree for accelerator
	accel.nextFreeNode = 0
	accel.nAllocedNodes = 0
	if accel.maxDepth <= 0 {
		accel.maxDepth = int(Round2Int(8 + 1.3*float64(Log2Int(float64(len(accel.primitives))))))
	}
	// Compute bounds for kd-tree construction
	accel.bounds = *CreateEmptyBBox()
	primBounds := make([]BBox, 0, len(accel.primitives))
	for _, p := range accel.primitives {
		b := p.WorldBound()
		accel.bounds = *UnionBBoxes(&accel.bounds, b)
		primBounds = append(primBounds, *b)
	}

	// Allocate working memory for kd-tree construction
	var edges [3][]BoundEdge
	for i := 0; i < 3; i++ {
		edges[i] = make([]BoundEdge, 2*len(accel.primitives), 2*len(accel.primitives))
	}
	prims0 := make([]int, len(accel.primitives), len(accel.primitives))
	prims1 := make([]int, (accel.maxDepth+1)*len(accel.primitives), (accel.maxDepth+1)*len(accel.primitives))

	// Initialize _primNums_ for kd-tree construction
	primNums := make([]int, len(accel.primitives), len(accel.primitives))
	for i := 0; i < len(accel.primitives); i++ {
		primNums[i] = i
	}

	// Start recursive construction of kd-tree
	accel.buildTree(0, accel.bounds, primBounds, primNums, len(accel.primitives),
		accel.maxDepth, &edges, prims0, prims1, 0)

	Info("KdTreeAccel: %s", accel)
	for i := 0; i < accel.nextFreeNode-1; i++ {
		Info("Node[%d]: %s", i, &accel.nodes[i])
	}
	//PBRT_KDTREE_FINISHED_CONSTRUCTION(this);

	return accel
}

func (accel *KdTreeAccel) WorldBound() *BBox { return &accel.bounds }

func (*KdTreeAccel) CanIntersect() bool { return true }

const (
	MAX_TODO = 64
)

func (accel *KdTreeAccel) Intersect(ray *RayDifferential) (hit bool, isect *Intersection) {
	//PBRT_KDTREE_INTERSECTION_TEST(const_cast<KdTreeAccel *>(this), const_cast<Ray *>(&ray));
	// Compute initial parametric range of ray inside kd-tree extent
	var tmin, tmax float64
	var boxHit bool
	if boxHit, tmin, tmax = accel.bounds.IntersectP(CreateRayFromRayDifferential(ray)); !boxHit {
		//PBRT_KDTREE_RAY_MISSED_BOUNDS();
		return false, nil
	}

	// Prepare to traverse kd-tree for ray
	invDir := CreateVector(1.0/ray.Dir.X, 1.0/ray.Dir.Y, 1.0/ray.Dir.Z)

	var todo [MAX_TODO]KdToDo
	todoPos := 0

	// Traverse kd-tree nodes in order for ray
	prevNodeIdx := -1
	curNodeIdx := 0
	for prevNodeIdx != curNodeIdx {
		node := &accel.nodes[curNodeIdx]
		prevNodeIdx = curNodeIdx

		// Bail out if we found a hit closer than the current node
		if ray.Maxt < tmin {
			break
		}
		if node.axis != SPLIT_LEAF {
			//PBRT_KDTREE_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<KdAccelNode *>(node));
			// Process kd-tree interior node

			// Compute parametric distance along ray to split plane
			axis := int(node.axis)
			tplane := (node.split - ray.Origin.At(axis)) * invDir.At(axis)

			// Get node children pointers for ray
			var firstChild, secondChild int
			belowFirst := (ray.Origin.At(axis) < node.split) || (ray.Origin.At(axis) == node.split && ray.Dir.At(axis) <= 0)
			if belowFirst {
				firstChild = curNodeIdx + 1
				secondChild = node.aboveChild
			} else {
				firstChild = node.aboveChild
				secondChild = curNodeIdx + 1
			}

			// Advance to next child node, possibly enqueue other child
			if tplane > tmax || tplane <= 0 {
				curNodeIdx = firstChild
			} else if tplane < tmin {
				curNodeIdx = secondChild
			} else {
				// Enqueue _secondChild_ in todo list
				todo[todoPos].nodeIdx = secondChild
				todo[todoPos].tmin = tplane
				todo[todoPos].tmax = tmax
				todoPos++
				curNodeIdx = firstChild
				tmax = tplane
			}
		} else {
			//PBRT_KDTREE_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<KdAccelNode *>(node), node->nPrimitives());
			// Check for intersections inside leaf node
			nPrimitives := node.nPrims
			if nPrimitives == 1 {
				prim := accel.primitives[node.onePrimitive]
				// Check one primitive inside leaf node
				//PBRT_KDTREE_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
				if rayHit, rayIsect := prim.Intersect(ray); rayHit {
					//PBRT_KDTREE_INTERSECTION_HIT(const_cast<Primitive *>(prim.GetPtr()));
					hit = true
					isect = rayIsect
				}
			} else {
				for _, ip := range node.primitives {
					prim := accel.primitives[ip]
					// Check one primitive inside leaf node
					//PBRT_KDTREE_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
					if rayHit, rayIsect := prim.Intersect(ray); rayHit {
						//PBRT_KDTREE_INTERSECTION_HIT(const_cast<Primitive *>(prim.GetPtr()));
						hit = true
						isect = rayIsect
					}
				}
			}

			// Grab next node to process from todo list
			if todoPos > 0 {
				todoPos--
				curNodeIdx = todo[todoPos].nodeIdx
				tmin = todo[todoPos].tmin
				tmax = todo[todoPos].tmax
			} else {
				break
			}
		}
	}
	//PBRT_KDTREE_INTERSECTION_FINISHED();
	//if hit { Assert(isect != nil) }
	return hit, isect
}

func (accel *KdTreeAccel) IntersectP(ray *Ray) bool {
	//PBRT_KDTREE_INTERSECTIONP_TEST(const_cast<KdTreeAccel *>(this), const_cast<Ray *>(&ray));
	// Compute initial parametric range of ray inside kd-tree extent
	var tmin, tmax float64
	var hit bool
	if hit, tmin, tmax = accel.bounds.IntersectP(ray); !hit {
		//PBRT_KDTREE_RAY_MISSED_BOUNDS();
		return false
	}

	// Prepare to traverse kd-tree for ray
	invDir := CreateVector(1.0/ray.Dir.X, 1.0/ray.Dir.Y, 1.0/ray.Dir.Z)

	var todo [MAX_TODO]KdToDo
	todoPos := 0

	prevNodeIdx := -1
	curNodeIdx := 0
	for prevNodeIdx != curNodeIdx {
		node := &accel.nodes[curNodeIdx]
		prevNodeIdx = curNodeIdx

		if node.axis == SPLIT_LEAF {
			//PBRT_KDTREE_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<KdAccelNode *>(node), node->nPrimitives());
			// Check for shadow ray intersections inside leaf node
			if node.nPrims == 1 {
				prim := accel.primitives[node.onePrimitive]
				//PBRT_KDTREE_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
				if hit = prim.IntersectP(ray); hit {
					//PBRT_KDTREE_INTERSECTIONP_HIT(const_cast<Primitive *>(prim.GetPtr()));
					return true
				}
			} else {
				for _, ip := range node.primitives {
					prim := accel.primitives[ip]
					//PBRT_KDTREE_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
					if hit = prim.IntersectP(ray); hit {
						//PBRT_KDTREE_INTERSECTIONP_HIT(const_cast<Primitive *>(prim.GetPtr()));
						return true
					}
				}
			}

			// Grab next node to process from todo list
			if todoPos > 0 {
				todoPos--
				curNodeIdx = todo[todoPos].nodeIdx
				tmin = todo[todoPos].tmin
				tmax = todo[todoPos].tmax
			} else {
				break
			}
		} else {
			//PBRT_KDTREE_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<KdAccelNode *>(node));
			// Process kd-tree interior node

			// Compute parametric distance along ray to split plane
			axis := int(node.axis)
			tplane := (node.split - ray.Origin.At(axis)) * invDir.At(axis)

			// Get node children pointers for ray
			var firstChild, secondChild int
			belowFirst := (ray.Origin.At(axis) < node.split) || (ray.Origin.At(axis) == node.split && ray.Dir.At(axis) <= 0)
			if belowFirst {
				firstChild = curNodeIdx + 1
				secondChild = node.aboveChild
			} else {
				firstChild = node.aboveChild
				secondChild = curNodeIdx + 1
			}

			// Advance to next child node, possibly enqueue other child
			if tplane > tmax || tplane <= 0 {
				curNodeIdx = firstChild
			} else if tplane < tmin {
				curNodeIdx = secondChild
			} else {
				// Enqueue _secondChild_ in todo list
				todo[todoPos].nodeIdx = secondChild
				todo[todoPos].tmin = tplane
				todo[todoPos].tmax = tmax
				todoPos++
				curNodeIdx = firstChild
				tmax = tplane
			}
		}
	}
	//PBRT_KDTREE_INTERSECTIONP_MISSED();
	return false
}

func (*KdTreeAccel) Refine(refined []Primitive) []Primitive {
	Severe("Unimplemented KdTreeAccel::Refine() method called!")
	return nil
}

func (accel *KdTreeAccel) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(accel, refined)
}

func (*KdTreeAccel) GetAreaLight() AreaLight {
	Severe("KdTreeAccel::GetAreaLight() method called; should have gone to GeometricPrimitive")
	return nil
}

func (*KdTreeAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	Severe("KdTreeAccel::GetBSDF() method called; should have gone to GeometricPrimitive")
	return nil
}

func (p *KdTreeAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	Severe("KdTreeAccel::GetBSSRDF() method called; should have gone to GeometricPrimitive")
	return nil
}

func (accel *KdTreeAccel) PrimitiveId() uint32 { return accel.primitiveId }

func CreateKdTreeAccelerator(prims []Primitive, ps *ParamSet) *KdTreeAccel {
	isectCost := ps.FindIntParam("intersectcost", 80)
	travCost := ps.FindIntParam("traversalcost", 1)
	emptyBonus := ps.FindFloatParam("emptybonus", 0.5)
	maxPrims := ps.FindIntParam("maxprims", 1)
	maxDepth := ps.FindIntParam("maxdepth", -1)
	return NewKdTreeAccel(prims, isectCost, travCost, emptyBonus, maxPrims, maxDepth)
}

func (accel *KdTreeAccel) String() string {
	return fmt.Sprintf("kdtree[icost:%d tcost:%d maxp:%d maxd:%d empty:%f numprims:%d allocnodes:%d nextnode:%d bounds:%s]",
		accel.isectCost, accel.traversalCost, accel.maxPrims, accel.maxDepth, accel.emptyBonus, len(accel.primitives), accel.nAllocedNodes, accel.nextFreeNode, &accel.bounds)
}

func (accel *KdTreeAccel) buildTree(nodeNum int, nodeBounds BBox, allPrimBounds []BBox, primNums []int, nPrimitives, depth int, edges *[3][]BoundEdge, prims0, prims1 []int, badRefines int) {
	Assert(nodeNum == accel.nextFreeNode)
	// Get next free node from _nodes_ array
	if accel.nextFreeNode == accel.nAllocedNodes {
		nAlloc := Maxi(2*accel.nAllocedNodes, 512)
		newNodes := make([]KdAccelNode, nAlloc, nAlloc)
		if accel.nAllocedNodes > 0 {
			copy(newNodes, accel.nodes)
		}
		accel.nodes = newNodes
		accel.nAllocedNodes = nAlloc
	}
	accel.nextFreeNode++

	// Initialize leaf node if termination criteria met
	if nPrimitives <= accel.maxPrims || depth == 0 {
		//PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth-depth);
		accel.nodes[nodeNum].initLeaf(primNums, nPrimitives, accel.arena)
		return
	}

	// Initialize interior node and continue recursion

	// Choose split axis position for interior node
	bestAxis, bestOffset := -1, -1
	bestCost := INFINITY
	oldCost := float64(accel.isectCost) * float64(nPrimitives)
	totalSA := nodeBounds.SurfaceArea()
	invTotalSA := 1.0 / totalSA
	d := nodeBounds.PMax.Sub(&nodeBounds.PMin)

	// Choose which axis to split along
	axis := nodeBounds.MaximumExtent()
	retries := 0

retrySplit:

	// Initialize edges for _axis_
	for i := 0; i < nPrimitives; i++ {
		pn := primNums[i]
		bbox := allPrimBounds[pn]
		edges[axis][2*i] = BoundEdge{bbox.PMin.At(axis), pn, EDGE_START}
		edges[axis][2*i+1] = BoundEdge{bbox.PMax.At(axis), pn, EDGE_END}
	}
	edgeSorter := &boundEdgeSorter{axis, edges[axis][:2*nPrimitives]}
	sort.Sort(edgeSorter)

	// Compute cost of all splits for _axis_ to find best
	nBelow, nAbove := 0, nPrimitives
	for i := 0; i < 2*nPrimitives; i++ {
		if edges[axis][i].edgeType == EDGE_END {
			nAbove--
		}
		edget := edges[axis][i].t
		if edget > nodeBounds.PMin.At(axis) && edget < nodeBounds.PMax.At(axis) {
			// Compute cost for split at _i_th edge
			otherAxis0 := (axis + 1) % 3
			otherAxis1 := (axis + 2) % 3
			belowSA := 2 * (d.At(otherAxis0)*d.At(otherAxis1) + (edget-nodeBounds.PMin.At(axis))*(d.At(otherAxis0)+d.At(otherAxis1)))
			aboveSA := 2 * (d.At(otherAxis0)*d.At(otherAxis1) + (nodeBounds.PMax.At(axis)-edget)*(d.At(otherAxis0)+d.At(otherAxis1)))
			pBelow := belowSA * invTotalSA
			pAbove := aboveSA * invTotalSA
			eb := 0.0
			if nAbove == 0 || nBelow == 0 {
				eb = accel.emptyBonus
			}
			cost := float64(accel.traversalCost) + float64(accel.isectCost)*(1.0-eb)*(pBelow*float64(nBelow)+pAbove*float64(nAbove))

			// Update best split if this is lowest cost so far
			if cost < bestCost {
				bestCost = cost
				bestAxis = axis
				bestOffset = i
			}
		}
		if edges[axis][i].edgeType == EDGE_START {
			nBelow++
		}
	}
	Assert(nBelow == nPrimitives && nAbove == 0)

	// Create leaf if no good splits were found
	if bestAxis == -1 && retries < 2 {
		retries++
		axis = (axis + 1) % 3
		goto retrySplit
	}
	if bestCost > oldCost {
		badRefines++
	}
	if (bestCost > 4.0*oldCost && nPrimitives < 16) ||
		bestAxis == -1 || badRefines == 3 {
		//PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth-depth);
		accel.nodes[nodeNum].initLeaf(primNums, nPrimitives, accel.arena)
		return
	}

	// Classify primitives with respect to split
	n0, n1 := 0, 0
	for i := 0; i < bestOffset; i++ {
		if edges[bestAxis][i].edgeType == EDGE_START {
			prims0[n0] = edges[bestAxis][i].primNum
			n0++
		}
	}
	for i := bestOffset + 1; i < 2*nPrimitives; i++ {
		if edges[bestAxis][i].edgeType == EDGE_END {
			prims1[n1] = edges[bestAxis][i].primNum
			n1++
		}
	}

	// Recursively initialize children nodes
	tsplit := edges[bestAxis][bestOffset].t
	//PBRT_KDTREE_CREATED_INTERIOR_NODE(bestAxis, tsplit);
	bounds0 := nodeBounds
	bounds1 := nodeBounds
	bounds0.PMax.Set(bestAxis, tsplit)
	bounds1.PMin.Set(bestAxis, tsplit)
	accel.buildTree(nodeNum+1, bounds0,
		allPrimBounds, prims0, n0, depth-1, edges,
		prims0, prims1[nPrimitives:], badRefines)
	aboveChild := accel.nextFreeNode
	accel.nodes[nodeNum].initInterior(SplitAxis(bestAxis), aboveChild, tsplit)
	accel.buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1,
		depth-1, edges, prims0, prims1[nPrimitives:], badRefines)
}

func (node *KdAccelNode) initLeaf(primNums []int, numPrims int, arena *MemoryArena) {
	node.axis = SPLIT_LEAF
	node.nPrims = numPrims
	if numPrims == 0 {
		node.onePrimitive = 0
	} else if numPrims == 1 {
		node.onePrimitive = primNums[0]
	} else {
		node.primitives = make([]int, numPrims, numPrims)
		copy(node.primitives, primNums)
	}
}

func (node *KdAccelNode) initInterior(axis SplitAxis, aboveChild int, split float64) {
	node.split = split
	node.axis = axis
	node.aboveChild = aboveChild
}

func (b BoundEdge) String() string {
	return fmt.Sprintf("edge[t:%f]", b.t)
}

type boundEdgeSorter struct {
	axis  int
	nodes []BoundEdge
}

func (s *boundEdgeSorter) Len() int {
	return len(s.nodes)
}
func (s *boundEdgeSorter) Swap(i, j int) {
	s.nodes[i], s.nodes[j] = s.nodes[j], s.nodes[i]
}
func (s *boundEdgeSorter) Less(i, j int) bool {
	if s.nodes[i].t == s.nodes[j].t {
		return s.nodes[i].edgeType < s.nodes[j].edgeType
	}
	return s.nodes[i].t < s.nodes[j].t
}

func (node *KdAccelNode) String() string {
	if node.axis == SPLIT_LEAF {
		return fmt.Sprintf("kdnode[leaf np:%d]", node.nPrims)
	} else {
		return fmt.Sprintf("kdnode[axis:%s split:%f child:%d]", SplitAxisNames[node.axis], node.split, node.aboveChild)
	}
}
