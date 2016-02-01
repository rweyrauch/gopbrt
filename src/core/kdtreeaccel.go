package core

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
		flags      SplitAxis
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
		node       *KdAccelNode
		tmin, tmax float64
	}
)

func (p *KdTreeAccel) WorldBound() *BBox                                  { return nil }
func (p *KdTreeAccel) CanIntersect() bool                                 { return false }
func (p *KdTreeAccel) Intersect(r *RayDifferential) (bool, *Intersection) { return false, nil }
func (p *KdTreeAccel) IntersectP(r *Ray) bool                             { return false }
func (p *KdTreeAccel) Refine(refined []Primitive) []Primitive             { return refined }
func (p *KdTreeAccel) FullyRefine(refined []Primitive) []Primitive        { return refined }
func (p *KdTreeAccel) GetAreaLight() AreaLight                            { return nil }
func (p *KdTreeAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (p *KdTreeAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (p *KdTreeAccel) PrimitiveId() uint32 { return p.primitiveId }

func CreateKdTreeAccelerator(prims []Primitive, ps *ParamSet) *KdTreeAccel { return nil }

func (accel *KdTreeAccel) buildTree(nodeNum int, nodeBounds *BBox, allPrimBounds []BBox, primNums []int, nPrimitives, depth int, edges [3][]BoundEdge, prims0, prims1 []int, badRefines int) {
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
	d := nodeBounds.pMax.Sub(&nodeBounds.pMin)

	// Choose which axis to split along
	axis := nodeBounds.MaximumExtent()
	retries := 0

retrySplit:

	// Initialize edges for _axis_
	for i := 0; i < nPrimitives; i++ {
		pn := primNums[i]
		bbox := allPrimBounds[pn]
		edges[axis][2*i] = BoundEdge{bbox.pMin.At(axis), pn, EDGE_START}
		edges[axis][2*i+1] = BoundEdge{bbox.pMax.At(axis), pn, EDGE_END}
	}
	// TODO sort edges by axis
	//sort(&edges[axis][0], &edges[axis][2*nPrimitives]);

	// Compute cost of all splits for _axis_ to find best
	nBelow, nAbove := 0, nPrimitives
	for i := 0; i < 2*nPrimitives; i++ {
		if edges[axis][i].edgeType == EDGE_END {
			nAbove--
		}
		edget := edges[axis][i].t
		if edget > nodeBounds.pMin.At(axis) && edget < nodeBounds.pMax.At(axis) {
			// Compute cost for split at _i_th edge
			otherAxis0 := (axis + 1) % 3
			otherAxis1 := (axis + 2) % 3
			belowSA := 2 * (d.At(otherAxis0)*d.At(otherAxis1) + (edget-nodeBounds.pMin.At(axis))*(d.At(otherAxis0)+d.At(otherAxis1)))
			aboveSA := 2 * (d.At(otherAxis0)*d.At(otherAxis1) + (nodeBounds.pMax.At(axis)-edget)*(d.At(otherAxis0)+d.At(otherAxis1)))
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
	bounds0.pMax.Set(bestAxis, tsplit)
	bounds1.pMin.Set(bestAxis, tsplit)
	accel.buildTree(nodeNum+1, bounds0,
		allPrimBounds, prims0, n0, depth-1, edges,
		prims0, prims1[nPrimitives:], badRefines)
	aboveChild := accel.nextFreeNode
	accel.nodes[nodeNum].initInterior(SplitAxis(bestAxis), aboveChild, tsplit)
	accel.buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1,
		depth-1, edges, prims0, prims1[nPrimitives:], badRefines)
}

func (node *KdAccelNode) initLeaf(primNums []int, numPrims int, arena *MemoryArena) {
	node.flags = SPLIT_LEAF
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
	node.flags = axis
	node.aboveChild = aboveChild
}
