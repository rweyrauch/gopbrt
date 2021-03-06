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
	"sort"
	"strings"
)

const (
	SPLIT_MIDDLE SplitMethod = iota
	SPLIT_EQUAL_COUNTS
	SPLIT_SAH
)
const (
	SAH_NUM_BUCKETS = 12
)

type SplitMethod int

type (
	BVHAccel struct {
		PrimitiveData
		maxPrimsInNode int
		splitMethod    SplitMethod
		primitives     []Primitive
		nodes          []linearBVHNode
	}

	sahBucketInfo struct {
		count  int
		bounds BBox
	}

	bvhPrimitiveInfo struct {
		primitiveNumber int
		centroid        Point
		bounds          BBox
	}

	bvhBuildNode struct {
		bounds                                  BBox
		children                                [2]*bvhBuildNode
		splitAxis, firstPrimOffset, nPrimitives int
	}

	linearBVHNode struct {
		bounds            BBox
		offset            int // primitive or second child
		nPrimitives, axis int
	}
)

func createBVHPrimitiveInfo(pn int, b *BBox) bvhPrimitiveInfo {
	centroid := &Point{b.PMax.X + b.PMin.X, b.PMax.Y + b.PMin.Y, b.PMax.Z + b.PMin.Z}
	return bvhPrimitiveInfo{pn, *centroid.Scale(0.5), *b}
}

func newBVHBuildNode() *bvhBuildNode {
	node := new(bvhBuildNode)
	node.bounds = *CreateEmptyBBox()
	node.children[0] = nil
	node.children[1] = nil
	node.splitAxis = 0
	node.firstPrimOffset = 0
	node.nPrimitives = 0
	return node
}
func (node *bvhBuildNode) initLeaf(first, n int, b *BBox) {
	node.firstPrimOffset = first
	node.nPrimitives = n
	node.bounds = *b
}
func (node *bvhBuildNode) InitInterior(axis int, c0, c1 *bvhBuildNode) {
	node.children[0] = c0
	node.children[1] = c1
	node.bounds = *UnionBBoxes(&c0.bounds, &c1.bounds)
	node.splitAxis = axis
	node.nPrimitives = 0
}

func intersectP(bounds *BBox, ray RayBase, invDir *Vector, dirIsNeg *[3]int) bool {
	// Check for ray intersection against $x$ and $y$ slabs
	tmin := (bounds.PointAtIndex(dirIsNeg[0]).X - ray.Origin().X) * invDir.X
	tmax := (bounds.PointAtIndex(1-dirIsNeg[0]).X - ray.Origin().X) * invDir.X
	tymin := (bounds.PointAtIndex(dirIsNeg[1]).Y - ray.Origin().Y) * invDir.Y
	tymax := (bounds.PointAtIndex(1-dirIsNeg[1]).Y - ray.Origin().Y) * invDir.Y
	if (tmin > tymax) || (tymin > tmax) {
		return false
	}
	if tymin > tmin {
		tmin = tymin
	}
	if tymax < tmax {
		tmax = tymax
	}

	// Check for ray intersection against $z$ slab
	tzmin := (bounds.PointAtIndex(dirIsNeg[2]).Z - ray.Origin().Z) * invDir.Z
	tzmax := (bounds.PointAtIndex(1-dirIsNeg[2]).Z - ray.Origin().Z) * invDir.Z
	if (tmin > tzmax) || (tzmin > tmax) {
		return false
	}
	if tzmin > tmin {
		tmin = tzmin
	}
	if tzmax < tmax {
		tmax = tzmax
	}
	return (tmin < ray.Maxt()) && (tmax > ray.Mint())
}

func NewBVHAccel(prims []Primitive, maxPrims int, sm string) *BVHAccel {
	bvh := new(BVHAccel)
	bvh.primitiveId = GeneratePrimitiveId()
	bvh.maxPrimsInNode = Mini(255, maxPrims)
	bvh.primitives = make([]Primitive, 0, len(prims))
	for _, p := range prims {
		bvh.primitives = p.FullyRefine(bvh.primitives)
	}
	if strings.Compare(sm, "sah") == 0 {
		bvh.splitMethod = SPLIT_SAH
	} else if strings.Compare(sm, "middle") == 0 {
		bvh.splitMethod = SPLIT_MIDDLE
	} else if strings.Compare(sm, "equal") == 0 {
		bvh.splitMethod = SPLIT_EQUAL_COUNTS
	} else {
		Warning("BVH split method \"%s\" unknown.  Using \"sah\".", sm)
		bvh.splitMethod = SPLIT_SAH
	}

	// early out when no primitives
	if len(bvh.primitives) == 0 {
		bvh.nodes = nil
		return bvh
	}

	// Build BVH from _primitives_

	// Initialize _buildData_ array for primitives
	buildData := make([]bvhPrimitiveInfo, 0, len(bvh.primitives))
	for i, p := range bvh.primitives {
		bbox := p.WorldBound()
		buildData = append(buildData, createBVHPrimitiveInfo(i, bbox))
	}

	// Recursively build BVH tree for primitives
	var buildArena *MemoryArena = nil
	orderedPrims := make([]Primitive, 0, len(bvh.primitives))
	root, totalNodes := bvh.recursiveBuild(buildArena, buildData, 0, len(bvh.primitives), &orderedPrims)
	Assert(len(bvh.primitives) == len(orderedPrims))
	copy(bvh.primitives, orderedPrims)
	Debug("BVH created with %d nodes for %d primitives", totalNodes, len(bvh.primitives))

	// Compute representation of depth-first traversal of BVH tree
	bvh.nodes = make([]linearBVHNode, totalNodes, totalNodes)
	var offset int = 0
	_, offset = bvh.flattenBVHTree(root, offset)
	if offset != totalNodes {
		Severe("Incorrect number of BVH detected.  Expected: %d  Got: %d", totalNodes, offset)
	}

	return bvh
}

func (bvh *BVHAccel) flattenBVHTree(node *bvhBuildNode, offset int) (myOffset, nextOffset int) {
	linearNode := &bvh.nodes[offset]
	linearNode.bounds = node.bounds
	myOffset = offset
	nextOffset = offset + 1
	if node.nPrimitives > 0 {
		Assert(node.children[0] == nil && node.children[1] == nil)
		linearNode.offset = node.firstPrimOffset
		linearNode.nPrimitives = node.nPrimitives
	} else {
		// Create interior flattened BVH node
		linearNode.axis = node.splitAxis
		linearNode.nPrimitives = 0
		_, nextOffset = bvh.flattenBVHTree(node.children[0], nextOffset)
		linearNode.offset, nextOffset = bvh.flattenBVHTree(node.children[1], nextOffset)
	}
	return myOffset, nextOffset
}

func (bvh *BVHAccel) WorldBound() *BBox {
	if bvh.nodes != nil {
		return &bvh.nodes[0].bounds
	} else {
		return CreateEmptyBBox()
	}
}

func partition(buildData []bvhPrimitiveInfo, first, last int, comp func(info *bvhPrimitiveInfo) bool) int {
	if first == last {
		return first
	}
	first++
	part := first
	if first == last {
		if comp(&buildData[part]) {
			return first
		} else {
			return part
		}
	}
	for first != last {
		if comp(&buildData[part]) {
			part++
		} else if comp(&buildData[first]) {
			buildData[part], buildData[first] = buildData[first], buildData[part]
			part++
		}
		first++
	}
	return part
}

type bvhInfoSorter struct {
	buildData []bvhPrimitiveInfo
	by        func(info0, info1 *bvhPrimitiveInfo) bool
}

func (s *bvhInfoSorter) Len() int {
	return len(s.buildData)
}
func (s *bvhInfoSorter) Swap(i, j int) {
	s.buildData[i], s.buildData[j] = s.buildData[j], s.buildData[i]
}
func (s *bvhInfoSorter) Less(i, j int) bool {
	return s.by(&s.buildData[i], &s.buildData[j])
}

type By func(i0, i1 *bvhPrimitiveInfo) bool

func (by By) Sort(buildData []bvhPrimitiveInfo) {
	infoSorter := &bvhInfoSorter{buildData, by}
	sort.Sort(infoSorter)
}

// TODO: implement a real nth_element function rather than a full sort on the range
func nth_element(buildData []bvhPrimitiveInfo, first, nth, last int, comp func(info0, info1 *bvhPrimitiveInfo) bool) {
	By(comp).Sort(buildData[first:last])
}

func (bvh *BVHAccel) recursiveBuild(buildArena *MemoryArena, buildData []bvhPrimitiveInfo, start, end int, orderedPrims *[]Primitive) (node *bvhBuildNode, totalNodes int) {
	Assert(start != end)
	totalNodes++
	node = newBVHBuildNode()
	// Compute bounds of all primitives in BVH node
	bbox := CreateEmptyBBox()
	for i := start; i < end; i++ {
		bbox = UnionBBoxes(bbox, &buildData[i].bounds)
	}

	nPrimitives := end - start
	if nPrimitives == 1 {
		// Create leaf _BVHBuildNode_
		firstPrimOffset := len(*orderedPrims)
		for i := start; i < end; i++ {
			primNum := buildData[i].primitiveNumber
			*orderedPrims = append(*orderedPrims, bvh.primitives[primNum])
		}
		node.initLeaf(firstPrimOffset, nPrimitives, bbox)
	} else {
		// Compute bound of primitive centroids, choose split dimension _dim_
		centroidBounds := CreateEmptyBBox()
		for i := start; i < end; i++ {
			centroidBounds = UnionBBoxPoint(centroidBounds, &buildData[i].centroid)
		}
		dim := centroidBounds.MaximumExtent()

		// Partition primitives into two sets and build children
		mid := (start + end) / 2
		if centroidBounds.PMax.At(dim) == centroidBounds.PMin.At(dim) {
			// If nPrimitives is no greater than maxPrimsInNode,
			// then all the nodes can be stored in a compact bvh node.
			if nPrimitives <= bvh.maxPrimsInNode {
				// Create leaf _BVHBuildNode_
				firstPrimOffset := len(*orderedPrims)
				for i := start; i < end; i++ {
					primNum := buildData[i].primitiveNumber
					*orderedPrims = append(*orderedPrims, bvh.primitives[primNum])
				}
				node.initLeaf(firstPrimOffset, nPrimitives, bbox)

				return node, totalNodes
			} else {
				// else if nPrimitives is greater than maxPrimsInNode, we
				// need to split it further to guarantee each node contains
				// no more than maxPrimsInNode primitives.
				c0, totalC0Nodes := bvh.recursiveBuild(buildArena, buildData, start, mid, orderedPrims)
				totalNodes += totalC0Nodes
				c1, totalC1Nodes := bvh.recursiveBuild(buildArena, buildData, mid, end, orderedPrims)
				totalNodes += totalC1Nodes
				node.InitInterior(dim, c0, c1)

				return node, totalNodes
			}
		}

		// Partition primitives based on _splitMethod_
		switch bvh.splitMethod {
		case SPLIT_MIDDLE:
			// Partition primitives through node's midpoint
			pmid := 0.5 * (centroidBounds.PMin.At(dim) + centroidBounds.PMax.At(dim))
			compareToMid := func(info *bvhPrimitiveInfo) bool {
				if info.centroid.At(dim) < pmid {
					return true
				} else {
					return false
				}
			}
			mid := partition(buildData, start, end-1, compareToMid)
			if mid == start || mid == end {
				// for lots of prims with large overlapping bounding boxes, this
				// may fail to partition; in that case don't break and fall through
				// to SPLIT_EQUAL_COUNTS
				mid = (start + end) / 2
				comparePoints := func(info0, info1 *bvhPrimitiveInfo) bool {
					return info0.centroid.At(dim) < info0.centroid.At(dim)
				}
				nth_element(buildData, start, mid, end-1, comparePoints)
			}
		case SPLIT_EQUAL_COUNTS:
			// Partition primitives into equally-sized subsets
			mid = (start + end) / 2
			comparePoints := func(info0, info1 *bvhPrimitiveInfo) bool {
				return info0.centroid.At(dim) < info0.centroid.At(dim)
			}
			nth_element(buildData, start, mid, end-1, comparePoints)
		case SPLIT_SAH:
		default:
			// Partition primitives using approximate SAH
			if nPrimitives <= 4 {
				// Partition primitives into equally-sized subsets
				mid = (start + end) / 2
				comparePoints := func(info0, info1 *bvhPrimitiveInfo) bool {
					return info0.centroid.At(dim) < info0.centroid.At(dim)
				}
				nth_element(buildData, start, mid, end-1, comparePoints)
			} else {
				// Allocate _BucketInfo_ for SAH partition buckets
				buckets := make([]sahBucketInfo, SAH_NUM_BUCKETS, SAH_NUM_BUCKETS)

				// Initialize _BucketInfo_ for SAH partition buckets
				for i := start; i < end; i++ {
					b := int(float64(SAH_NUM_BUCKETS) *
						((buildData[i].centroid.At(dim) - centroidBounds.PMin.At(dim)) /
							(centroidBounds.PMax.At(dim) - centroidBounds.PMin.At(dim))))
					if b == SAH_NUM_BUCKETS {
						b = SAH_NUM_BUCKETS - 1
					}
					Assert(b >= 0 && b < SAH_NUM_BUCKETS)
					buckets[b].count++
					buckets[b].bounds = *UnionBBoxes(&buckets[b].bounds, &buildData[i].bounds)
				}

				// Compute costs for splitting after each bucket
				cost := make([]float64, SAH_NUM_BUCKETS-1, SAH_NUM_BUCKETS-1)
				for i := 0; i < len(cost); i++ {
					var b0, b1 *BBox
					count0 := 0
					count1 := 0
					for j := 0; j <= i; j++ {
						b0 = UnionBBoxes(b0, &buckets[j].bounds)
						count0 += buckets[j].count
					}
					for j := i + 1; j < SAH_NUM_BUCKETS; j++ {
						b1 = UnionBBoxes(b1, &buckets[j].bounds)
						count1 += buckets[j].count
					}
					cost[i] = 0.125 + (float64(count0)*b0.SurfaceArea()+float64(count1)*b1.SurfaceArea())/bbox.SurfaceArea()
				}

				// Find bucket to split at that minimizes SAH metric
				minCost := cost[0]
				minCostSplit := 0
				for i := 1; i < len(cost); i++ {
					if cost[i] < minCost {
						minCost = cost[i]
						minCostSplit = i
					}
				}

				// Either create leaf or split primitives at selected SAH bucket
				if nPrimitives > bvh.maxPrimsInNode || minCost < float64(nPrimitives) {
					compareToBucket := func(info *bvhPrimitiveInfo) bool {
						b := int(float64(SAH_NUM_BUCKETS) * ((info.centroid.At(dim) - centroidBounds.PMin.At(dim)) / (centroidBounds.PMax.At(dim) - centroidBounds.PMin.At(dim))))
						if b == SAH_NUM_BUCKETS {
							b = SAH_NUM_BUCKETS - 1
						}
						Assert(b >= 0 && b < SAH_NUM_BUCKETS)
						return b <= minCostSplit
					}
					mid = partition(buildData, start, end-1, compareToBucket)
				} else {
					// Create leaf _BVHBuildNode_
					firstPrimOffset := len(*orderedPrims)
					for i := start; i < end; i++ {
						primNum := buildData[i].primitiveNumber
						*orderedPrims = append(*orderedPrims, bvh.primitives[primNum])
					}
					node.initLeaf(firstPrimOffset, nPrimitives, bbox)
					return node, totalNodes
				}
			}
		}
		c0, totalC0Nodes := bvh.recursiveBuild(buildArena, buildData, start, mid, orderedPrims)
		totalNodes += totalC0Nodes
		c1, totalC1Nodes := bvh.recursiveBuild(buildArena, buildData, mid, end, orderedPrims)
		totalNodes += totalC1Nodes
		node.InitInterior(dim, c0, c1)
	}
	return node, totalNodes
}

func (p *BVHAccel) CanIntersect() bool { return true }

func (bvh *BVHAccel) Intersect(ray RayBase) (bool, *Intersection) {
	if bvh.nodes == nil {
		return false, nil
	}
	hitSomething := false
	var isect *Intersection

	invDir := CreateVector(1.0/ray.Dir().X, 1.0/ray.Dir().Y, 1.0/ray.Dir().Z)
	dirIsNeg := [3]int{0, 0, 0}
	if invDir.X < 0 {
		dirIsNeg[0] = 1
	}
	if invDir.Y < 0 {
		dirIsNeg[1] = 1
	}
	if invDir.Z < 0 {
		dirIsNeg[2] = 1
	}
	// Follow ray through BVH nodes to find primitive intersections
	todo := make([]int, 64, 64)
	todoOffset, nodeNum := 0, 0
	for {
		node := &bvh.nodes[nodeNum]
		// Check ray against BVH node
		if intersectP(&node.bounds, ray, invDir, &dirIsNeg) {
			if node.nPrimitives > 0 {
				// Intersect ray with primitives in leaf BVH node
				for i := 0; i < node.nPrimitives; i++ {
					if hit, nisect := bvh.primitives[node.offset+i].Intersect(ray); hit {
						hitSomething = true
						isect = nisect
					} else {
					}
				}
				if todoOffset == 0 {
					break
				}
				todoOffset--
				nodeNum = todo[todoOffset]
			} else {
				// Put far BVH node on _todo_ stack, advance to near node
				if dirIsNeg[node.axis] != 0 {
					todo[todoOffset] = nodeNum + 1
					todoOffset++
					nodeNum = node.offset
				} else {
					todo[todoOffset] = node.offset
					todoOffset++
					nodeNum = nodeNum + 1
				}
			}
		} else {
			if todoOffset == 0 {
				break
			}
			todoOffset--
			nodeNum = todo[todoOffset]
		}
	}
	return hitSomething, isect
}

func (bvh *BVHAccel) IntersectP(ray RayBase) bool {
	if bvh.nodes == nil {
		return false
	}
	invDir := CreateVector(1.0/ray.Dir().X, 1.0/ray.Dir().Y, 1.0/ray.Dir().Z)
	dirIsNeg := [3]int{0, 0, 0}
	if invDir.X < 0 {
		dirIsNeg[0] = 1
	}
	if invDir.Y < 0 {
		dirIsNeg[1] = 1
	}
	if invDir.Z < 0 {
		dirIsNeg[2] = 1
	}
	todo := make([]int, 64, 64)
	todoOffset, nodeNum := 0, 0
	for {
		node := &bvh.nodes[nodeNum]
		if intersectP(&node.bounds, ray, invDir, &dirIsNeg) {
			// Process BVH node _node_ for traversal
			if node.nPrimitives > 0 {
				for i := 0; i < node.nPrimitives; i++ {
					if bvh.primitives[node.offset+i].IntersectP(ray) {
						return true
					} else {
					}
				}
				if todoOffset == 0 {
					break
				}
				todoOffset--
				nodeNum = todo[todoOffset]
			} else {
				if dirIsNeg[node.axis] != 0 {
					/// second child first
					todo[todoOffset] = nodeNum + 1
					todoOffset++
					nodeNum = node.offset
				} else {
					todo[todoOffset] = node.offset
					todoOffset++
					nodeNum = nodeNum + 1
				}
			}
		} else {
			if todoOffset == 0 {
				break
			}
			todoOffset--
			nodeNum = todo[todoOffset]
		}
	}
	return false
}

func (bvh *BVHAccel) Refine(refined []Primitive) []Primitive {
	Severe("Unimplemented BVHAccel::Refine() method called!")
	return nil
}

func (bvh *BVHAccel) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(bvh, refined)
}

func (p *BVHAccel) GetAreaLight() AreaLight {
	Severe("BVHAccel::GetAreaLight() method called; should have gone to GeometricPrimitive")
	return nil
}

func (p *BVHAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	Severe("BVHAccel::GetBSDF() method called; should have gone to GeometricPrimitive")
	return nil
}

func (p *BVHAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	Severe("BVHAccel::GetBSSRDF() method called; should have gone to GeometricPrimitive")
	return nil
}

func (bvh *BVHAccel) PrimitiveId() uint32 { return bvh.primitiveId }

func CreateBVHAccelerator(prims []Primitive, ps *ParamSet) *BVHAccel {
	splitMethod := ps.FindStringParam("splitmethod", "sah")
	maxPrimsInNode := ps.FindIntParam("maxnodeprims", 4)
	return NewBVHAccel(prims, maxPrimsInNode, splitMethod)
}
