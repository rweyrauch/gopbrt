package pbrt

import (
	"strings"
)

const (
	SPLIT_MIDDLE = iota
	SPLIT_EQUAL_COUNTS
	SPLIT_SAH
)
const (
	SAH_NUM_BUCKETS = 12
)

type SplitMethod int

type (
	SAHBucketInfo struct {
		count  int
		bounds BBox
	}

	BVHPrimitiveInfo struct {
		primitiveNumber int
		centroid        Point
		bounds          BBox
	}

	BVHBuildNode struct {
		bounds                                  BBox
		children                                [2]*BVHBuildNode
		splitAxis, firstPrimOffset, nPrimitives int
	}

	CompareToMid struct {
		dim int
		mid float64
	}
	ComparePoints struct {
		dim int
	}
	CompareToBucket struct {
		splitBucket, nBuckets, dim int
		centroidBounds             BBox
	}

	LinearBVHNode struct {
		bounds            BBox
		offset            int // primitive or second child
		nPrimitives, axis int
	}

	BVHAccel struct {
		PrimitiveData
		maxPrimsInNode int
		splitMethod    SplitMethod
		primitives     []Primitive
		nodes          []LinearBVHNode
	}
)

func CreateBVHPrimitiveInfo(pn int, b *BBox) BVHPrimitiveInfo {
	centroid := &Point{b.pMax.x + b.pMin.x, b.pMax.y + b.pMin.y, b.pMax.z + b.pMin.z}
	return BVHPrimitiveInfo{pn, *centroid.Scale(0.5), *b}
}

func CreateBVHBuildNode() *BVHBuildNode {
	node := new(BVHBuildNode)
	node.bounds = *CreateEmptyBBox()
	node.children[0] = nil
	node.children[1] = nil
	node.splitAxis = 0
	node.firstPrimOffset = 0
	node.nPrimitives = 0
	return node
}
func (node *BVHBuildNode) InitLeaf(first, n int, b *BBox) {
	node.firstPrimOffset = first
	node.nPrimitives = n
	node.bounds = *b
}
func (node *BVHBuildNode) InitInterior(axis int, c0, c1 *BVHBuildNode) {
	node.children[0] = c0
	node.children[1] = c1
	node.bounds = *UnionBBoxes(&c0.bounds, &c1.bounds)
	node.splitAxis = axis
	node.nPrimitives = 0
}

func IntersectP(bounds *BBox, ray *Ray, invDir *Vector, dirIsNeg [3]int) bool {
	// Check for ray intersection against $x$ and $y$ slabs
	tmin := (bounds.PointAtIndex(dirIsNeg[0]).x - ray.origin.x) * invDir.x
	tmax := (bounds.PointAtIndex(1-dirIsNeg[0]).x - ray.origin.x) * invDir.x
	tymin := (bounds.PointAtIndex(dirIsNeg[1]).y - ray.origin.y) * invDir.y
	tymax := (bounds.PointAtIndex(1-dirIsNeg[1]).y - ray.origin.y) * invDir.y
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
	tzmin := (bounds.PointAtIndex(dirIsNeg[2]).z - ray.origin.z) * invDir.z
	tzmax := (bounds.PointAtIndex(1-dirIsNeg[2]).z - ray.origin.z) * invDir.z
	if (tmin > tzmax) || (tzmin > tmax) {
		return false
	}
	if tzmin > tmin {
		tmin = tzmin
	}
	if tzmax < tmax {
		tmax = tzmax
	}
	return (tmin < ray.maxt) && (tmax > ray.mint)
}

func CreateBVHAccel(prims []Primitive, maxPrims int, sm string) *BVHAccel {
	bvh := new(BVHAccel)
	bvh.primitiveId = GeneratePrimitiveId()
	bvh.maxPrimsInNode = Mini(255, maxPrims)
	for _, p := range prims {
		bvh.primitives = p.FullyRefine(bvh.primitives)
	}
	if strings.Compare(sm, "shah") == 0 {
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
	//PBRT_BVH_STARTED_CONSTRUCTION(bvh, len(bvh.primitives))

	// Initialize _buildData_ array for primitives
	buildData := make([]BVHPrimitiveInfo, 0, len(bvh.primitives))
	for i, p := range bvh.primitives {
		bbox := p.WorldBound()
		buildData = append(buildData, CreateBVHPrimitiveInfo(i, bbox))
	}

	// Recursively build BVH tree for primitives
	var buildArena *MemoryArena = nil
	orderedPrims := make([]Primitive, 0, len(bvh.primitives))
	root, totalNodes := bvh.recursiveBuild(buildArena, buildData, 0, len(bvh.primitives), orderedPrims)
	bvh.primitives, orderedPrims = orderedPrims, bvh.primitives
	Info("BVH created with %d nodes for %d primitives", totalNodes, len(bvh.primitives))

	// Compute representation of depth-first traversal of BVH tree
	bvh.nodes = make([]LinearBVHNode, totalNodes, totalNodes)
	var offset int = 0
	offset, _ = bvh.flattenBVHTree(root, offset)
	if offset != totalNodes {
		Severe("Incorrect number of BVH detected.")
	}
	//PBRT_BVH_FINISHED_CONSTRUCTION(bvh)

	return bvh
}

func (bvh *BVHAccel) flattenBVHTree(node *BVHBuildNode, offset int) (myOffset, nextOffset int) {
	linearNode := &bvh.nodes[offset]
	linearNode.bounds = node.bounds
	myOffset = offset
	nextOffset = offset + 1
	if node.nPrimitives > 0 {
		//Assert(!node->children[0] && !node->children[1])
		linearNode.offset = node.firstPrimOffset
		linearNode.nPrimitives = node.nPrimitives
	} else {
		// Creater interior flattened BVH node
		linearNode.axis = node.splitAxis
		linearNode.nPrimitives = 0
		bvh.flattenBVHTree(node.children[0], nextOffset)
		linearNode.offset, nextOffset = bvh.flattenBVHTree(node.children[1], nextOffset)
	}
	return myOffset, offset
}

func (bvh *BVHAccel) WorldBound() *BBox {
	if bvh.nodes != nil {
		return &bvh.nodes[0].bounds
	} else {
		return CreateEmptyBBox()
	}
}

func partition(buildData []BVHPrimitiveInfo, first, last int, comp CompareToMid) int {
	return 0
}
func partitionB(buildData []BVHPrimitiveInfo, first, last int, comp CompareToBucket) int {
	return 0
}
func nth_element(buildData []BVHPrimitiveInfo, first, nth, last int, comp ComparePoints) {

}

func (bvh *BVHAccel) recursiveBuild(buildArena *MemoryArena, buildData []BVHPrimitiveInfo, start, end int, orderedPrims []Primitive) (node *BVHBuildNode, totalNodes int) {
	//Assert(start != end);
	totalNodes++
	node = CreateBVHBuildNode()
	// Compute bounds of all primitives in BVH node
	bbox := CreateEmptyBBox()
	for i := start; i < end; i++ {
		bbox = UnionBBoxes(bbox, &buildData[i].bounds)
	}

	nPrimitives := end - start
	if nPrimitives == 1 {
		// Create leaf _BVHBuildNode_
		firstPrimOffset := len(orderedPrims)
		for i := start; i < end; i++ {
			primNum := buildData[i].primitiveNumber
			orderedPrims = append(orderedPrims, bvh.primitives[primNum])
		}
		node.InitLeaf(firstPrimOffset, nPrimitives, bbox)
	} else {
		// Compute bound of primitive centroids, choose split dimension _dim_
		centroidBounds := CreateEmptyBBox()
		for i := start; i < end; i++ {
			centroidBounds = UnionBBoxPoint(centroidBounds, &buildData[i].centroid)
		}
		dim := centroidBounds.MaximumExtent()

		// Partition primitives into two sets and build children
		mid := (start + end) / 2
		if centroidBounds.pMax.At(dim) == centroidBounds.pMin.At(dim) {
			// If nPrimitives is no greater than maxPrimsInNode,
			// then all the nodes can be stored in a compact bvh node.
			if nPrimitives <= bvh.maxPrimsInNode {
				// Create leaf _BVHBuildNode_
				firstPrimOffset := len(orderedPrims)
				for i := start; i < end; i++ {
					primNum := buildData[i].primitiveNumber
					orderedPrims = append(orderedPrims, bvh.primitives[primNum])
				}
				node.InitLeaf(firstPrimOffset, nPrimitives, bbox)

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
			pmid := 0.5 * (centroidBounds.pMin.At(dim) + centroidBounds.pMax.At(dim))
			mid := partition(buildData, start, end-1, CompareToMid{dim, pmid})
			if mid == start || mid == end {
				// for lots of prims with large overlapping bounding boxes, this
				// may fail to partition; in that case don't break and fall through
				// to SPLIT_EQUAL_COUNTS
				mid = (start + end) / 2
				nth_element(buildData, start, mid, end-1, ComparePoints{dim})
			}
		case SPLIT_EQUAL_COUNTS:
			// Partition primitives into equally-sized subsets
			mid = (start + end) / 2
			nth_element(buildData, start, mid, end-1, ComparePoints{dim})
		case SPLIT_SAH:
		default:
			// Partition primitives using approximate SAH
			if nPrimitives <= 4 {
				// Partition primitives into equally-sized subsets
				mid = (start + end) / 2
				nth_element(buildData, start, mid, end-1, ComparePoints{dim})
			} else {
				// Allocate _BucketInfo_ for SAH partition buckets
				buckets := make([]SAHBucketInfo, SAH_NUM_BUCKETS, SAH_NUM_BUCKETS)

				// Initialize _BucketInfo_ for SAH partition buckets
				for i := start; i < end; i++ {
					b := int(float64(SAH_NUM_BUCKETS) *
						((buildData[i].centroid.At(dim) - centroidBounds.pMin.At(dim)) /
							(centroidBounds.pMax.At(dim) - centroidBounds.pMin.At(dim))))
					if b == SAH_NUM_BUCKETS {
						b = SAH_NUM_BUCKETS - 1
					}
					//Assert(b >= 0 && b < SAH_NUM_BUCKETS);
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
					mid = partitionB(buildData, start, end-1, CompareToBucket{minCostSplit, SAH_NUM_BUCKETS, dim, *centroidBounds})
				} else {
					// Create leaf _BVHBuildNode_
					firstPrimOffset := len(orderedPrims)
					for i := start; i < end; i++ {
						primNum := buildData[i].primitiveNumber
						orderedPrims = append(orderedPrims, bvh.primitives[primNum])
					}
					node.InitLeaf(firstPrimOffset, nPrimitives, bbox)
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

func (bvh *BVHAccel) Intersect(ray *Ray) (hit bool, isect *Intersection) {
	if bvh.nodes == nil {
		return false, nil
	}
	//PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray))
	hit = false
	invDir := CreateVector(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z)
	dirIsNeg := [3]int{0, 0, 0}
	if invDir.x < 0 {
		dirIsNeg[0] = 1
	}
	if invDir.y < 0 {
		dirIsNeg[1] = 1
	}
	if invDir.z < 0 {
		dirIsNeg[2] = 1
	}
	// Follow ray through BVH nodes to find primitive intersections
	todo := make([]int, 64, 64)
	todoOffset, nodeNum := 0, 0
	for {
		node := &bvh.nodes[nodeNum]
		// Check ray against BVH node
		if IntersectP(&node.bounds, ray, invDir, dirIsNeg) {
			if node.nPrimitives > 0 {
				// Intersect ray with primitives in leaf BVH node
				//PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node))
				for i := 0; i < node.nPrimitives; i++ {
					//PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()))
					if hit, isect = bvh.primitives[node.offset+i].Intersect(ray); hit {
						//PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()))
						hit = true
					} else {
						//PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()))
					}
				}
				if todoOffset == 0 {
					break
				}
				todoOffset--
				nodeNum = todo[todoOffset]
			} else {
				// Put far BVH node on _todo_ stack, advance to near node
				//PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node))
				if dirIsNeg[node.axis] != 0 {
					todoOffset++
					todo[todoOffset] = nodeNum + 1
					nodeNum = node.offset
				} else {
					todoOffset++
					todo[todoOffset] = node.offset
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
	//PBRT_BVH_INTERSECTION_FINISHED()
	return hit, isect
}

func (bvh *BVHAccel) IntersectP(ray *Ray) bool {
	if bvh.nodes == nil {
		return false
	}
	//PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray))
	invDir := CreateVector(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z)
	dirIsNeg := [3]int{0, 0, 0}
	if invDir.x < 0 {
		dirIsNeg[0] = 1
	}
	if invDir.y < 0 {
		dirIsNeg[1] = 1
	}
	if invDir.z < 0 {
		dirIsNeg[2] = 1
	}
	todo := make([]int, 64, 64)
	todoOffset, nodeNum := 0, 0
	for {
		node := &bvh.nodes[nodeNum]
		if IntersectP(&node.bounds, ray, invDir, dirIsNeg) {
			// Process BVH node _node_ for traversal
			if node.nPrimitives > 0 {
				//PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node))
				for i := 0; i < node.nPrimitives; i++ {
					//PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()))
					if bvh.primitives[node.offset+i].IntersectP(ray) {
						//PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()))
						return true
					} else {
						//PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()))
					}
				}
				if todoOffset == 0 {
					break
				}
				todoOffset--
				nodeNum = todo[todoOffset]
			} else {
				//PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node))
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
	//PBRT_BVH_INTERSECTIONP_FINISHED()
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
	return CreateBVHAccel(prims, maxPrimsInNode, splitMethod)
}
