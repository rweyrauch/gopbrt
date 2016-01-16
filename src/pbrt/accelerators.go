package pbrt

const (
	SPLIT_MIDDLE = iota
	SPLIT_EQUAL_COUNTS
	SPLIT_SAH
)

type SplitMethod int

type (
	Voxel         struct{}
	KdAccelNode   struct{}
	LinearBVHNode struct{}

	BVHAccel struct {
		PrimitiveData
		maxPrimsInNode int
		splitMethod    SplitMethod
		primitives     []Primitive
		nodes          []LinearBVHNode
	}

	GridAccel struct {
		PrimitiveData
		primitives      []Primitive
		nVoxels         [3]int
		bounds          BBox
		width, invWidth Vector
		voxels          [][]Voxel
		voxelArena      *MemoryArena
	}

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
)

func CreateBVHAccel(prims []Primitive) *BVHAccel                { return nil }
func (p *BVHAccel) WorldBound() *BBox                           { return nil }
func (p *BVHAccel) CanIntersect() bool                          { return false }
func (p *BVHAccel) Intersect(r *Ray) (bool, *Intersection)      { return false, nil }
func (p *BVHAccel) IntersectP(r *Ray) bool                      { return false }
func (p *BVHAccel) Refine(refined []Primitive) []Primitive      { return refined }
func (p *BVHAccel) FullyRefine(refined []Primitive) []Primitive { return refined }
func (p *BVHAccel) GetAreaLight() AreaLight                     { return nil }
func (p *BVHAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (p *BVHAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (p *BVHAccel) PrimitiveId() uint32 { return p.primitiveId }

func (p *GridAccel) WorldBound() *BBox                           { return nil }
func (p *GridAccel) CanIntersect() bool                          { return false }
func (p *GridAccel) Intersect(r *Ray) (bool, *Intersection)      { return false, nil }
func (p *GridAccel) IntersectP(r *Ray) bool                      { return false }
func (p *GridAccel) Refine(refined []Primitive) []Primitive      { return refined }
func (p *GridAccel) FullyRefine(refined []Primitive) []Primitive { return refined }
func (p *GridAccel) GetAreaLight() AreaLight                     { return nil }
func (p *GridAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (p *GridAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (p *GridAccel) PrimitiveId() uint32 { return p.primitiveId }

func (p *KdTreeAccel) WorldBound() *BBox                           { return nil }
func (p *KdTreeAccel) CanIntersect() bool                          { return false }
func (p *KdTreeAccel) Intersect(r *Ray) (bool, *Intersection)      { return false, nil }
func (p *KdTreeAccel) IntersectP(r *Ray) bool                      { return false }
func (p *KdTreeAccel) Refine(refined []Primitive) []Primitive      { return refined }
func (p *KdTreeAccel) FullyRefine(refined []Primitive) []Primitive { return refined }
func (p *KdTreeAccel) GetAreaLight() AreaLight                     { return nil }
func (p *KdTreeAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (p *KdTreeAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (p *KdTreeAccel) PrimitiveId() uint32 { return p.primitiveId }

func CreateBVHAccelerator(prims []Primitive, ps *ParamSet) *BVHAccel       { return nil }
func CreateGridAccelerator(prims []Primitive, ps *ParamSet) *GridAccel     { return nil }
func CreateKdTreeAccelerator(prims []Primitive, ps *ParamSet) *KdTreeAccel { return nil }
