package pbrt

var (
	nextPrimitiveId uint32 = 0
)

func GeneratePrimitiveId() uint32 {
	nextPrimitiveId++
	return nextPrimitiveId
}

type Primitive interface {
	WorldBound() *BBox
	CanIntersect() bool
	Intersect(r *Ray) (bool, *Intersection)
	IntersectP(r *Ray) bool
	Refine() (refined []Primitive)
	FullyRefine() (refined []Primitive)
	GetAreaLight() *AreaLight
	GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF
	GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF
	PrimitiveId() uint32	
}

type PrimitiveData struct {
	primitiveId uint32
}

type GeometricPrimitive struct {
	PrimitiveData
	shape     Shape
	material  Material
	areaLight *AreaLight
}

func (p *GeometricPrimitive) WorldBound() *BBox {
	return p.shape.WorldBound()
}

func (p *GeometricPrimitive) CanIntersect() bool {
	return p.shape.CanIntersect()
}

func (p *GeometricPrimitive) Intersect(r *Ray) (hit bool, isect *Intersection) {
    var thit, rayEpsilon float64
    
    if hit, thit, rayEpsilon, _ = p.shape.Intersect(r); !hit {
        return false, nil
    }
    
    isect = new(Intersection)
    isect.primitive = p
    isect.WorldToObject = p.shape.WorldToObject()
    isect.ObjectToWorld = p.shape.ObjectToWorld()
    isect.shapeId = p.shape.ShapeId()
    isect.primitiveId = p.primitiveId
    isect.rayEpsilon = rayEpsilon
    r.maxt = thit
    
    return true, isect
}

func (p *GeometricPrimitive) IntersectP(r *Ray) bool {
	return p.shape.IntersectP(r)
}

func (p *GeometricPrimitive) Refine() (refined []Primitive) {
	return nil
}

func (p *GeometricPrimitive) FullyRefine() (refined []Primitive) {
	return nil
}

func (p *GeometricPrimitive) GetAreaLight() *AreaLight {
	return p.areaLight
}

func (p *GeometricPrimitive) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
    dgs := p.shape.GetShadingGeometry(objectToWorld, dg)
    return p.material.GetBSDF(dg, dgs, arena)
}

func (p *GeometricPrimitive) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
    dgs := p.shape.GetShadingGeometry(objectToWorld, dg)
    return p.material.GetBSSRDF(dg, dgs, arena)
}

func (p *GeometricPrimitive) PrimitiveId() uint32 {
	return p.primitiveId
}

type TransformedPrimitive struct {
	PrimitiveData
	primitive     Primitive
	worldToPrimitive  *AnimatedTransform
}

func (p *TransformedPrimitive) WorldBound() *BBox {
	return nil
}

func (p *TransformedPrimitive) CanIntersect() bool {
	return true
}

func (p *TransformedPrimitive) Intersect(r *Ray) (hit bool, isect *Intersection) {
    w2p := p.worldToPrimitive.Interpolate(r.time)
    ray := RayTransform(w2p, r)
    if hit, isect = p.primitive.Intersect(ray); !hit {
        return false, nil
    }
    
    r.maxt = ray.maxt
    isect.primitiveId = p.primitiveId
    if !IsIdentityTransform(w2p) {
        // Compute world-to-object transformation for instance
        isect.WorldToObject = isect.WorldToObject.MultTransform(w2p)
        isect.ObjectToWorld = InverseTransform(isect.WorldToObject)

        // Transform instance's differential geometry to world space
        primitiveToWorld := InverseTransform(w2p)
        isect.dg.p = PointTransform(primitiveToWorld, isect.dg.p)
        isect.dg.nn = NormalizeNormal(NormalTransform(primitiveToWorld, isect.dg.nn))
        isect.dg.dpdu = VectorTransform(primitiveToWorld, isect.dg.dpdu)
        isect.dg.dpdv = VectorTransform(primitiveToWorld, isect.dg.dpdv)
        isect.dg.dndu = NormalTransform(primitiveToWorld, isect.dg.dndu)
        isect.dg.dndv = NormalTransform(primitiveToWorld, isect.dg.dndv)
    }
    return true, isect
}

func (p *TransformedPrimitive) IntersectP(r *Ray) bool {
	return p.primitive.IntersectP(RayAnimatedTransform(p.worldToPrimitive, r))
}

func (p *TransformedPrimitive) Refine() (refined []Primitive) {
	return nil
}

func (p *TransformedPrimitive) FullyRefine() (refined []Primitive) {
	return nil
}

func (p *TransformedPrimitive) GetAreaLight() *AreaLight {
	return nil
}

func (p *TransformedPrimitive) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}

func (p *TransformedPrimitive) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}

func (p *TransformedPrimitive) PrimitiveId() uint32 {
	return p.primitiveId
}

type Aggregate struct {
	PrimitiveData
}

func (p *Aggregate) WorldBound() *BBox {
	return nil
}

func (p *Aggregate) CanIntersect() bool {
	return true
}

func (p *Aggregate) Intersect(r *Ray) (bool, *Intersection) {
	return false, nil
}

func (p *Aggregate) IntersectP(r *Ray) bool {
	return true
}

func (p *Aggregate) Refine() (refined []Primitive) {
	return nil
}

func (p *Aggregate) FullyRefine() (refined []Primitive) {
	return nil
}

func (p *Aggregate) GetAreaLight() *AreaLight {
	return nil
}

func (p *Aggregate) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}

func (p *Aggregate) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}

func (p *Aggregate) PrimitiveId() uint32 {
	return p.primitiveId
}
