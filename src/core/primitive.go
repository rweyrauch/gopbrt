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
	Intersect(r RayBase) (bool, *Intersection)
	IntersectP(r RayBase) bool
	Refine(refined []Primitive) []Primitive
	FullyRefine(refined []Primitive) []Primitive
	GetAreaLight() AreaLight
	GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF
	GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF
	PrimitiveId() uint32
}

func PrimitiveFullyRefine(p Primitive, refined []Primitive) []Primitive {
	todo := make([]Primitive, 0, 8)
	todo = append(todo, p)

	for len(todo) != 0 {
		// Refine last primitive in todo list
		prim := todo[len(todo)-1]
		todo = todo[:len(todo)-1] // pop off of stack

		if prim.CanIntersect() {
			if refined == nil {
				refined = make([]Primitive, 0, 2)
			}
			refined = append(refined, prim)
		} else {
			todo = prim.Refine(todo)
		}
	}
	return refined
}

type PrimitiveData struct {
	primitiveId uint32
}

type GeometricPrimitive struct {
	PrimitiveData
	shape     Shape
	material  Material
	areaLight AreaLight
}

func NewGeometricPrimitive(s Shape, mtl Material, arealight AreaLight) *GeometricPrimitive {
	gp := new(GeometricPrimitive)
	gp.primitiveId = GeneratePrimitiveId()
	gp.shape = s
	gp.material = mtl
	gp.areaLight = arealight

	return gp
}

func (p *GeometricPrimitive) WorldBound() *BBox {
	if p.shape != nil {
		return p.shape.WorldBound()
	}
	return nil
}

func (p *GeometricPrimitive) CanIntersect() bool {
	if p.shape != nil {
		return p.shape.CanIntersect()
	}
	return false
}

func (p *GeometricPrimitive) Intersect(r RayBase) (hit bool, isect *Intersection) {
	var thit, rayEpsilon float64
	var dg *DifferentialGeometry

	if hit, thit, rayEpsilon, dg = p.shape.Intersect(r); !hit {
		return false, nil
	}

	isect = new(Intersection)
	isect.dg = dg
	isect.primitive = p
	isect.WorldToObject = p.shape.WorldToObject()
	isect.ObjectToWorld = p.shape.ObjectToWorld()
	isect.shapeId = p.shape.ShapeId()
	isect.primitiveId = p.primitiveId
	isect.rayEpsilon = rayEpsilon
	r.SetMaxt(thit)

	return true, isect
}

func (p *GeometricPrimitive) IntersectP(r RayBase) bool {
	return p.shape.IntersectP(r)
}

func (p *GeometricPrimitive) Refine(refined []Primitive) []Primitive {
	r := p.shape.Refine(nil)
	for _, prim := range r {
		gp := NewGeometricPrimitive(prim, p.material, p.areaLight)
		if refined == nil {
			refined = make([]Primitive, 0, 2)
		}
		refined = append(refined, gp)
	}
	return refined
}

func (p *GeometricPrimitive) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(p, refined)
}

func (p *GeometricPrimitive) GetAreaLight() AreaLight {
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
	primitive        Primitive
	worldToPrimitive *AnimatedTransform
}

func NewTransformedPrimitive(prim Primitive, worldToPrimitive *AnimatedTransform) *TransformedPrimitive {
	tp := new(TransformedPrimitive)
	tp.primitiveId = GeneratePrimitiveId()
	tp.primitive = prim
	tp.worldToPrimitive = worldToPrimitive

	return tp
}

func (p *TransformedPrimitive) WorldBound() *BBox {
	return MotionBoundsAnimatedTransform(p.worldToPrimitive, p.primitive.WorldBound(), true)
}

func (p *TransformedPrimitive) CanIntersect() bool {
	return true
}

func (p *TransformedPrimitive) Intersect(r RayBase) (hit bool, isect *Intersection) {
	w2p := p.worldToPrimitive.Interpolate(r.Time())
	ray := r.Transform(w2p)
	if hit, isect = p.primitive.Intersect(ray); !hit {
		return false, nil
	}

	r.SetMaxt(ray.Maxt())
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

func (p *TransformedPrimitive) IntersectP(r RayBase) bool {
	return p.primitive.IntersectP(r.AnimatedTransform(p.worldToPrimitive))
}

func (p *TransformedPrimitive) Refine(refined []Primitive) []Primitive {
	return refined
}

func (p *TransformedPrimitive) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(p, refined)
}

func (p *TransformedPrimitive) GetAreaLight() AreaLight {
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

func (p *Aggregate) Intersect(r RayBase) (bool, *Intersection) {
	return false, nil
}

func (p *Aggregate) IntersectP(r RayBase) bool {
	return true
}

func (p *Aggregate) Refine(refined []Primitive) []Primitive {
	return refined
}

func (p *Aggregate) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(p, refined)
}

func (p *Aggregate) GetAreaLight() AreaLight {
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

type NoneAccel struct {
	PrimitiveData
	primitives []Primitive
	bounds     BBox
}

func NewNoneAccelerator(prims []Primitive) *NoneAccel {
	none := new(NoneAccel)
	none.primitiveId = GeneratePrimitiveId()
	none.primitives = make([]Primitive, 0, 16)

	// Initialize _primitives_ with primitives for none
	for _, p := range prims {
		none.primitives = p.FullyRefine(none.primitives)
	}

	// Compute bounds and choose grid resolution
	none.bounds = *CreateEmptyBBox()
	for _, p := range none.primitives {
		none.bounds = *UnionBBoxes(&none.bounds, p.WorldBound())
	}
	return none

}
func (g *NoneAccel) WorldBound() *BBox { return &g.bounds }

func (g *NoneAccel) CanIntersect() bool { return true }

func (g *NoneAccel) Intersect(ray RayBase) (bool, *Intersection) {
	var hit bool
	var isect *Intersection

	for _, p := range g.primitives {
		if phit, pisect := p.Intersect(ray); phit {
			hit = phit
			isect = pisect
		}
	}
	return hit, isect
}

func (g *NoneAccel) IntersectP(ray RayBase) bool {
	for _, p := range g.primitives {
		if p.IntersectP(ray) {
			return true
		}
	}
	return false
}

func (g *NoneAccel) Refine(refined []Primitive) []Primitive { return refined }
func (g *NoneAccel) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(g, refined)
}

func (g *NoneAccel) GetAreaLight() AreaLight { return nil }
func (g *NoneAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (g *NoneAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (g *NoneAccel) PrimitiveId() uint32 { return g.primitiveId }

func CreateNoneAccelerator(prims []Primitive, ps *ParamSet) *NoneAccel {
	return NewNoneAccelerator(prims)
}
