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
)

type (
	Voxel struct {
		primitives      []Primitive
		allCanIntersect bool
	}

	GridAccel struct {
		PrimitiveData
		primitives      []Primitive
		nVoxels         [3]int
		bounds          BBox
		width, invWidth Vector
		voxels          []*Voxel
		voxelArena      *MemoryArena
	}
)

func NewGridAccel(prims []Primitive, refineImmediately bool) *GridAccel {
	grid := new(GridAccel)
	grid.primitiveId = GeneratePrimitiveId()
	grid.primitives = make([]Primitive, 0, 16)

	// Initialize _primitives_ with primitives for grid
	if refineImmediately {
		for _, p := range prims {
			grid.primitives = p.FullyRefine(grid.primitives)
		}
	} else {
		grid.primitives = append(grid.primitives, prims...)
	}

	// Compute bounds and choose grid resolution
	grid.bounds = *CreateEmptyBBox()
	for _, p := range grid.primitives {
		grid.bounds = *UnionBBoxes(&grid.bounds, p.WorldBound())
	}
	delta := grid.bounds.PMax.Sub(&grid.bounds.PMin)

	// Find _voxelsPerUnitDist_ for grid
	maxAxis := grid.bounds.MaximumExtent()
	invMaxWidth := 1.0 / delta.At(maxAxis)
	Assert(invMaxWidth > 0.0)
	cubeRoot := 3.0 * math.Pow(float64(len(grid.primitives)), 1.0/3.0)
	voxelsPerUnitDist := cubeRoot * invMaxWidth
	for axis := 0; axis < 3; axis++ {
		grid.nVoxels[axis] = Round2Int(delta.At(axis) * voxelsPerUnitDist)
		grid.nVoxels[axis] = Clampi(grid.nVoxels[axis], 1, 64)
	}

	// Compute voxel widths and allocate voxels
	grid.width.X = delta.X / float64(grid.nVoxels[0])
	if grid.width.X == 0 {
		grid.invWidth.X = 0.0
	} else {
		grid.invWidth.X = 1.0 / grid.width.X
	}
	grid.width.Y = delta.Y / float64(grid.nVoxels[1])
	if grid.width.Y == 0 {
		grid.invWidth.Y = 0.0
	} else {
		grid.invWidth.Y = 1.0 / grid.width.Y
	}
	grid.width.Z = delta.Z / float64(grid.nVoxels[2])
	if grid.width.Z == 0 {
		grid.invWidth.Z = 0.0
	} else {
		grid.invWidth.Z = 1.0 / grid.width.Z
	}

	nv := grid.nVoxels[0] * grid.nVoxels[1] * grid.nVoxels[2]
	grid.voxels = make([]*Voxel, nv, nv)

	// Add primitives to grid voxels
	for _, p := range grid.primitives {
		// Find voxel extent of primitive
		pb := p.WorldBound()
		var vmin, vmax [3]int
		for axis := 0; axis < 3; axis++ {
			vmin[axis] = grid.posToVoxel(&pb.PMin, axis)
			vmax[axis] = grid.posToVoxel(&pb.PMax, axis)
		}

		// Add primitive to overlapping voxels
		for z := vmin[2]; z <= vmax[2]; z++ {
			for y := vmin[1]; y <= vmax[1]; y++ {
				for x := vmin[0]; x <= vmax[0]; x++ {
					o := grid.offset(x, y, z)
					if grid.voxels[o] == nil {
						// Allocate new voxel and store primitive in it
						grid.voxels[o] = newVoxel(p)
					} else {
						// Add primitive to already-allocated voxel
						grid.voxels[o].addPrimitive(p)
					}
				}
			}
		}
	}

	return grid
}

func (g *GridAccel) WorldBound() *BBox { return &g.bounds }

func (g *GridAccel) CanIntersect() bool { return true }

func (g *GridAccel) Intersect(ray RayBase) (bool, *Intersection) {
	// Check ray against overall grid bounds
	var rayT float64
	var ok bool
	if g.bounds.Inside(ray.PointAt(ray.Mint())) {
		rayT = ray.Mint()
	} else if ok, rayT, _ = g.bounds.IntersectP(ray); !ok {
		return false, nil
	}
	gridIntersect := ray.PointAt(rayT)

	// Set up 3D DDA for ray
	var NextCrossingT, DeltaT [3]float64
	var Step, Out, Pos [3]int
	for axis := 0; axis < 3; axis++ {
		// Compute current voxel for axis
		Pos[axis] = g.posToVoxel(gridIntersect, axis)
		if ray.Dir().At(axis) >= 0 {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(g.voxelToPos(Pos[axis]+1, axis)-gridIntersect.At(axis))/ray.Dir().At(axis)
			DeltaT[axis] = g.width.At(axis) / ray.Dir().At(axis)
			Step[axis] = 1
			Out[axis] = g.nVoxels[axis]
		} else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(g.voxelToPos(Pos[axis], axis)-gridIntersect.At(axis))/ray.Dir().At(axis)
			DeltaT[axis] = -g.width.At(axis) / ray.Dir().At(axis)
			Step[axis] = -1
			Out[axis] = -1
		}
	}

	// Walk ray through voxel grid
	hitSomething := false
	var isect *Intersection
	for {
		// Check for intersection in current voxel and advance to next
		voxel := g.voxels[g.offset(Pos[0], Pos[1], Pos[2])]
		if voxel != nil {
			//var hit bool
			if hit, visect := voxel.intersect(ray); hit {
				hitSomething = true
				isect = visect
			}
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		bits := 0
		if NextCrossingT[0] < NextCrossingT[1] {
			bits += 1 << 2
		}
		if NextCrossingT[0] < NextCrossingT[2] {
			bits += 1 << 1
		}
		if NextCrossingT[1] < NextCrossingT[2] {
			bits += 1 << 0
		}

		cmpToAxis := [8]int{2, 1, 2, 1, 2, 2, 0, 0}
		stepAxis := cmpToAxis[bits]
		if ray.Maxt() < NextCrossingT[stepAxis] {
			break
		}
		Pos[stepAxis] += Step[stepAxis]
		if Pos[stepAxis] == Out[stepAxis] {
			break
		}
		NextCrossingT[stepAxis] += DeltaT[stepAxis]
	}
	return hitSomething, isect
}

func (g *GridAccel) IntersectP(ray RayBase) bool {
	// Check ray against overall grid bounds
	var rayT float64
	var ok bool
	if g.bounds.Inside(ray.PointAt(ray.Mint())) {
		rayT = ray.Mint()
	} else if ok, rayT, _ = g.bounds.IntersectP(ray); !ok {
		return false
	}
	gridIntersect := ray.PointAt(rayT)

	// Set up 3D DDA for ray
	var NextCrossingT, DeltaT [3]float64
	var Step, Out, Pos [3]int
	for axis := 0; axis < 3; axis++ {
		// Compute current voxel for axis
		Pos[axis] = g.posToVoxel(gridIntersect, axis)
		if ray.Dir().At(axis) >= 0 {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(g.voxelToPos(Pos[axis]+1, axis)-gridIntersect.At(axis))/ray.Dir().At(axis)
			DeltaT[axis] = g.width.At(axis) / ray.Dir().At(axis)
			Step[axis] = 1
			Out[axis] = g.nVoxels[axis]
		} else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(g.voxelToPos(Pos[axis], axis)-gridIntersect.At(axis))/ray.Dir().At(axis)
			DeltaT[axis] = -g.width.At(axis) / ray.Dir().At(axis)
			Step[axis] = -1
			Out[axis] = -1
		}
	}

	// Walk grid for shadow ray
	for {
		o := g.offset(Pos[0], Pos[1], Pos[2])
		voxel := g.voxels[o]
		if voxel != nil && voxel.intersectP(ray) {
			return true
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		bits := 0
		if NextCrossingT[0] < NextCrossingT[1] {
			bits += 1 << 2
		}
		if NextCrossingT[0] < NextCrossingT[2] {
			bits += 1 << 1
		}
		if NextCrossingT[1] < NextCrossingT[2] {
			bits += 1 << 0
		}
		cmpToAxis := [8]int{2, 1, 2, 1, 2, 2, 0, 0}
		stepAxis := cmpToAxis[bits]
		if ray.Maxt() < NextCrossingT[stepAxis] {
			break
		}
		Pos[stepAxis] += Step[stepAxis]
		if Pos[stepAxis] == Out[stepAxis] {
			break
		}
		NextCrossingT[stepAxis] += DeltaT[stepAxis]
	}
	return false
}

func (g *GridAccel) Refine(refined []Primitive) []Primitive { return refined }
func (g *GridAccel) FullyRefine(refined []Primitive) []Primitive {
	return PrimitiveFullyRefine(g, refined)
}

func (g *GridAccel) GetAreaLight() AreaLight { return nil }
func (g *GridAccel) GetBSDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSDF {
	return nil
}
func (g *GridAccel) GetBSSRDF(dg *DifferentialGeometry, objectToWorld *Transform, arena *MemoryArena) *BSSRDF {
	return nil
}
func (g *GridAccel) PrimitiveId() uint32 { return g.primitiveId }

func (g *GridAccel) posToVoxel(P *Point, axis int) int {
	v := Float2Int((P.At(axis) - g.bounds.PMin.At(axis)) * g.invWidth.At(axis))
	return Clampi(v, 0, g.nVoxels[axis]-1)
}
func (g *GridAccel) voxelToPos(p, axis int) float64 {
	return g.bounds.PMin.At(axis) + float64(p)*g.width.At(axis)
}
func (g *GridAccel) offset(x, y, z int) int {
	return z*g.nVoxels[0]*g.nVoxels[1] + y*g.nVoxels[0] + x
}

func CreateGridAccelerator(prims []Primitive, ps *ParamSet) *GridAccel {
	refineImmediately := ps.FindBoolParam("refineimmediately", false)
	return NewGridAccel(prims, refineImmediately)
}

func newVoxel(op Primitive) *Voxel {
	v := new(Voxel)
	v.allCanIntersect = false
	v.primitives = make([]Primitive, 1, 8)
	v.primitives[0] = op
	return v
}

func (v *Voxel) addPrimitive(prim Primitive) {
	v.primitives = append(v.primitives, prim)
}

func (v *Voxel) intersect(ray RayBase) (bool, *Intersection) {
	// Refine primitives in voxel if needed
	if !v.allCanIntersect {
		for i, prim := range v.primitives {
			// Refine primitive _prim_ if it's not intersectable
			if !prim.CanIntersect() {
				p := make([]Primitive, 0, 8)
				p = prim.FullyRefine(p)
				Assert(len(p) > 0)
				if len(p) == 1 {
					v.primitives[i] = p[0]
				} else {
					v.primitives[i] = NewGridAccel(p, false)
				}
			}
		}
		v.allCanIntersect = true
	}

	// Loop over primitives in voxel and find intersections
	hitSomething := false
	var isect *Intersection
	for _, prim := range v.primitives {
		if hit, visect := prim.Intersect(ray); hit {
			hitSomething = true
			isect = visect
		}
	}
	return hitSomething, isect
}

func (v *Voxel) intersectP(ray RayBase) bool {
	// Refine primitives in voxel if needed
	if !v.allCanIntersect {
		for i, prim := range v.primitives {
			// Refine primitive _prim_ if it's not intersectable
			if !prim.CanIntersect() {
				p := make([]Primitive, 0, 8)
				p = prim.FullyRefine(p)
				Assert(len(p) > 0)
				if len(p) == 1 {
					v.primitives[i] = p[0]
				} else {
					v.primitives[i] = NewGridAccel(p, false)
				}
			}
		}
		v.allCanIntersect = true
	}
	for _, prim := range v.primitives {
		if prim.IntersectP(ray) {
			return true
		}
	}
	return false
}
