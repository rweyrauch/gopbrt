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
	"math"
)

type (
	Vector struct {
		X, Y, Z float64
	}

	Point struct {
		X, Y, Z float64
	}

	Normal struct {
		X, Y, Z float64
	}

	Ray struct {
		origin     Point
		dir        Vector
		mint, maxt float64
		time       float64
		depth      int
	}

	RayDifferential struct {
		Ray
		hasDifferentials         bool
		rxOrigin, ryOrigin       Point
		rxDirection, ryDirection Vector
	}

	BBox struct {
		pMin, pMax Point
	}
)

const (
	INFINITY  = math.MaxFloat64
	INFINITYF = math.MaxFloat32
)

func CreateVector(x, y, z float64) *Vector {
	return &Vector{x, y, z}
}
func CreateVectorFromNormal(n *Normal) *Vector {
	return &Vector{n.X, n.Y, n.Z}
}
func CreateVectorFromPoint(p *Point) *Vector {
	return &Vector{p.X, p.Y, p.Z}
}
func (v *Vector) Negate() *Vector {
	return &Vector{-v.X, -v.Y, -v.Z}
}
func (v1 *Vector) Add(v2 *Vector) *Vector {
	return &Vector{v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z}
}
func (v1 *Vector) Sub(v2 *Vector) *Vector {
	return &Vector{v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z}
}
func (v *Vector) Scale(s float64) *Vector {
	return &Vector{s * v.X, s * v.Y, s * v.Z}
}
func (v *Vector) InvScale(invs float64) *Vector {
	s := 1.0 / invs
	return &Vector{s * v.X, s * v.Y, s * v.Z}
}
func (v *Vector) LengthSquared() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}
func (v *Vector) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}
func (v *Vector) At(i int) float64 {
	if i == 0 {
		return v.X
	} else if i == 1 {
		return v.Y
	} else if i == 2 {
		return v.Z
	}
	panic(fmt.Errorf("Vector.At: Invalid index %d", i))
}
func (v *Vector) Set(i int, val float64) {
	if i == 0 {
		v.X = val
	} else if i == 1 {
		v.Y = val
	} else if i == 2 {
		v.Z = val
	} else {
		panic(fmt.Errorf("Vector.Set: Invalid index %d", i))
	}
}
func EqualVector(v1, v2 *Vector) bool {
	return v1.X == v2.X && v1.Y == v2.Y && v1.Z == v2.Z
}
func NotEqualVector(v1, v2 *Vector) bool {
	return v1.X != v2.X || v1.Y != v2.Y || v1.Z != v2.Z
}
func DotVector(v1, v2 *Vector) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}
func AbsDotVector(v1, v2 *Vector) float64 {
	return math.Abs(DotVector(v1, v2))
}
func CrossVector(v1, v2 *Vector) *Vector {
	return &Vector{(v1.Y * v2.Z) - (v1.Z * v2.Y),
		(v1.Z * v2.X) - (v1.X * v2.Z),
		(v1.X * v2.Y) - (v1.Y * v2.X)}
}
func CrossVectorNormal(v1 *Vector, v2 *Normal) *Vector {
	return &Vector{(v1.Y * v2.Z) - (v1.Z * v2.Y),
		(v1.Z * v2.X) - (v1.X * v2.Z),
		(v1.X * v2.Y) - (v1.Y * v2.X)}
}
func CrossNormalVector(v1 *Normal, v2 *Vector) *Vector {
	return &Vector{(v1.Y * v2.Z) - (v1.Z * v2.Y),
		(v1.Z * v2.X) - (v1.X * v2.Z),
		(v1.X * v2.Y) - (v1.Y * v2.X)}
}
func NormalizeVector(v *Vector) *Vector {
	len := v.Length()
	return &Vector{v.X / len, v.Y / len, v.Z / len}
}
func CoordinateSystem(v1 *Vector) (v2, v3 *Vector) {
	if math.Abs(v1.X) > math.Abs(v1.Y) {
		invLen := 1.0 / math.Sqrt(v1.X*v1.X+v1.Z*v1.Z)
		v2 = &Vector{-v1.Z * invLen, 0.0, v1.X * invLen}
	} else {
		invLen := 1.0 / math.Sqrt(v1.Y*v1.Y+v1.Z*v1.Z)
		v2 = &Vector{0.0, v1.Z * invLen, -v1.Y * invLen}
	}
	v3 = CrossVector(v1, v2)
	return v2, v3
}
func (v *Vector) HasNaNs() bool {
	return math.IsNaN(v.X) || math.IsNaN(v.Y) || math.IsNaN(v.Z)
}

func (v *Vector) String() string {
	return fmt.Sprintf("v[ %f, %f, %f ]", v.X, v.Y, v.Z)
}

func SphericalDirection(sintheta, costheta, phi float64) *Vector {
	return &Vector{sintheta * math.Cos(phi), sintheta * math.Sin(phi), costheta}
}
func SphericalDirectionVectors(sintheta, costheta, phi float64, x, y, z *Vector) *Vector {
	return x.Scale(sintheta * math.Cos(phi)).Add(y.Scale(sintheta * math.Sin(phi)).Add(z.Scale(costheta)))
}
func SphericalTheta(v *Vector) float64 {
	return math.Acos(Clamp(v.Z, -1.0, 1.0))
}
func SphericalPhi(v *Vector) float64 {
	p := math.Atan2(v.Y, v.X)
	if p < 0.0 {
		return p + 2.0*math.Pi
	}
	return p
}

func CreatePoint(x, y, z float64) *Point {
	return &Point{x, y, z}
}
func (p1 *Point) Add(v2 *Vector) *Point {
	return &Point{p1.X + v2.X, p1.Y + v2.Y, p1.Z + v2.Z}
}
func (p1 *Point) AddPoint(p2 *Point) *Point {
	return &Point{p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z}
}
func (p1 *Point) Sub(p2 *Point) *Vector {
	return &Vector{p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z}
}
func (p1 *Point) SubVector(v2 *Vector) *Point {
	return &Point{p1.X - v2.X, p1.Y - v2.Y, p1.Z - v2.Z}
}
func (p *Point) Scale(s float64) *Point {
	return &Point{s * p.X, s * p.Y, s * p.Z}
}
func (p *Point) InvScale(invs float64) *Point {
	s := 1.0 / invs
	return &Point{s * p.X, s * p.Y, s * p.Z}
}
func (p *Point) At(i int) float64 {
	if i == 0 {
		return p.X
	} else if i == 1 {
		return p.Y
	} else if i == 2 {
		return p.Z
	}
	panic(fmt.Errorf("Point.At: Invalid index %d", i))
}
func (p *Point) Set(i int, val float64) {
	if i == 0 {
		p.X = val
	} else if i == 1 {
		p.Y = val
	} else if i == 2 {
		p.Z = val
	} else {
		panic(fmt.Errorf("Point.Set: Invalid index %d", i))
	}
}

func EqualPoint(p1, p2 *Point) bool {
	return p1.X == p2.X && p1.Y == p2.Y && p1.Z == p2.Z
}
func NotEqualPoint(p1, p2 *Point) bool {
	return p1.X != p2.X || p1.Y != p2.Y || p1.Z != p2.Z
}
func DistancePoint(p1, p2 *Point) float64 {
	return p1.Sub(p2).Length()
}
func DistanceSquaredPoint(p1, p2 *Point) float64 {
	return p1.Sub(p2).LengthSquared()
}
func (p *Point) HasNaNs() bool {
	return math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsNaN(p.Z)
}

func (p *Point) String() string {
	return fmt.Sprintf("p[ %f, %f, %f ]", p.X, p.Y, p.Z)
}

func CreateNormal(x, y, z float64) *Normal {
	return &Normal{x, y, z}
}
func CreateNormalFromVector(v *Vector) *Normal {
	return &Normal{v.X, v.Y, v.Z}
}
func (n *Normal) Negate() *Normal {
	return &Normal{-n.X, -n.Y, -n.Z}
}
func (n1 *Normal) Add(n2 *Normal) *Normal {
	return &Normal{n1.X + n2.X, n1.Y + n2.Y, n1.Z + n2.Z}
}
func (n1 *Normal) Sub(n2 *Normal) *Normal {
	return &Normal{n1.X - n2.X, n1.Y - n2.Y, n1.Z - n2.Z}
}
func (n *Normal) Scale(s float64) *Normal {
	return &Normal{s * n.X, s * n.Y, s * n.Z}
}
func (n *Normal) InvScale(invs float64) *Normal {
	s := 1.0 / invs
	return &Normal{s * n.X, s * n.Y, s * n.Z}
}
func (n *Normal) LengthSquared() float64 {
	return n.X*n.X + n.Y*n.Y + n.Z*n.Z
}
func (n *Normal) Length() float64 {
	return math.Sqrt(n.LengthSquared())
}
func (n *Normal) At(i int) float64 {
	if i == 0 {
		return n.X
	} else if i == 1 {
		return n.Y
	} else if i == 2 {
		return n.Z
	}
	panic(fmt.Errorf("Normal.At: Invalid index %d", i))
}

func EqualNormal(n1, n2 *Normal) bool {
	return n1.X == n2.X && n1.Y == n2.Y && n1.Z == n2.Z
}
func NotEqualNormal(n1, n2 *Normal) bool {
	return n1.X != n2.X || n1.Y != n2.Y || n1.Z != n2.Z
}
func NormalizeNormal(n *Normal) *Normal {
	len := n.Length()
	return &Normal{n.X / len, n.Y / len, n.Z / len}
}
func DotNormalVector(v1 *Normal, v2 *Vector) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}
func DotVectorNormal(v1 *Vector, v2 *Normal) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}
func DotNormal(v1, v2 *Normal) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}
func AbsDotNormalVector(v1 *Normal, v2 *Vector) float64 {
	return math.Abs(DotNormalVector(v1, v2))
}
func AbsDotVectorNormal(v1 *Vector, v2 *Normal) float64 {
	return math.Abs(DotVectorNormal(v1, v2))
}
func AbsDotNormal(v1, v2 *Normal) float64 {
	return math.Abs(DotNormal(v1, v2))
}
func FaceforwardNormalVector(n *Normal, v *Vector) *Normal {
	if DotNormalVector(n, v) < 0.0 {
		return n.Negate()
	}
	return n
}
func FaceforwardNormal(n, n2 *Normal) *Normal {
	if DotNormal(n, n2) < 0.0 {
		return n.Negate()
	}
	return n
}
func FaceforwardVector(v, v2 *Vector) *Vector {
	if DotVector(v, v2) < 0.0 {
		return v.Negate()
	}
	return v
}
func FaceforwardVectorNormal(v *Vector, n2 *Normal) *Vector {
	if DotVectorNormal(v, n2) < 0.0 {
		return v.Negate()
	}
	return v
}
func CrossNormal(v1 *Normal, v2 *Normal) *Vector {
	return &Vector{(v1.Y * v2.Z) - (v1.Z * v2.Y),
		(v1.Z * v2.X) - (v1.X * v2.Z),
		(v1.X * v2.Y) - (v1.Y * v2.X)}
}

func (n *Normal) HasNaNs() bool {
	return math.IsNaN(n.X) || math.IsNaN(n.Y) || math.IsNaN(n.Z)
}

func (n *Normal) String() string {
	return fmt.Sprintf("n[ %f, %f, %f ]", n.X, n.Y, n.Z)
}

func CreateRay(origin *Point, direction *Vector, start, end, t float64, d int) *Ray {
	return &Ray{*origin, *direction, start, end, t, d}
}
func CreateChildRay(origin *Point, direction *Vector, parent *Ray, start, end float64) *Ray {
	return &Ray{*origin, *direction, start, end, parent.time, parent.depth + 1}
}

func (r *Ray) PointAt(t float64) *Point {
	return r.origin.Add(r.dir.Scale(t))
}

func (ray *Ray) HasNaNs() bool {
	return ray.dir.HasNaNs() || ray.origin.HasNaNs() || math.IsNaN(ray.mint) || math.IsNaN(ray.maxt)
}

func (r *Ray) String() string {
	return fmt.Sprintf("ray[origin:%v dir:%v start:%f end:%f time:%f depth:%d]", &r.origin, &r.dir, r.mint, r.maxt, r.time, r.depth)
}
func (r *RayDifferential) String() string {
	return fmt.Sprintf("raydiff[origin:%v dir:%v start:%f end:%f time:%f depth:%d diffs:%s]", &r.origin, &r.dir, r.mint, r.maxt, r.time, r.depth, r.hasDifferentials)
}

func CreateRayDifferential(origin *Point, direction *Vector, start, end, t float64, d int) *RayDifferential {
	return &RayDifferential{Ray{*origin, *direction, start, end, t, d}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func CreateChildRayDifferential(origin *Point, direction *Vector, parent *Ray, start, end float64) *RayDifferential {
	return &RayDifferential{Ray{*origin, *direction, start, end, parent.time, parent.depth + 1}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func CreateRayDifferentialFromRay(ray *Ray) *RayDifferential {
	return &RayDifferential{Ray{ray.origin, ray.dir, ray.mint, ray.maxt, ray.time, ray.depth}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func CreateRayFromRayDifferential(ray *RayDifferential) *Ray {
	return &Ray{ray.origin, ray.dir, ray.mint, ray.maxt, ray.time, ray.depth}
}

func (r *RayDifferential) ScaleDifferentials(s float64) {
	r.rxOrigin = *r.origin.Add(r.rxOrigin.Sub(&r.origin).Scale(s))
	r.ryOrigin = *r.origin.Add(r.ryOrigin.Sub(&r.origin).Scale(s))
	r.rxDirection = *r.dir.Add(r.rxDirection.Sub(&r.dir).Scale(s))
	r.ryDirection = *r.dir.Add(r.ryDirection.Sub(&r.dir).Scale(s))
}

func (ray *RayDifferential) HasNaNs() bool {
	return ray.dir.HasNaNs() || ray.origin.HasNaNs() || math.IsNaN(ray.mint) || math.IsNaN(ray.maxt) ||
		(ray.hasDifferentials && ray.rxOrigin.HasNaNs() || ray.ryOrigin.HasNaNs() || ray.rxDirection.HasNaNs() || ray.ryDirection.HasNaNs())
}

func CreateEmptyBBox() *BBox {
	return &BBox{Point{INFINITY, INFINITY, INFINITY}, Point{-INFINITY, -INFINITY, -INFINITY}}
}
func CreateBBoxFromPoint(p *Point) *BBox {
	return &BBox{*p, *p}
}
func CreateBBoxFromPoints(p1, p2 *Point) *BBox {
	bbox := new(BBox)
	bbox.pMin.X, bbox.pMin.Y, bbox.pMin.Z = math.Min(p1.X, p2.X), math.Min(p1.Y, p2.Y), math.Min(p1.Z, p2.Z)
	bbox.pMax.X, bbox.pMax.Y, bbox.pMax.Z = math.Max(p1.X, p2.X), math.Max(p1.Y, p2.Y), math.Max(p1.Z, p2.Z)
	return bbox
}
func UnionBBoxPoint(b *BBox, p *Point) *BBox {
	return &BBox{Point{math.Min(b.pMin.X, p.X), math.Min(b.pMin.Y, p.Y), math.Min(b.pMin.Z, p.Z)},
		Point{math.Max(b.pMax.X, p.X), math.Max(b.pMax.Y, p.Y), math.Max(b.pMax.Z, p.Z)}}
}

func UnionBBoxes(b1, b2 *BBox) *BBox {
	return &BBox{Point{math.Min(b1.pMin.X, b2.pMin.X), math.Min(b1.pMin.Y, b2.pMin.Y), math.Min(b1.pMin.Z, b2.pMin.Z)},
		Point{math.Max(b1.pMax.X, b2.pMax.X), math.Max(b1.pMax.Y, b2.pMax.Y), math.Max(b1.pMax.Z, b2.pMax.Z)}}
}
func (b1 *BBox) Overlaps(b2 *BBox) bool {
	x := (b1.pMax.X >= b2.pMin.X) && (b1.pMin.X <= b2.pMax.X)
	y := (b1.pMax.Y >= b2.pMin.Y) && (b1.pMin.Y <= b2.pMax.Y)
	z := (b1.pMax.Z >= b2.pMin.Z) && (b1.pMin.Z <= b2.pMax.Z)
	return x && y && z
}
func (b *BBox) Inside(pt *Point) bool {
	return (pt.X >= b.pMin.X && pt.X <= b.pMax.X &&
		pt.Y >= b.pMin.Y && pt.Y <= b.pMax.Y &&
		pt.Z >= b.pMin.Z && pt.Z <= b.pMax.Z)
}
func (b *BBox) Expand(delta float64) {
	b.pMin = Point{b.pMin.X - delta, b.pMin.Y - delta, b.pMin.Z - delta}
	b.pMax = Point{b.pMax.X + delta, b.pMax.Y + delta, b.pMax.Z + delta}
}
func (b *BBox) SurfaceArea() float64 {
	d := b.pMax.Sub(&b.pMin)
	return 2.0 * (d.X*d.Y + d.X*d.Z + d.Y*d.Z)
}
func (b *BBox) Volume() float64 {
	d := b.pMax.Sub(&b.pMin)
	return d.X * d.Y * d.Z
}
func (b *BBox) MaximumExtent() int {
	diag := b.pMax.Sub(&b.pMin)
	if diag.X > diag.Y && diag.X > diag.Z {
		return 0
	} else if diag.Y > diag.Z {
		return 1
	} else {
		return 2
	}
}
func (b *BBox) PointAtIndex(i int) *Point {
	if i == 0 {
		return &b.pMin
	} else if i == 1 {
		return &b.pMax
	}
	panic(fmt.Errorf("PointAtIndex: Invalid index, %d", i))
}

func (b *BBox) Lerp(tx, ty, tz float64) *Point {
	return &Point{Lerp(tx, b.pMin.X, b.pMax.X), Lerp(ty, b.pMin.Y, b.pMax.Y), Lerp(tz, b.pMin.Z, b.pMax.Z)}
}
func (b *BBox) Offset(p *Point) *Vector {
	return &Vector{(p.X - b.pMin.X) / (b.pMax.X - b.pMin.X),
		(p.Y - b.pMin.Y) / (b.pMax.Y - b.pMin.Y),
		(p.Z - b.pMin.Z) / (b.pMax.Z - b.pMin.Z)}
}
func (b *BBox) BoundingSphere() (c *Point, rad float64) {
	c = b.pMin.Add(&Vector{b.pMax.X, b.pMax.Y, b.pMax.Z}).Scale(0.5)
	if b.Inside(c) {
		rad = DistancePoint(c, &b.pMax)
	} else {
		rad = 0.0
	}
	return c, rad
}
func (b *BBox) IntersectP(ray *Ray) (bool, float64, float64) {
	t0, t1 := ray.mint, ray.maxt
	for i := 0; i < 3; i++ {
		// Update interval for _i_th bounding box slab
		invRayDir := 1.0 / ray.dir.At(i)
		tNear := (b.pMin.At(i) - ray.origin.At(i)) * invRayDir
		tFar := (b.pMax.At(i) - ray.origin.At(i)) * invRayDir

		// Update parametric interval from slab intersection $t$s
		if tNear > tFar {
			tNear, tFar = tFar, tNear
		}
		if tNear > t0 {
			t0 = tNear
		}
		if tFar < t1 {
			t1 = tFar
		}
		if t0 > t1 {
			return false, 0.0, 0.0
		}
	}
	return true, t0, t1
}
func EqualBBox(b1, b2 *BBox) bool {
	return EqualPoint(&b1.pMin, &b2.pMin) && EqualPoint(&b1.pMax, &b2.pMax)
}
func NotEqualBBox(b1, b2 *BBox) bool {
	return NotEqualPoint(&b1.pMin, &b2.pMin) || NotEqualPoint(&b1.pMax, &b2.pMax)
}
func (b *BBox) String() string {
	return fmt.Sprintf("bbox[ min: %v max: %v ]", &b.pMin, &b.pMax)
}
