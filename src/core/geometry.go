package pbrt

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
		Origin          Point
		Dir          Vector
		Mint, Maxt float64
		Time       float64
		Depth      int
	}

	RayDifferential struct {
		Ray
		HasDifferentials         bool
		RxOrigin, RyOrigin       Point
		RxDirection, RyDirection Vector
	}

	BBox struct {
		PMin, PMax Point
	}
)

const (
	INFINITY = math.MaxFloat64
)

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

func (p1 *Point) Add(v2 *Vector) *Point {
	return &Point{p1.X + v2.X, p1.Y + v2.Y, p1.Z + v2.Z}
}
func (p1 *Point) Sub(p2 *Point) *Vector {
	return &Vector{p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z}
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

func CreateRay(origin *Point, direction *Vector, start, end, t float64, d int) *Ray {
	return &Ray{*origin, *direction, start, end, t, d}
}
func CreateChildRay(origin *Point, direction *Vector, parent *Ray, start, end float64) *Ray {
	return &Ray{*origin, *direction, start, end, parent.Time, parent.Depth + 1}
}
func (r *Ray) PointAt(t float64) *Point {
	return r.Origin.Add(r.Dir.Scale(t))
}

func CreateRayDifferential(origin *Point, direction *Vector, start, end, t float64, d int) *RayDifferential {
	return &RayDifferential{Ray{*origin, *direction, start, end, t, d}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func CreateChildRayDifferential(origin *Point, direction *Vector, parent *Ray, start, end float64) *RayDifferential {
	return &RayDifferential{Ray{*origin, *direction, start, end, parent.Time, parent.Depth + 1}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func CreateRayDifferentialFromRay(ray *Ray) *RayDifferential {
	return &RayDifferential{Ray{ray.Origin, ray.Dir, ray.Mint, ray.Maxt, ray.Time, ray.Depth}, false, Point{0, 0, 0}, Point{0, 0, 0}, Vector{0, 0, 0}, Vector{0, 0, 0}}
}
func (r *RayDifferential) ScaleDifferentials(s float64) {
	r.RxOrigin = *r.Origin.Add(r.RxOrigin.Sub(&r.Origin).Scale(s))
	r.RyOrigin = *r.Origin.Add(r.RyOrigin.Sub(&r.Origin).Scale(s))
	r.RxDirection = *r.Dir.Add(r.RxDirection.Sub(&r.Dir).Scale(s))
	r.RyDirection = *r.Dir.Add(r.RyDirection.Sub(&r.Dir).Scale(s))
}

func CreateEmptyBBox() *BBox {
	return &BBox{Point{INFINITY, INFINITY, INFINITY}, Point{-INFINITY, -INFINITY, -INFINITY}}
}
func CreateBBoxFromPoint(p *Point) *BBox {
	return &BBox{*p, *p}
}
func CreateBBoxFromPoints(p1, p2 *Point) *BBox {
	return &BBox{*p1, *p2}
}
func UnionBBoxPoint(b *BBox, p *Point) *BBox {
	return &BBox{Point{math.Min(b.PMin.X, p.X), math.Min(b.PMin.Y, p.Y), math.Min(b.PMin.Z, p.Z)},
		Point{math.Max(b.PMax.X, p.X), math.Max(b.PMax.Y, p.Y), math.Max(b.PMax.Z, p.Z)}}
}

func UnionBBoxes(b1, b2 *BBox) *BBox {
	return &BBox{Point{math.Min(b1.PMin.X, b2.PMin.X), math.Min(b1.PMin.Y, b2.PMin.Y), math.Min(b1.PMin.Z, b2.PMin.Z)},
		Point{math.Max(b1.PMax.X, b2.PMax.X), math.Max(b1.PMax.Y, b2.PMax.Y), math.Max(b1.PMax.Z, b2.PMax.Z)}}
}
func (b1 *BBox) Overlaps(b2 *BBox) bool {
	x := (b1.PMax.X >= b2.PMin.X) && (b1.PMin.X <= b2.PMax.X)
	y := (b1.PMax.Y >= b2.PMin.Y) && (b1.PMin.Y <= b2.PMax.Y)
	z := (b1.PMax.Z >= b2.PMin.Z) && (b1.PMin.Z <= b2.PMax.Z)
	return x && y && z
}
func (b *BBox) Inside(pt *Point) bool {
	return (pt.X >= b.PMin.X && pt.X <= b.PMax.X &&
		pt.Y >= b.PMin.Y && pt.Y <= b.PMax.Y &&
		pt.Z >= b.PMin.Z && pt.Z <= b.PMax.Z)
}
func (b *BBox) Expand(delta float64) {
	b.PMin = Point{b.PMin.X - delta, b.PMin.Y - delta, b.PMin.Z - delta}
	b.PMax = Point{b.PMax.X + delta, b.PMax.Y + delta, b.PMax.Z + delta}
}
func (b *BBox) SurfaceArea() float64 {
	d := b.PMax.Sub(&b.PMin)
	return 2.0 * (d.X*d.Y + d.X*d.Z + d.Y*d.Z)
}
func (b *BBox) Volume() float64 {
	d := b.PMax.Sub(&b.PMin)
	return d.X * d.Y * d.Z
}
func (b *BBox) MaximumExtent() int {
	diag := b.PMax.Sub(&b.PMin)
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
		return &b.PMin
	} else if i == 1 {
		return &b.PMax
	}
	panic(fmt.Errorf("PointAtIndex: Invalid index, %d", i))
}

func (b *BBox) Lerp(tx, ty, tz float64) *Point {
	return &Point{Lerp(tx, b.PMin.X, b.PMax.X), Lerp(ty, b.PMin.Y, b.PMax.Y), Lerp(tz, b.PMin.Z, b.PMax.Z)}
}
func (b *BBox) Offset(p *Point) *Vector {
	return &Vector{(p.X - b.PMin.X) / (b.PMax.X - b.PMin.X),
		(p.Y - b.PMin.Y) / (b.PMax.Y - b.PMin.Y),
		(p.Z - b.PMin.Z) / (b.PMax.Z - b.PMin.Z)}
}
func (b *BBox) BoundingSphere() (c *Point, rad float64) {
	c = b.PMin.Add(&Vector{b.PMax.X, b.PMax.Y, b.PMax.Z}).Scale(0.5)
	if b.Inside(c) {
		rad = DistancePoint(c, &b.PMax)
	} else {
		rad = 0.0
	}
	return c, rad
}
func (b *BBox) IntersectP(ray *Ray) (bool, float64, float64) {
	t0, t1 := ray.Mint, ray.Maxt
	for i := 0; i < 3; i++ {
		// Update interval for _i_th bounding box slab
		invRayDir := 1.0 / ray.Dir.At(i)
		tNear := (b.PMin.At(i) - ray.Origin.At(i)) * invRayDir
		tFar := (b.PMax.At(i) - ray.Origin.At(i)) * invRayDir

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
	return EqualPoint(&b1.PMin, &b2.PMin) && EqualPoint(&b1.PMax, &b2.PMax)
}
func NotEqualBBox(b1, b2 *BBox) bool {
	return NotEqualPoint(&b1.PMin, &b2.PMin) || NotEqualPoint(&b1.PMax, &b2.PMax)
}
