package core

import (
	"fmt"
	"math"
)

type (
	Vector struct {
		x, y, z float64
	}

	Point struct {
		x, y, z float64
	}

	Normal struct {
		x, y, z float64
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
	return &Vector{n.x, n.y, n.z}
}
func CreateVectorFromPoint(p *Point) *Vector {
	return &Vector{p.x, p.y, p.z}
}
func (v *Vector) Negate() *Vector {
	return &Vector{-v.x, -v.y, -v.z}
}
func (v1 *Vector) Add(v2 *Vector) *Vector {
	return &Vector{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}
}
func (v1 *Vector) Sub(v2 *Vector) *Vector {
	return &Vector{v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}
}
func (v *Vector) Scale(s float64) *Vector {
	return &Vector{s * v.x, s * v.y, s * v.z}
}
func (v *Vector) InvScale(invs float64) *Vector {
	s := 1.0 / invs
	return &Vector{s * v.x, s * v.y, s * v.z}
}
func (v *Vector) LengthSquared() float64 {
	return v.x*v.x + v.y*v.y + v.z*v.z
}
func (v *Vector) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}
func (v *Vector) At(i int) float64 {
	if i == 0 {
		return v.x
	} else if i == 1 {
		return v.y
	} else if i == 2 {
		return v.z
	}
	panic(fmt.Errorf("Vector.At: Invalid index %d", i))
}
func (v *Vector) Set(i int, val float64) {
	if i == 0 {
		v.x = val
	} else if i == 1 {
		v.y = val
	} else if i == 2 {
		v.z = val
	}
	panic(fmt.Errorf("Vector.Set: Invalid index %d", i))	
}
func EqualVector(v1, v2 *Vector) bool {
	return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z
}
func NotEqualVector(v1, v2 *Vector) bool {
	return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z
}
func DotVector(v1, v2 *Vector) float64 {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
}
func AbsDotVector(v1, v2 *Vector) float64 {
	return math.Abs(DotVector(v1, v2))
}
func CrossVector(v1, v2 *Vector) *Vector {
	return &Vector{(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)}
}
func CrossVectorNormal(v1 *Vector, v2 *Normal) *Vector {
	return &Vector{(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)}
}
func CrossNormalVector(v1 *Normal, v2 *Vector) *Vector {
	return &Vector{(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)}
}
func NormalizeVector(v *Vector) *Vector {
	len := v.Length()
	return &Vector{v.x / len, v.y / len, v.z / len}
}
func CoordinateSystem(v1 *Vector) (v2, v3 *Vector) {
	if math.Abs(v1.x) > math.Abs(v1.y) {
		invLen := 1.0 / math.Sqrt(v1.x*v1.x+v1.z*v1.z)
		v2 = &Vector{-v1.z * invLen, 0.0, v1.x * invLen}
	} else {
		invLen := 1.0 / math.Sqrt(v1.y*v1.y+v1.z*v1.z)
		v2 = &Vector{0.0, v1.z * invLen, -v1.y * invLen}
	}
	v3 = CrossVector(v1, v2)
	return v2, v3
}
func (v *Vector) HasNaNs() bool {
	return math.IsNaN(v.x) || math.IsNaN(v.y) || math.IsNaN(v.z)
}

func (v *Vector) String() string {
	return fmt.Sprintf("v[ %f, %f, %f ]", v.x, v.y, v.z)
}

func SphericalDirection(sintheta, costheta, phi float64) *Vector {
	return &Vector{sintheta * math.Cos(phi), sintheta * math.Sin(phi), costheta}
}
func SphericalDirectionVectors(sintheta, costheta, phi float64, x, y, z *Vector) *Vector {
	return x.Scale(sintheta * math.Cos(phi)).Add(y.Scale(sintheta * math.Sin(phi)).Add(z.Scale(costheta)))
}
func SphericalTheta(v *Vector) float64 {
	return math.Acos(Clamp(v.z, -1.0, 1.0))
}
func SphericalPhi(v *Vector) float64 {
	p := math.Atan2(v.y, v.x)
	if p < 0.0 {
		return p + 2.0*math.Pi
	}
	return p
}

func CreatePoint(x, y, z float64) *Point {
	return &Point{x, y, z}
}
func (p1 *Point) Add(v2 *Vector) *Point {
	return &Point{p1.x + v2.x, p1.y + v2.y, p1.z + v2.z}
}
func (p1 *Point) AddPoint(p2 *Point) *Point {
	return &Point{p1.x + p2.x, p1.y + p2.y, p1.z + p2.z}
}
func (p1 *Point) Sub(p2 *Point) *Vector {
	return &Vector{p1.x - p2.x, p1.y - p2.y, p1.z - p2.z}
}
func (p1 *Point) SubVector(v2 *Vector) *Point {
	return &Point{p1.x - v2.x, p1.y - v2.y, p1.z - v2.z}	
}
func (p *Point) Scale(s float64) *Point {
	return &Point{s * p.x, s * p.y, s * p.z}
}
func (p *Point) InvScale(invs float64) *Point {
	s := 1.0 / invs
	return &Point{s * p.x, s * p.y, s * p.z}
}
func (p *Point) At(i int) float64 {
	if i == 0 {
		return p.x
	} else if i == 1 {
		return p.y
	} else if i == 2 {
		return p.z
	}
	panic(fmt.Errorf("Point.At: Invalid index %d", i))
}
func (p *Point) Set(i int, val float64) {
	if i == 0 {
		p.x = val
	} else if i == 1 {
		p.y = val
	} else if i == 2 {
		p.z = val
	}
	panic(fmt.Errorf("Point.Set: Invalid index %d", i))	
}

func EqualPoint(p1, p2 *Point) bool {
	return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z
}
func NotEqualPoint(p1, p2 *Point) bool {
	return p1.x != p2.x || p1.y != p2.y || p1.z != p2.z
}
func DistancePoint(p1, p2 *Point) float64 {
	return p1.Sub(p2).Length()
}
func DistanceSquaredPoint(p1, p2 *Point) float64 {
	return p1.Sub(p2).LengthSquared()
}
func (p *Point) HasNaNs() bool {
	return math.IsNaN(p.x) || math.IsNaN(p.y) || math.IsNaN(p.z)
}

func (p *Point) String() string {
	return fmt.Sprintf("p[ %f, %f, %f ]", p.x, p.y, p.z)
}

func CreateNormal(x, y, z float64) *Normal {
	return &Normal{x, y, z}
}
func CreateNormalFromVector(v *Vector) *Normal {
	return &Normal{v.x, v.y, v.z}
}
func (n *Normal) Negate() *Normal {
	return &Normal{-n.x, -n.y, -n.z}
}
func (n1 *Normal) Add(n2 *Normal) *Normal {
	return &Normal{n1.x + n2.x, n1.y + n2.y, n1.z + n2.z}
}
func (n1 *Normal) Sub(n2 *Normal) *Normal {
	return &Normal{n1.x - n2.x, n1.y - n2.y, n1.z - n2.z}
}
func (n *Normal) Scale(s float64) *Normal {
	return &Normal{s * n.x, s * n.y, s * n.z}
}
func (n *Normal) InvScale(invs float64) *Normal {
	s := 1.0 / invs
	return &Normal{s * n.x, s * n.y, s * n.z}
}
func (n *Normal) LengthSquared() float64 {
	return n.x*n.x + n.y*n.y + n.z*n.z
}
func (n *Normal) Length() float64 {
	return math.Sqrt(n.LengthSquared())
}
func (n *Normal) At(i int) float64 {
	if i == 0 {
		return n.x
	} else if i == 1 {
		return n.y
	} else if i == 2 {
		return n.z
	}
	panic(fmt.Errorf("Normal.At: Invalid index %d", i))
}

func EqualNormal(n1, n2 *Normal) bool {
	return n1.x == n2.x && n1.y == n2.y && n1.z == n2.z
}
func NotEqualNormal(n1, n2 *Normal) bool {
	return n1.x != n2.x || n1.y != n2.y || n1.z != n2.z
}
func NormalizeNormal(n *Normal) *Normal {
	len := n.Length()
	return &Normal{n.x / len, n.y / len, n.z / len}
}
func DotNormalVector(v1 *Normal, v2 *Vector) float64 {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
}
func DotVectorNormal(v1 *Vector, v2 *Normal) float64 {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
}
func DotNormal(v1, v2 *Normal) float64 {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
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
	return &Vector{(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)}
}

func (n *Normal) HasNaNs() bool {
	return math.IsNaN(n.x) || math.IsNaN(n.y) || math.IsNaN(n.z)
}

func (n *Normal) String() string {
	return fmt.Sprintf("n[ %f, %f, %f ]", n.x, n.y, n.z)
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
	return fmt.Sprintf("ray[origin: %v dir: %v start: %f end: %f time: %f", r.origin, r.dir, r.mint, r.maxt, r.time)
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
	bbox.pMin.x, bbox.pMin.y, bbox.pMin.z = math.Min(p1.x, p2.x), math.Min(p1.y, p2.y), math.Min(p1.z, p2.z)
    bbox.pMax.x, bbox.pMax.y, bbox.pMax.z = math.Max(p1.x, p2.x), math.Max(p1.y, p2.y), math.Max(p1.z, p2.z)
	return bbox
}
func UnionBBoxPoint(b *BBox, p *Point) *BBox {
	return &BBox{Point{math.Min(b.pMin.x, p.x), math.Min(b.pMin.y, p.y), math.Min(b.pMin.z, p.z)},
		Point{math.Max(b.pMax.x, p.x), math.Max(b.pMax.y, p.y), math.Max(b.pMax.z, p.z)}}
}

func UnionBBoxes(b1, b2 *BBox) *BBox {
	return &BBox{Point{math.Min(b1.pMin.x, b2.pMin.x), math.Min(b1.pMin.y, b2.pMin.y), math.Min(b1.pMin.z, b2.pMin.z)},
		Point{math.Max(b1.pMax.x, b2.pMax.x), math.Max(b1.pMax.y, b2.pMax.y), math.Max(b1.pMax.z, b2.pMax.z)}}
}
func (b1 *BBox) Overlaps(b2 *BBox) bool {
	x := (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x)
	y := (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y)
	z := (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z)
	return x && y && z
}
func (b *BBox) Inside(pt *Point) bool {
	return (pt.x >= b.pMin.x && pt.x <= b.pMax.x &&
		pt.y >= b.pMin.y && pt.y <= b.pMax.y &&
		pt.z >= b.pMin.z && pt.z <= b.pMax.z)
}
func (b *BBox) Expand(delta float64) {
	b.pMin = Point{b.pMin.x - delta, b.pMin.y - delta, b.pMin.z - delta}
	b.pMax = Point{b.pMax.x + delta, b.pMax.y + delta, b.pMax.z + delta}
}
func (b *BBox) SurfaceArea() float64 {
	d := b.pMax.Sub(&b.pMin)
	return 2.0 * (d.x*d.y + d.x*d.z + d.y*d.z)
}
func (b *BBox) Volume() float64 {
	d := b.pMax.Sub(&b.pMin)
	return d.x * d.y * d.z
}
func (b *BBox) MaximumExtent() int {
	diag := b.pMax.Sub(&b.pMin)
	if diag.x > diag.y && diag.x > diag.z {
		return 0
	} else if diag.y > diag.z {
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
	return &Point{Lerp(tx, b.pMin.x, b.pMax.x), Lerp(ty, b.pMin.y, b.pMax.y), Lerp(tz, b.pMin.z, b.pMax.z)}
}
func (b *BBox) Offset(p *Point) *Vector {
	return &Vector{(p.x - b.pMin.x) / (b.pMax.x - b.pMin.x),
		(p.y - b.pMin.y) / (b.pMax.y - b.pMin.y),
		(p.z - b.pMin.z) / (b.pMax.z - b.pMin.z)}
}
func (b *BBox) BoundingSphere() (c *Point, rad float64) {
	c = b.pMin.Add(&Vector{b.pMax.x, b.pMax.y, b.pMax.z}).Scale(0.5)
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
	return fmt.Sprintf("bbox[ min: %v max: %v ]", b.pMin, b.pMax)
}
