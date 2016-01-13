package pbrt

import (
	"fmt"
	"math"
)

type Matrix4x4 struct {
	M [4][4]float64
}

func CreateIdentityMatrix4x4() *Matrix4x4 {
	mat := new(Matrix4x4)
	mat.M[0] = [4]float64{1.0, 0.0, 0.0, 0.0}
	mat.M[1] = [4]float64{0, 1, 0, 0}
	mat.M[2] = [4]float64{0, 0, 1, 0}
	mat.M[3] = [4]float64{0, 0, 0, 1}
	return mat
}

func CreateMatrix4x4(t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33 float64) *Matrix4x4 {
	mat := new(Matrix4x4)
	mat.M[0] = [4]float64{t00, t01, t02, t03}
	mat.M[1] = [4]float64{t10, t11, t12, t13}
	mat.M[2] = [4]float64{t20, t21, t22, t23}
	mat.M[3] = [4]float64{t30, t31, t32, t33}
	return mat
}

func CopyMatrix4x4(m *Matrix4x4) *Matrix4x4 {
	r := new(Matrix4x4)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r.M[i][j] = m.M[i][j]
		}
	}
	return r
}

func EqualMatrix4x4(m1, m2 *Matrix4x4) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if m1.M[i][j] != m2.M[i][j] {
				return false
			}
		}
	}
	return true
}

func NotEqualMatrix4x4(m1, m2 *Matrix4x4) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if m1.M[i][j] != m2.M[i][j] {
				return true
			}
		}
	}
	return false
}

func TransposeMatrix4x4(m *Matrix4x4) *Matrix4x4 {
	return CreateMatrix4x4(m.M[0][0], m.M[1][0], m.M[2][0], m.M[3][0],
		m.M[0][1], m.M[1][1], m.M[2][1], m.M[3][1],
		m.M[0][2], m.M[1][2], m.M[2][2], m.M[3][2],
		m.M[0][3], m.M[1][3], m.M[2][3], m.M[3][3])
}

func (m *Matrix4x4) String() string {
	s := "[ "
	for i := 0; i < 4; i++ {
		s = s + "  [ "
		for j := 0; j < 4; j++ {
			s = s + fmt.Sprint(s, "%f", m.M[i][j])
			if j != 3 {
				s = s + ", "
			}
		}
		s = s + " ]\n"
	}
	s = s + " ] "
	return s
}

func MulMatrix4x4(m1, m2 *Matrix4x4) *Matrix4x4 {
	r := new(Matrix4x4)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r.M[i][j] = m1.M[i][0]*m2.M[0][j] +
				m1.M[i][1]*m2.M[1][j] +
				m1.M[i][2]*m2.M[2][j] +
				m1.M[i][3]*m2.M[3][j]
		}
	}
	return r
}

func InverseMatrix4x4(m *Matrix4x4) (*Matrix4x4, error) {
	var indxc, indxr [4]int
	ipiv := []int{0, 0, 0, 0}
	minv := CopyMatrix4x4(m)
	for i := 0; i < 4; i++ {
		irow, icol := -1, -1
		big := 0.0
		// Choose pivot
		for j := 0; j < 4; j++ {
			if ipiv[j] != 1 {
				for k := 0; k < 4; k++ {
					if ipiv[k] == 0 {
						if math.Abs(minv.M[j][k]) >= big {
							big = math.Abs(minv.M[j][k])
							irow = j
							icol = k
						}
					} else if ipiv[k] > 1 {
						return nil, fmt.Errorf("Singular matrix in InverseMatrix4x4")
					}
				}
			}
		}
		ipiv[icol]++
		// Swap rows _irow_ and _icol_ for pivot
		if irow != icol {
			for k := 0; k < 4; k++ {
				minv.M[irow][k], minv.M[icol][k] = minv.M[icol][k], minv.M[irow][k]
			}
		}
		indxr[i] = irow
		indxc[i] = icol
		if minv.M[icol][icol] == 0.0 {
			return nil, fmt.Errorf("Singular matrix in InverseMatrix4x4")
		}

		// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
		pivinv := 1.0 / minv.M[icol][icol]
		minv.M[icol][icol] = 1.0
		for j := 0; j < 4; j++ {
			minv.M[icol][j] *= pivinv
		}
		// Subtract this row from others to zero out their columns
		for j := 0; j < 4; j++ {
			if j != icol {
				save := minv.M[j][icol]
				minv.M[j][icol] = 0
				for k := 0; k < 4; k++ {
					minv.M[j][k] -= minv.M[icol][k] * save
				}
			}
		}
	}
	// Swap columns to reflect permutation
	for j := 3; j >= 0; j-- {
		if indxr[j] != indxc[j] {
			for k := 0; k < 4; k++ {
				minv.M[k][indxr[j]], minv.M[k][indxc[j]] = minv.M[k][indxc[j]], minv.M[k][indxr[j]]
			}
		}
	}
	return minv, nil
}

type Transform struct {
	M, MInv *Matrix4x4
}

func CreateTransform(mat *Matrix4x4) (*Transform, error) {
	t := new(Transform)
	t.M = mat
	var err error
	t.MInv, err = InverseMatrix4x4(mat)
	if err != nil {
		return nil, err
	}
	return t, nil
}

func CreateTransformExplicit(mat, minv *Matrix4x4) *Transform {
	return &Transform{mat, minv}
}

func (t *Transform) String() string {
	return t.M.String()
}

func InverseTransform(t *Transform) *Transform {
	return &Transform{t.MInv, t.M}
}

func TransposeTransform(t *Transform) *Transform {
	return &Transform{TransposeMatrix4x4(t.M), TransposeMatrix4x4(t.MInv)}
}

func EqualTransform(t1, t2 *Transform) bool {
	return EqualMatrix4x4(t1.M, t2.M) && EqualMatrix4x4(t1.MInv, t2.MInv)
}

func NotEqualTransform(t1, t2 *Transform) bool {
	return NotEqualMatrix4x4(t1.M, t2.M) || NotEqualMatrix4x4(t1.MInv, t2.MInv)
}

func LessTransform(t1, t2 *Transform) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if t1.M.M[i][j] < t2.M.M[i][j] {
				return true
			}
			if t1.M.M[i][j] > t2.M.M[i][j] {
				return false
			}
		}
	}
	return false
}

func IsIdentityTransform(t *Transform) bool {
	return (t.M.M[0][0] == 1.0 && t.M.M[0][1] == 0.0 &&
		t.M.M[0][2] == 0.0 && t.M.M[0][3] == 0.0 &&
		t.M.M[1][0] == 0.0 && t.M.M[1][1] == 1.0 &&
		t.M.M[1][2] == 0.0 && t.M.M[1][3] == 0.0 &&
		t.M.M[2][0] == 0.0 && t.M.M[2][1] == 0.0 &&
		t.M.M[2][2] == 1.0 && t.M.M[2][3] == 0.0 &&
		t.M.M[3][0] == 0.0 && t.M.M[3][1] == 0.0 &&
		t.M.M[3][2] == 0.0 && t.M.M[3][3] == 1.0)
}

func HasScaleTransform(t *Transform) bool {
	la2 := VectorTransform(t, &Vector{1, 0, 0}).LengthSquared()
	lb2 := VectorTransform(t, &Vector{0, 1, 0}).LengthSquared()
	lc2 := VectorTransform(t, &Vector{0, 0, 1}).LengthSquared()
	return (la2 < 0.999 || la2 > 1.001 || lb2 < 0.999 || lb2 > 1.001 || lc2 < 0.999 || lc2 > 1.001)
}

func PointTransform(t *Transform, pt *Point) *Point {
	x, y, z := pt.X, pt.Y, pt.Z
	xp := t.M.M[0][0]*x + t.M.M[0][1]*y + t.M.M[0][2]*z + t.M.M[0][3]
	yp := t.M.M[1][0]*x + t.M.M[1][1]*y + t.M.M[1][2]*z + t.M.M[1][3]
	zp := t.M.M[2][0]*x + t.M.M[2][1]*y + t.M.M[2][2]*z + t.M.M[2][3]
	wp := t.M.M[3][0]*x + t.M.M[3][1]*y + t.M.M[3][2]*z + t.M.M[3][3]

	if wp == 1.0 {
		return &Point{xp, yp, zp}
	} else {
		return &Point{xp / wp, yp / wp, zp / wp}
	}
}

func VectorTransform(t *Transform, v *Vector) *Vector {
	x, y, z := v.X, v.Y, v.Z
	return &Vector{t.M.M[0][0]*x + t.M.M[0][1]*y + t.M.M[0][2]*z,
		t.M.M[1][0]*x + t.M.M[1][1]*y + t.M.M[1][2]*z,
		t.M.M[2][0]*x + t.M.M[2][1]*y + t.M.M[2][2]*z}
}

func NormalTransform(t *Transform, n *Normal) *Normal {
	x, y, z := n.X, n.Y, n.Z
	return &Normal{t.MInv.M[0][0]*x + t.MInv.M[1][0]*y + t.MInv.M[2][0]*z,
		t.MInv.M[0][1]*x + t.MInv.M[1][1]*y + t.MInv.M[2][1]*z,
		t.MInv.M[0][2]*x + t.MInv.M[1][2]*y + t.MInv.M[2][2]*z}
}

func RayTransform(t *Transform, r *Ray) *Ray {
	// TODO: implement this
	return nil
}

func RayDifferentialTransform(t *Transform, r *RayDifferential) *RayDifferential {
	// TODO: implement this
	return nil
}

func BBoxTransform(t *Transform, b *BBox) *BBox {
	// TODO: implement this
	return nil
}

func (t1 *Transform) MultTransform(t2 *Transform) *Transform {
	m1 := MulMatrix4x4(t1.M, t2.M)
	m2 := MulMatrix4x4(t2.MInv, t1.MInv)
	return &Transform{m1, m2}
}

func SwapsHandednessTransform(t *Transform) bool {
	det := ((t.M.M[0][0] *
		(t.M.M[1][1]*t.M.M[2][2] -
			t.M.M[1][2]*t.M.M[2][1])) -
		(t.M.M[0][1] *
			(t.M.M[1][0]*t.M.M[2][2] -
				t.M.M[1][2]*t.M.M[2][0])) +
		(t.M.M[0][2] *
			(t.M.M[1][0]*t.M.M[2][1] -
				t.M.M[1][1]*t.M.M[2][0])))
	return det < 0.0
}

func TranslateTransform(delta *Vector) *Transform {
	m := CreateMatrix4x4(1, 0, 0, delta.X,
		0, 1, 0, delta.Y,
		0, 0, 1, delta.Z,
		0, 0, 0, 1)
	minv := CreateMatrix4x4(1, 0, 0, -delta.X,
		0, 1, 0, -delta.Y,
		0, 0, 1, -delta.Z,
		0, 0, 0, 1)
	return &Transform{m, minv}
}

func ScaleTransform(x, y, z float64) *Transform {
	m := CreateMatrix4x4(x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1)
	minv := CreateMatrix4x4(1.0/x, 0, 0, 0,
		0, 1.0/y, 0, 0,
		0, 0, 1.0/z, 0,
		0, 0, 0, 1)
	return &Transform{m, minv}
}

func RotateXTransform(angle float64) *Transform {
	sin_t := math.Sin(Radians(angle))
	cos_t := math.Cos(Radians(angle))
	m := CreateMatrix4x4(1, 0, 0, 0,
		0, cos_t, -sin_t, 0,
		0, sin_t, cos_t, 0,
		0, 0, 0, 1)
	return &Transform{m, TransposeMatrix4x4(m)}
}

func RotateYTransform(angle float64) *Transform {
	sin_t := math.Sin(Radians(angle))
	cos_t := math.Cos(Radians(angle))
	m := CreateMatrix4x4(cos_t, 0, sin_t, 0,
		0, 1, 0, 0,
		-sin_t, 0, cos_t, 0,
		0, 0, 0, 1)
	return &Transform{m, TransposeMatrix4x4(m)}
}

func RotateZTransform(angle float64) *Transform {
	sin_t := math.Sin(Radians(angle))
	cos_t := math.Cos(Radians(angle))
	m := CreateMatrix4x4(cos_t, -sin_t, 0, 0,
		sin_t, cos_t, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1)
	return &Transform{m, TransposeMatrix4x4(m)}
}

func RotateTransform(angle float64, axis *Vector) *Transform {
	a := NormalizeVector(axis)
	s := math.Sin(Radians(angle))
	c := math.Cos(Radians(angle))

	mat := CreateMatrix4x4((a.X*a.X + (1.0-a.X*a.X)*c), (a.X*a.Y*(1.0-c) - a.Z*s), (a.X*a.Z*(1.0-c) + a.Y*s), 0,
		(a.X*a.Y*(1.0-c) + a.Z*s), (a.Y*a.Y + (1.0-a.Y*a.Y)*c), (a.Y*a.Z*(1.0-c) - a.X*s), 0,
		(a.X*a.Z*(1.0-c) - a.Y*s), (a.Y*a.Z*(1.0-c) + a.X*s), (a.Z*a.Z + (1.0-a.Z*a.Z)*c), 0,
		0, 0, 0, 1)

	return &Transform{mat, TransposeMatrix4x4(mat)}
}

func LookAtTransform(pos, look *Point, up *Vector) (*Transform, error) {
	// Initialize first three columns of viewing matrix
	dir := NormalizeVector(look.Sub(pos))
	if CrossVector(NormalizeVector(up), dir).Length() == 0 {
		return nil, fmt.Errorf("\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) passed to LookAt are pointing in the same direction.  Using the identity transformation.", up.X, up.Y, up.Z, dir.X, dir.Y, dir.Z)
	}

	camToWorld := new(Matrix4x4)
	// Initialize fourth column of viewing matrix
	camToWorld.M[0][3] = pos.X
	camToWorld.M[1][3] = pos.Y
	camToWorld.M[2][3] = pos.Z
	camToWorld.M[3][3] = 1

	left := NormalizeVector(CrossVector(NormalizeVector(up), dir))
	newUp := CrossVector(dir, left)
	camToWorld.M[0][0] = left.X
	camToWorld.M[1][0] = left.Y
	camToWorld.M[2][0] = left.Z
	camToWorld.M[3][0] = 0.0
	camToWorld.M[0][1] = newUp.X
	camToWorld.M[1][1] = newUp.Y
	camToWorld.M[2][1] = newUp.Z
	camToWorld.M[3][1] = 0.0
	camToWorld.M[0][2] = dir.X
	camToWorld.M[1][2] = dir.Y
	camToWorld.M[2][2] = dir.Z
	camToWorld.M[3][2] = 0.0
	worldToCam, err := InverseMatrix4x4(camToWorld)
	if err != nil {
		return nil, fmt.Errorf("CameraToWorld matrix is not invertable.")
	}
	return &Transform{worldToCam, camToWorld}, nil
}

func OrthographicTransform(znear, zfar float64) *Transform {
	return ScaleTransform(1.0, 1.0, 1.0/(zfar-znear)).MultTransform(TranslateTransform(&Vector{0.0, 0.0, -znear}))
}

func PerspectiveTransform(fov, znear, zfar float64) *Transform {
	// Perform projective divide
	persp := CreateMatrix4x4(1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, zfar/(zfar-znear), -zfar*znear/(zfar-znear),
		0, 0, 1, 0)

	// Scale to canonical viewing volume
	invTanAng := 1.0 / math.Tan(Radians(fov)/2.0)
	perspTrans, _ := CreateTransform(persp)
	return ScaleTransform(invTanAng, invTanAng, 1.0).MultTransform(perspTrans)
}

type AnimatedTransform struct {
	startTime, endTime           float64
	startTransform, endTransform *Transform
	actuallyAnimated             bool
	t                            [2]*Vector
	r                            [2]*Quaternion
	s                            [2]*Matrix4x4
}

func CreateAnimatedTransform(trans1 *Transform, time1 float64, trans2 *Transform, time2 float64) *AnimatedTransform {
	t := new(AnimatedTransform)
	t.startTime = time1
	t.startTransform = trans1
	t.endTime = time2
	t.endTransform = trans2
	t.actuallyAnimated = NotEqualTransform(trans1, trans2)
	t.t[0], t.r[0], t.s[0] = DecomposeMatrix4x4(t.startTransform.M)
	t.t[1], t.r[1], t.s[1] = DecomposeMatrix4x4(t.endTransform.M)

	return t
}

func DecomposeMatrix4x4(m *Matrix4x4) (*Vector, *Quaternion, *Matrix4x4) {

	// Extract translation _T_ from transformation matrix
	T := &Vector{m.M[0][3], m.M[1][3], m.M[2][3]}

	// Compute new transformation matrix _M_ without translation
	M := CopyMatrix4x4(m)
	for i := 0; i < 3; i++ {
		M.M[i][3] = 0.0
		M.M[3][i] = 0.0
	}
	M.M[3][3] = 1.0

	// Extract rotation _R_ from transformation matrix
	norm := 1.0
	count := 0
	R := CopyMatrix4x4(M)
	for count < 100 && norm > .0001 {
		// Compute next matrix _Rnext_ in series
		Rnext := new(Matrix4x4)
		Rit, _ := InverseMatrix4x4(TransposeMatrix4x4(R))
		for i := 0; i < 4; i++ {
			for j := 0; j < 4; j++ {
				Rnext.M[i][j] = 0.5 * (R.M[i][j] + Rit.M[i][j])
			}
		}
		// Compute norm of difference between _R_ and _Rnext_
		norm = 0.0
		for i := 0; i < 3; i++ {
			n := math.Abs(R.M[i][0]-Rnext.M[i][0]) +
				math.Abs(R.M[i][1]-Rnext.M[i][1]) +
				math.Abs(R.M[i][2]-Rnext.M[i][2])
			norm = math.Max(norm, n)
		}
		R = Rnext
		count++
	}
	// XXX TODO FIXME deal with flip...
	Rquat := CreateQuaternionFromMatrix4x4(R)

	// Compute scale _S_ using rotation and original matrix
	Rinv, _ := InverseMatrix4x4(R)
	S := MulMatrix4x4(Rinv, M)

	return T, Rquat, S
}

func (t *AnimatedTransform) Interpolate(time float64) *Transform {
	// Handle boundary conditions for matrix interpolation
	if !t.actuallyAnimated || time <= t.startTime {
		return t.startTransform
	}
	if time >= t.endTime {
		return t.endTransform
	}
	dt := (time - t.startTime) / (t.endTime - t.startTime)
	// Interpolate translation at _dt_
	trans := t.t[0].Scale(1.0 - dt).Add(t.t[1].Scale(dt))

	// Interpolate rotation at _dt_
	rotate := Slerp(dt, t.r[0], t.r[1])

	// Interpolate scale at _dt_
	scale := new(Matrix4x4)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			scale.M[i][j] = Lerp(dt, t.s[0].M[i][j], t.s[1].M[i][j])
		}
	}
	// Compute interpolated matrix as product of interpolated components
	scaleTrans, _ := CreateTransform(scale)
	rotTrans := rotate.ToTransform()
	return TranslateTransform(trans).MultTransform(rotTrans.MultTransform(scaleTrans))
}

func VectorAnimatedTransform(t *AnimatedTransform, time float64, v *Vector) *Vector {
	// TODO: implement this
	return nil
}

func PointAnimatedTransform(t *AnimatedTransform, time float64, p *Point) *Point {
	// TODO: implement this
	return nil
}

func RayAnimatedTransform(t *AnimatedTransform, r *Ray) *Ray {
	// TODO: implement this
	return nil
}

func RayDifferentialAnimatedTransform(t *AnimatedTransform, r *RayDifferential) *RayDifferential {
	// TODO: implement this
	return nil
}

func MotionBoundsAnimatedTransform(t *AnimatedTransform, b *BBox, useInverse bool) *BBox {
	if !t.actuallyAnimated {
		if useInverse {
			return BBoxTransform(InverseTransform(t.startTransform), b)
		} else {
			return BBoxTransform(t.startTransform, b)
		}
	}
	ret := CreateEmptyBBox()
	nSteps := 128
	for i := 0; i < nSteps; i++ {
		time := Lerp(float64(i)/float64(nSteps-1), t.startTime, t.endTime)
		trans := t.Interpolate(time)
		if useInverse {
			trans = InverseTransform(trans)
		}
		ret = UnionBBoxes(ret, BBoxTransform(trans, b))
	}
	return ret
}

func HasScaleAnimatedTransform(t *AnimatedTransform) bool {
	return HasScaleTransform(t.startTransform) || HasScaleTransform(t.endTransform)
}

// Matrix4x4 Method Definitions
func SolveLinearSystem2x2(A [2][2]float64, B [2]float64) (ok bool, x0, x1 float64) {
	det := A[0][0]*A[1][1] - A[0][1]*A[1][0]
	if math.Abs(det) < 1.0e-10 {
		return false, 0.0, 0.0
	}
	x0 = (A[1][1]*B[0] - A[0][1]*B[1]) / det
	x1 = (A[0][0]*B[1] - A[1][0]*B[0]) / det
	if math.IsNaN(x0) || math.IsNaN(x1) {
		return false, 0.0, 0.0
	}
	return true, x0, x1
}
