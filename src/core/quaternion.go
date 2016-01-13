package pbrt

import (
	"math"
)

type Quaternion struct {
	V Vector
	W float64
}

func (q1 *Quaternion) Add(q2 *Quaternion) *Quaternion {
	return &Quaternion{*q1.V.Add(&q2.V), q1.W + q2.W}
}

func (q1 *Quaternion) Sub(q2 *Quaternion) *Quaternion {
	return &Quaternion{*q1.V.Sub(&q2.V), q1.W - q2.W}
}

func (q1 *Quaternion) Scale(s float64) *Quaternion {
	return &Quaternion{*q1.V.Scale(s), q1.W * s}
}
func (q1 *Quaternion) InvScale(invs float64) *Quaternion {
	s := 1.0 / invs
	return &Quaternion{*q1.V.Scale(s), q1.W * s}
}
func (q1 *Quaternion) ToTransform() *Transform {
	xx, yy, zz := q1.V.X*q1.V.X, q1.V.Y*q1.V.Y, q1.V.Z*q1.V.Z
	xy, xz, yz := q1.V.X*q1.V.Y, q1.V.X*q1.V.Z, q1.V.Y*q1.V.Z
	wx, wy, wz := q1.V.X*q1.W, q1.V.Y*q1.W, q1.V.Z*q1.W

	m := new(Matrix4x4)
	m.M[0][0] = 1.0 - 2.0*(yy+zz)
	m.M[0][1] = 2.0 * (xy + wz)
	m.M[0][2] = 2.0 * (xz - wy)
	m.M[1][0] = 2.0 * (xy - wz)
	m.M[1][1] = 1.0 - 2.0*(xx+zz)
	m.M[1][2] = 2.0 * (yz + wx)
	m.M[2][0] = 2.0 * (xz + wy)
	m.M[2][1] = 2.0 * (yz - wx)
	m.M[2][2] = 1.0 - 2.0*(xx+yy)

	// Transpose since we are left-handed.  Ugh.
	return CreateTransformExplicit(TransposeMatrix4x4(m), m)
}
func CreateQuaternionFromMatrix4x4(m *Matrix4x4) *Quaternion {
	q := new(Quaternion)
	trace := m.M[0][0] + m.M[1][1] + m.M[2][2]
	if trace > 0.0 {
		// Compute w from matrix trace, then xyz
		// 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
		s := math.Sqrt(trace + 1.0)
		q.W = s / 2.0
		s = 0.5 / s
		q.V.X = (m.M[2][1] - m.M[1][2]) * s
		q.V.Y = (m.M[0][2] - m.M[2][0]) * s
		q.V.Z = (m.M[1][0] - m.M[0][1]) * s
	} else {
		// Compute largest of $x$, $y$, or $z$, then remaining components
		nxt := []int{1, 2, 0}
		var qq [3]float64
		i := 0
		if m.M[1][1] > m.M[0][0] {
			i = 1
		}
		if m.M[2][2] > m.M[i][i] {
			i = 2
		}
		j := nxt[i]
		k := nxt[j]
		s := math.Sqrt((m.M[i][i] - (m.M[j][j] + m.M[k][k])) + 1.0)
		qq[i] = s * 0.5
		if s != 0.0 {
			s = 0.5 / s
		}
		q.W = (m.M[k][j] - m.M[j][k]) * s
		qq[j] = (m.M[j][i] + m.M[i][j]) * s
		qq[k] = (m.M[k][i] + m.M[i][k]) * s
		q.V.X = qq[0]
		q.V.Y = qq[1]
		q.V.Z = qq[2]
	}
	return q
}
func Slerp(t float64, q1, q2 *Quaternion) *Quaternion {
	cosTheta := DotQuaternion(q1, q2)
	if cosTheta > .9995 {
		return NormalizeQuaternion(q1.Scale(1.0 - t).Add(q2.Scale(t)))
	} else {
		theta := math.Acos(Clamp(cosTheta, -1.0, 1.0))
		thetap := theta * t
		qperp := NormalizeQuaternion(q2.Sub(q1.Scale(cosTheta)))
		return q1.Scale(math.Cos(thetap)).Add(qperp.Scale(math.Sin(thetap)))
	}

}
func DotQuaternion(q1, q2 *Quaternion) float64 {
	return DotVector(&q1.V, &q2.V) + q1.W*q2.W
}
func NormalizeQuaternion(q *Quaternion) *Quaternion {
	return q.InvScale(math.Sqrt(DotQuaternion(q, q)))
}
