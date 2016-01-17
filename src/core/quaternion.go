package pbrt

import (
	"math"
)

type Quaternion struct {
	v Vector
	w float64
}

func (q1 *Quaternion) Add(q2 *Quaternion) *Quaternion {
	return &Quaternion{*q1.v.Add(&q2.v), q1.w + q2.w}
}

func (q1 *Quaternion) Sub(q2 *Quaternion) *Quaternion {
	return &Quaternion{*q1.v.Sub(&q2.v), q1.w - q2.w}
}

func (q1 *Quaternion) Scale(s float64) *Quaternion {
	return &Quaternion{*q1.v.Scale(s), q1.w * s}
}
func (q1 *Quaternion) InvScale(invs float64) *Quaternion {
	s := 1.0 / invs
	return &Quaternion{*q1.v.Scale(s), q1.w * s}
}
func (q1 *Quaternion) ToTransform() *Transform {
	xx, yy, zz := q1.v.x*q1.v.x, q1.v.y*q1.v.y, q1.v.z*q1.v.z
	xy, xz, yz := q1.v.x*q1.v.y, q1.v.x*q1.v.z, q1.v.y*q1.v.z
	wx, wy, wz := q1.v.x*q1.w, q1.v.y*q1.w, q1.v.z*q1.w

	m := new(Matrix4x4)
	m.m[0][0] = 1.0 - 2.0*(yy+zz)
	m.m[0][1] = 2.0 * (xy + wz)
	m.m[0][2] = 2.0 * (xz - wy)
	m.m[1][0] = 2.0 * (xy - wz)
	m.m[1][1] = 1.0 - 2.0*(xx+zz)
	m.m[1][2] = 2.0 * (yz + wx)
	m.m[2][0] = 2.0 * (xz + wy)
	m.m[2][1] = 2.0 * (yz - wx)
	m.m[2][2] = 1.0 - 2.0*(xx+yy)

	// Transpose since we are left-handed.  Ugh.
	return CreateTransformExplicit(TransposeMatrix4x4(m), m)
}
func CreateQuaternionFromMatrix4x4(m *Matrix4x4) *Quaternion {
	q := new(Quaternion)
	trace := m.m[0][0] + m.m[1][1] + m.m[2][2]
	if trace > 0.0 {
		// Compute w from matrix trace, then xyz
		// 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
		s := math.Sqrt(trace + 1.0)
		q.w = s / 2.0
		s = 0.5 / s
		q.v.x = (m.m[2][1] - m.m[1][2]) * s
		q.v.y = (m.m[0][2] - m.m[2][0]) * s
		q.v.z = (m.m[1][0] - m.m[0][1]) * s
	} else {
		// Compute largest of $x$, $y$, or $z$, then remaining components
		nxt := []int{1, 2, 0}
		var qq [3]float64
		i := 0
		if m.m[1][1] > m.m[0][0] {
			i = 1
		}
		if m.m[2][2] > m.m[i][i] {
			i = 2
		}
		j := nxt[i]
		k := nxt[j]
		s := math.Sqrt((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + 1.0)
		qq[i] = s * 0.5
		if s != 0.0 {
			s = 0.5 / s
		}
		q.w = (m.m[k][j] - m.m[j][k]) * s
		qq[j] = (m.m[j][i] + m.m[i][j]) * s
		qq[k] = (m.m[k][i] + m.m[i][k]) * s
		q.v.x = qq[0]
		q.v.y = qq[1]
		q.v.z = qq[2]
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
	return DotVector(&q1.v, &q2.v) + q1.w*q2.w
}
func NormalizeQuaternion(q *Quaternion) *Quaternion {
	return q.InvScale(math.Sqrt(DotQuaternion(q, q)))
}
