package pbrt

import (
	"math"
)

type (
	TextureMapping2D interface {
		Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64)
	}

	UVMapping2D struct {
		su, sv, du, dv float64
	}

	SphericalMapping2D struct {
		worldToTexture *Transform
	}

	CylindricalMapping2D struct {
		worldToTexture *Transform
	}

	PlanarMapping2D struct {
		vs, vt *Vector
		ds, dt float64
	}

	TextureMapping3D interface {
		Map(dg *DifferentialGeometry) (p *Point, dpdx, dpdy *Vector)
	}

	IdentityMapping3D struct {
		worldToTexture *Transform
	}

	Texture interface {
		Evaluate(dg *DifferentialGeometry) float64
	}

	TextureFloat struct {
	}

	TextureSpectrum struct {
	}
)

func (m *UVMapping2D) Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64) {
	s = m.su*dg.u + m.du
	t = m.sv*dg.v + m.dv
	// Compute texture differentials for 2D identity mapping
	dsdx = m.su * dg.dudx
	dtdx = m.sv * dg.dvdx
	dsdy = m.su * dg.dudy
	dtdy = m.sv * dg.dvdy
	return s, t, dsdx, dtdx, dsdy, dtdy
}

func (m *SphericalMapping2D) Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64) {
	s, t = m.sphere(dg.p)
	// Compute texture coordinate differentials for sphere $(u,v)$ mapping
	delta := 0.1
	sx, tx := m.sphere(dg.p.Add(dg.dpdx.Scale(delta)))
	dsdx = (sx - s) / delta
	dtdx = (tx - t) / delta
	if dtdx > 0.5 {
		dtdx = 1.0 - dtdx
	} else if dtdx < -.50 {
		dtdx = -(dtdx + 1)
	}
	sy, ty := m.sphere(dg.p.Add(dg.dpdy.Scale(delta)))
	dsdy = (sy - s) / delta
	dtdy = (ty - t) / delta
	if dtdy > 0.5 {
		dtdy = 1.0 - dtdy
	} else if dtdy < -0.5 {
		dtdy = -(dtdy + 1)
	}
	return s, t, dsdx, dtdx, dsdy, dtdy
}

func (m *SphericalMapping2D) sphere(p *Point) (s, t float64) {
	vec := NormalizeVector(PointTransform(m.worldToTexture, p).Sub(CreatePoint(0, 0, 0)))
	theta := SphericalTheta(vec)
	phi := SphericalPhi(vec)
	s = theta / math.Pi
	t = phi / (2.0 * math.Pi)
	return s, t
}

func (m *CylindricalMapping2D) Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64) {
	s, t = m.cylinder(dg.p)
	// Compute texture coordinate differentials for cylinder $(u,v)$ mapping
	delta := 0.01
	sx, tx := m.cylinder(dg.p.Add(dg.dpdx.Scale(delta)))
	dsdx = (sx - s) / delta
	dtdx = (tx - t) / delta
	if dtdx > 0.5 {
		dtdx = 1.0 - dtdx
	} else if dtdx < -0.5 {
		dtdx = -(dtdx + 1)
	}
	sy, ty := m.cylinder(dg.p.Add(dg.dpdy.Scale(delta)))
	dsdy = (sy - s) / delta
	dtdy = (ty - t) / delta
	if dtdy > 0.5 {
		dtdy = 1.0 - dtdy
	} else if dtdy < -0.5 {
		dtdy = -(dtdy + 1)
	}
	return s, t, dsdx, dtdx, dsdy, dtdy
}

func (m *CylindricalMapping2D) cylinder(p *Point) (s, t float64) {
	vec := NormalizeVector(PointTransform(m.worldToTexture, p).Sub(CreatePoint(0, 0, 0)))
	s = (math.Pi + math.Atan2(vec.y, vec.x)) / (2.0 * math.Pi)
	t = vec.z
	return s, t
}

func (m *PlanarMapping2D) Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64) {
	vec := dg.p.Sub(CreatePoint(0, 0, 0))
	s = m.ds + DotVector(vec, m.vs)
	t = m.dt + DotVector(vec, m.vt)
	dsdx = DotVector(dg.dpdx, m.vs)
	dtdx = DotVector(dg.dpdx, m.vt)
	dsdy = DotVector(dg.dpdy, m.vs)
	dtdy = DotVector(dg.dpdy, m.vt)
	return s, t, dsdx, dtdx, dsdy, dtdy
}

func (m *IdentityMapping3D) Map(dg *DifferentialGeometry) (p *Point, dpdx, dpdy *Vector) {
    return PointTransform(m.worldToTexture, dg.p), VectorTransform(m.worldToTexture, dg.dpdx), VectorTransform(m.worldToTexture, dg.dpdy)
}

func (t *TextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *TextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	return CreateSpectrum1(0.0)
}

// Texture Inline Functions
func SmoothStep(min, max, value float64) float64 {
    v := Clamp((value - min) / (max - min), 0.0, 1.0)
    return v * v * (-2.0 * v  + 3.0)
}

func Lanczos(x, tau float64) float64 {
    x = math.Abs(x)
    if x < 1.0e-5 { return 1.0 }
    if x > 1.0  { return 0.0 }
    x *= math.Pi
    s := math.Sin(x * tau) / (x * tau)
    lanczos := math.Sin(x) / x
    return s * lanczos
}

func Noise(x, y, z float64) float64 {
	return 0.0
}

func NoiseP(p *Point) float64 {
	return 0.0
}

func FBm(p *Point, dpdx, dpdy *Vector, omega float64, maxOctaves int) float64 {
    // Compute number of octaves for antialiased FBm
    s2 := math.Max(dpdx.LengthSquared(), dpdy.LengthSquared())
    foctaves := math.Min(float64(maxOctaves), math.Max(0.0, -1.0 - 0.5 * Log2(s2)))
    octaves := Floor2Int(foctaves)

    // Compute sum of octaves of noise for FBm
    sum, lambda, o := 0.0, 1.0, 1.0
    for i := 0; i < octaves; i++ {
        sum += o * NoiseP(p.Scale(lambda))
        lambda *= 1.99
        o *= omega
    }
    partialOctave := foctaves - float64(octaves)
    sum += o * SmoothStep(0.3, 0.7, partialOctave) * NoiseP(p.Scale(lambda))
    return sum
}

func Turbulence(p *Point, dpdx, dpdy *Vector, omega float64, maxOctaves int) float64 {
    // Compute number of octaves for antialiased FBm
    s2 := math.Max(dpdx.LengthSquared(), dpdy.LengthSquared())
    foctaves := math.Min(float64(maxOctaves), math.Max(0.0, -1.0 - 0.5 * Log2(s2)))
    octaves := Floor2Int(foctaves)

    // Compute sum of octaves of noise for turbulence
    sum, lambda, o := 0.0, 1.0, 1.0
    for i := 0; i < octaves; i++ {
        sum += o * math.Abs(NoiseP(p.Scale(lambda)))
        lambda *= 1.99
        o *= omega
    }
    partialOctave := foctaves - float64(octaves)
    sum += o * SmoothStep(0.3, 0.7, partialOctave) * math.Abs(NoiseP(p.Scale(lambda)))

    // finally, add in value to account for average value of fabsf(Noise())
    // (~0.2) for the remaining octaves...
    sum += (float64(maxOctaves) - foctaves) * 0.2

    return sum
}
