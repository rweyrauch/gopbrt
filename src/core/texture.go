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
	"strings"
)

const (
	AA_NONE AntialiasMethod = iota
	AA_CLOSEDFORM
)

type AntialiasMethod int

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
		vs, vt Vector
		ds, dt float64
	}

	TextureMapping3D interface {
		Map(dg *DifferentialGeometry) (p *Point, dpdx, dpdy *Vector)
	}

	IdentityMapping3D struct {
		worldToTexture *Transform
	}
)

func NewUVMapping2D(su, sv, du, dv float64) *UVMapping2D {
	return &UVMapping2D{su, sv, du, dv}
}

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

func NewSphericalMapping2D(worldToTexture *Transform) *SphericalMapping2D {
	return &SphericalMapping2D{worldToTexture}
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

func NewCylindricalMapping2D(worldToTexture *Transform) *CylindricalMapping2D {
	return &CylindricalMapping2D{worldToTexture}
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
	s = (math.Pi + math.Atan2(vec.Y, vec.X)) / (2.0 * math.Pi)
	t = vec.Z
	return s, t
}

func NewPlanarMapping2D(vs, vt Vector, ds, dt float64) *PlanarMapping2D {
	return &PlanarMapping2D{vs, vt, ds, dt}
}

func (m *PlanarMapping2D) Map(dg *DifferentialGeometry) (s, t, dsdx, dtdx, dsdy, dtdy float64) {
	vec := dg.p.Sub(CreatePoint(0, 0, 0))
	s = m.ds + DotVector(vec, &m.vs)
	t = m.dt + DotVector(vec, &m.vt)
	dsdx = DotVector(dg.dpdx, &m.vs)
	dtdx = DotVector(dg.dpdx, &m.vt)
	dsdy = DotVector(dg.dpdy, &m.vs)
	dtdy = DotVector(dg.dpdy, &m.vt)
	return s, t, dsdx, dtdx, dsdy, dtdy
}

func NewIdentityMapping3D(worldToTexture *Transform) *IdentityMapping3D {
	return &IdentityMapping3D{worldToTexture}
}

func (m *IdentityMapping3D) Map(dg *DifferentialGeometry) (p *Point, dpdx, dpdy *Vector) {
	return PointTransform(m.worldToTexture, dg.p), VectorTransform(m.worldToTexture, dg.dpdx), VectorTransform(m.worldToTexture, dg.dpdy)
}

// Texture Inline Functions
func SmoothStep(min, max, value float64) float64 {
	v := Clamp((value-min)/(max-min), 0.0, 1.0)
	return v * v * (-2.0*v + 3.0)
}

func Lanczos(x, tau float64) float64 {
	x = math.Abs(x)
	if x < 1.0e-5 {
		return 1.0
	}
	if x > 1.0 {
		return 0.0
	}
	x *= math.Pi
	s := math.Sin(x*tau) / (x * tau)
	lanczos := math.Sin(x) / x
	return s * lanczos
}

const (
	NOISE_PERM_SIZE = 256
)

var NoisePerm [NOISE_PERM_SIZE * 2]int = [NOISE_PERM_SIZE * 2]int{
	151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96,
	53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
	// Remainder of the noise permutation table
	8, 99, 37, 240, 21, 10, 23,
	190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
	88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
	77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
	102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
	135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
	5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
	223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
	251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
	49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
	138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
	151, 160, 137, 91, 90, 15,
	131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
	190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
	88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
	77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
	102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
	135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
	5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
	223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
	251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
	49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
	138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
}

func Grad(x, y, z int, dx, dy, dz float64) float64 {
	h := NoisePerm[NoisePerm[NoisePerm[x]+y]+z]
	h &= 15
	u := dy
	if h < 8 || h == 12 || h == 13 {
		u = dx
	}
	v := dz
	if h < 4 || h == 12 || h == 13 {
		v = dy
	}
	var n float64
	if h&1 != 0 {
		n = -u
	} else {
		n = u
	}
	if h&2 != 0 {
		n -= v
	} else {
		n += v
	}
	return n
}

func NoiseWeight(t float64) float64 {
	t3 := t * t * t
	t4 := t3 * t
	return 6.0*t4*t - 15.0*t4 + 10.0*t3
}

func Noise(x, y, z float64) float64 {
	// Compute noise cell coordinates and offsets
	ix, iy, iz := Floor2Int(x), Floor2Int(y), Floor2Int(z)
	dx, dy, dz := x-float64(ix), y-float64(iy), z-float64(iz)

	// Compute gradient weights
	ix &= (NOISE_PERM_SIZE - 1)
	iy &= (NOISE_PERM_SIZE - 1)
	iz &= (NOISE_PERM_SIZE - 1)
	w000 := Grad(ix, iy, iz, dx, dy, dz)
	w100 := Grad(ix+1, iy, iz, dx-1, dy, dz)
	w010 := Grad(ix, iy+1, iz, dx, dy-1, dz)
	w110 := Grad(ix+1, iy+1, iz, dx-1, dy-1, dz)
	w001 := Grad(ix, iy, iz+1, dx, dy, dz-1)
	w101 := Grad(ix+1, iy, iz+1, dx-1, dy, dz-1)
	w011 := Grad(ix, iy+1, iz+1, dx, dy-1, dz-1)
	w111 := Grad(ix+1, iy+1, iz+1, dx-1, dy-1, dz-1)

	// Compute trilinear interpolation of weights
	wx, wy, wz := NoiseWeight(dx), NoiseWeight(dy), NoiseWeight(dz)
	x00 := Lerp(wx, w000, w100)
	x10 := Lerp(wx, w010, w110)
	x01 := Lerp(wx, w001, w101)
	x11 := Lerp(wx, w011, w111)
	y0 := Lerp(wy, x00, x10)
	y1 := Lerp(wy, x01, x11)
	return Lerp(wz, y0, y1)
}

func NoiseP(p *Point) float64 {
	return Noise(p.X, p.Y, p.Z)
}

func FBm(p *Point, dpdx, dpdy *Vector, omega float64, maxOctaves int) float64 {
	// Compute number of octaves for antialiased FBm
	s2 := math.Max(dpdx.LengthSquared(), dpdy.LengthSquared())
	foctaves := math.Min(float64(maxOctaves), math.Max(0.0, -1.0-0.5*Log2(s2)))
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
	foctaves := math.Min(float64(maxOctaves), math.Max(0.0, -1.0-0.5*Log2(s2)))
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

type (
	TextureFloat interface {
		Evaluate(dg *DifferentialGeometry) float64
	}

	TextureSpectrum interface {
		Evaluate(dg *DifferentialGeometry) *Spectrum
	}

	BilerpTextureFloat struct {
		mapping            TextureMapping2D
		v00, v01, v10, v11 float64
	}
	BilerpTextureSpectrum struct {
		mapping            TextureMapping2D
		v00, v01, v10, v11 Spectrum
	}

	Checkerboard2DTextureFloat struct {
		mapping    TextureMapping2D
		tex1, tex2 TextureFloat
		aaMethod   AntialiasMethod
	}
	Checkerboard2DTextureSpectrum struct {
		mapping    TextureMapping2D
		tex1, tex2 TextureSpectrum
		aaMethod   AntialiasMethod
	}
	Checkerboard3DTextureFloat struct {
		mapping    TextureMapping3D
		tex1, tex2 TextureFloat
	}
	Checkerboard3DTextureSpectrum struct {
		mapping    TextureMapping3D
		tex1, tex2 TextureSpectrum
	}

	ConstantTextureFloat struct {
		value float64
	}
	ConstantTextureSpectrum struct {
		value Spectrum
	}

	DotsTextureFloat struct {
		mapping               TextureMapping2D
		outsideDot, insideDot TextureFloat
	}
	DotsTextureSpectrum struct {
		mapping               TextureMapping2D
		outsideDot, insideDot TextureSpectrum
	}

	FBmTextureFloat struct {
		omega   float64
		octaves int
		mapping TextureMapping3D
	}
	FBmTextureSpectrum struct {
		omega   float64
		octaves int
		mapping TextureMapping3D
	}

	ImageTextureFloat struct {
		mapping TextureMapping2D
		mipmap  *MIPMapFloat
	}
	ImageTextureSpectrum struct {
		mapping TextureMapping2D
		mipmap  *MIPMapSpectrum
	}

	MarbleTextureFloat struct {
		mapping                 TextureMapping3D
		octaves                 int
		omega, scale, variation float64
	}
	MarbleTextureSpectrum struct {
		mapping                 TextureMapping3D
		octaves                 int
		omega, scale, variation float64
	}

	MixTextureFloat struct {
		tex1, tex2 TextureFloat
		amount     TextureFloat
	}
	MixTextureSpectrum struct {
		tex1, tex2 TextureSpectrum
		amount     TextureFloat
	}

	ScaleTextureFloat struct {
		tex1, tex2 TextureFloat
	}
	ScaleTextureSpectrum struct {
		tex1, tex2 TextureSpectrum
	}

	UVTextureFloat struct {
		mapping TextureMapping2D
	}
	UVTextureSpectrum struct {
		mapping TextureMapping2D
	}

	WindyTextureFloat struct {
		mapping TextureMapping3D
	}
	WindyTextureSpectrum struct {
		mapping TextureMapping3D
	}

	WrinkledTextureFloat struct {
		octaves int
		omega   float64
		mapping TextureMapping3D
	}
	WrinkledTextureSpectrum struct {
		octaves int
		omega   float64
		mapping TextureMapping3D
	}
)

func NewBilerpTextureFloat(mapping TextureMapping2D, v00, v01, v10, v11 float64) *BilerpTextureFloat {
	return &BilerpTextureFloat{mapping, v00, v01, v10, v11}
}
func (tex *BilerpTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	s, t, _, _, _, _ := tex.mapping.Map(dg)
	return (1-s)*(1-t)*tex.v00 + (1-s)*t*tex.v01 +
		s*(1-t)*tex.v10 + s*t*tex.v11
}

func NewBilerpTextureSpectrum(mapping TextureMapping2D, v00, v01, v10, v11 Spectrum) *BilerpTextureSpectrum {
	return &BilerpTextureSpectrum{mapping, v00, v01, v10, v11}
}
func (tex *BilerpTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	s, t, _, _, _, _ := tex.mapping.Map(dg)
	return tex.v00.Scale((1 - s) * (1 - t)).Add(tex.v01.Scale((1 - s) * t)).Add(tex.v10.Scale(s * (1 - t))).Add(tex.v11.Scale(s * t))
}

func NewCheckerboard2DTextureFloat(mapping TextureMapping2D, tex1, tex2 TextureFloat, aamode string) *Checkerboard2DTextureFloat {
	aaMethod := AA_NONE
	if strings.Compare(aamode, "closedform") == 0 {
		aaMethod = AA_CLOSEDFORM
	}
	return &Checkerboard2DTextureFloat{mapping, tex1, tex2, aaMethod}
}

func (tex *Checkerboard2DTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	s, t, dsdx, dtdx, dsdy, dtdy := tex.mapping.Map(dg)
	if tex.aaMethod == AA_NONE {
		// Point sample _Checkerboard2DTexture_
		if (Floor2Int(s)+Floor2Int(t))%2 == 0 {
			return tex.tex1.Evaluate(dg)
		}
		return tex.tex2.Evaluate(dg)
	} else {
		// Compute closed-form box-filtered _Checkerboard2DTexture_ value

		// Evaluate single check if filter is entirely inside one of them
		ds := math.Max(math.Abs(dsdx), math.Abs(dsdy))
		dt := math.Max(math.Abs(dtdx), math.Abs(dtdy))
		s0, s1 := s-ds, s+ds
		t0, t1 := t-dt, t+dt
		if Floor2Int(s0) == Floor2Int(s1) && Floor2Int(t0) == Floor2Int(t1) {
			// Point sample _Checkerboard2DTexture_
			if (Floor2Int(s)+Floor2Int(t))%2 == 0 {
				return tex.tex1.Evaluate(dg)
			}
			return tex.tex2.Evaluate(dg)
		}

		// Apply box filter to checkerboard region
		BUMPINT := func(x float64) float64 {
			return (float64(Floor2Int((x)/2)) + 2.0*math.Max((x/2)-float64(Floor2Int(x/2))-0.5, 0.0))
		}
		sint := (BUMPINT(s1) - BUMPINT(s0)) / (2.0 * ds)
		tint := (BUMPINT(t1) - BUMPINT(t0)) / (2.0 * dt)
		area2 := sint + tint - 2.0*sint*tint
		if ds > 1.0 || dt > 1.0 {
			area2 = 0.5
		}
		return (1.0-area2)*tex.tex1.Evaluate(dg) + area2*tex.tex2.Evaluate(dg)
	}
}

func NewCheckerboard2DTextureSpectrum(mapping TextureMapping2D, tex1, tex2 TextureSpectrum, aamode string) *Checkerboard2DTextureSpectrum {
	aaMethod := AA_NONE
	if strings.Compare(aamode, "closedform") == 0 {
		aaMethod = AA_CLOSEDFORM
	}
	return &Checkerboard2DTextureSpectrum{mapping, tex1, tex2, aaMethod}
}

func (tex *Checkerboard2DTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	s, t, dsdx, dtdx, dsdy, dtdy := tex.mapping.Map(dg)
	if tex.aaMethod == AA_NONE {
		// Point sample _Checkerboard2DTexture_
		if (Floor2Int(s)+Floor2Int(t))%2 == 0 {
			return tex.tex1.Evaluate(dg)
		}
		return tex.tex2.Evaluate(dg)
	} else {
		// Compute closed-form box-filtered _Checkerboard2DTexture_ value

		// Evaluate single check if filter is entirely inside one of them
		ds := math.Max(math.Abs(dsdx), math.Abs(dsdy))
		dt := math.Max(math.Abs(dtdx), math.Abs(dtdy))
		s0, s1 := s-ds, s+ds
		t0, t1 := t-dt, t+dt
		if Floor2Int(s0) == Floor2Int(s1) && Floor2Int(t0) == Floor2Int(t1) {
			// Point sample _Checkerboard2DTexture_
			if (Floor2Int(s)+Floor2Int(t))%2 == 0 {
				return tex.tex1.Evaluate(dg)
			}
			return tex.tex2.Evaluate(dg)
		}

		// Apply box filter to checkerboard region
		BUMPINT := func(x float64) float64 {
			return (float64(Floor2Int((x)/2)) + 2.0*math.Max((x/2)-float64(Floor2Int(x/2))-0.5, 0.0))
		}
		sint := (BUMPINT(s1) - BUMPINT(s0)) / (2.0 * ds)
		tint := (BUMPINT(t1) - BUMPINT(t0)) / (2.0 * dt)
		area2 := sint + tint - 2.0*sint*tint
		if ds > 1.0 || dt > 1.0 {
			area2 = 0.5
		}
		spec1 := tex.tex1.Evaluate(dg)
		spec2 := tex.tex2.Evaluate(dg)
		return spec1.Scale(1.0 - area2).Add(spec2.Scale(area2))
	}
}

func NewCheckerboard3DTextureFloat(mapping TextureMapping3D, tex1, tex2 TextureFloat) *Checkerboard3DTextureFloat {
	return &Checkerboard3DTextureFloat{mapping, tex1, tex2}
}
func (tex *Checkerboard3DTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	p, _, _ := tex.mapping.Map(dg)
	if (Floor2Int(p.X)+Floor2Int(p.Y)+Floor2Int(p.Z))%2 == 0 {
		return tex.tex1.Evaluate(dg)
	} else {
		return tex.tex2.Evaluate(dg)
	}
}

func NewCheckerboard3DTextureSpectrum(mapping TextureMapping3D, tex1, tex2 TextureSpectrum) *Checkerboard3DTextureSpectrum {
	return &Checkerboard3DTextureSpectrum{mapping, tex1, tex2}
}
func (tex *Checkerboard3DTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	p, _, _ := tex.mapping.Map(dg)
	if (Floor2Int(p.X)+Floor2Int(p.Y)+Floor2Int(p.Z))%2 == 0 {
		return tex.tex1.Evaluate(dg)
	} else {
		return tex.tex2.Evaluate(dg)
	}
}

func NewConstantTextureFloat(value float64) *ConstantTextureFloat {
	return &ConstantTextureFloat{value}
}
func (t *ConstantTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return t.value
}

func NewConstantTextureSpectrum(value Spectrum) *ConstantTextureSpectrum {
	return &ConstantTextureSpectrum{value}
}
func (t *ConstantTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	return &t.value
}

func NewDotsTextureFloat(mapping TextureMapping2D, tex1, tex2 TextureFloat) *DotsTextureFloat {
	return &DotsTextureFloat{mapping, tex1, tex2}
}

func (tex *DotsTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	// Compute cell indices for dots
	s, t, _, _, _, _ := tex.mapping.Map(dg)
	sCell, tCell := float64(Floor2Int(s+0.5)), float64(Floor2Int(t+0.5))

	// Return _insideDot_ result if point is inside dot
	if Noise(sCell+0.5, tCell+0.5, 0.0) > 0 {
		radius := 0.35
		maxShift := 0.5 - radius
		sCenter := sCell + maxShift*Noise(sCell+1.50, tCell+2.8, 0.0)
		tCenter := tCell + maxShift*Noise(sCell+4.5, tCell+9.8, 0.0)
		ds, dt := s-sCenter, t-tCenter
		if ds*ds+dt*dt < radius*radius {
			return tex.insideDot.Evaluate(dg)
		}
	}
	return tex.outsideDot.Evaluate(dg)
}

func NewDotsTextureSpectrum(mapping TextureMapping2D, tex1, tex2 TextureSpectrum) *DotsTextureSpectrum {
	return &DotsTextureSpectrum{mapping, tex1, tex2}
}
func (tex *DotsTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	// Compute cell indices for dots
	s, t, _, _, _, _ := tex.mapping.Map(dg)
	sCell, tCell := float64(Floor2Int(s+0.5)), float64(Floor2Int(t+0.5))

	// Return _insideDot_ result if point is inside dot
	if Noise(sCell+0.5, tCell+0.5, 0.0) > 0 {
		radius := 0.35
		maxShift := 0.5 - radius
		sCenter := sCell + maxShift*Noise(sCell+1.50, tCell+2.8, 0.0)
		tCenter := tCell + maxShift*Noise(sCell+4.5, tCell+9.8, 0.0)
		ds, dt := s-sCenter, t-tCenter
		if ds*ds+dt*dt < radius*radius {
			return tex.insideDot.Evaluate(dg)
		}
	}
	return tex.outsideDot.Evaluate(dg)
}

func NewFBmTextureFloat(oct int, roughness float64, mapping TextureMapping3D) *FBmTextureFloat {
	return &FBmTextureFloat{roughness, oct, mapping}
}
func (tex *FBmTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	return FBm(P, dpdx, dpdy, tex.omega, tex.octaves)
}

func NewFBmTextureSpectrum(oct int, roughness float64, mapping TextureMapping3D) *FBmTextureSpectrum {
	return &FBmTextureSpectrum{roughness, oct, mapping}
}
func (tex *FBmTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	return NewSpectrum1(FBm(P, dpdx, dpdy, tex.omega, tex.octaves))
}

func NewImageTextureFloat(mapping TextureMapping2D, filename string, trilinear bool, maxanisotropy float64, wrap WrapMode, scale, gamma float64) *ImageTextureFloat {
	Unimplemented()
	return nil
}

func (t *ImageTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	Unimplemented()
	return 0.0
}

func NewImageTextureSpectrum(mapping TextureMapping2D, filename string, trilinear bool, maxanisotropy float64, wrap WrapMode, scale, gamma float64) *ImageTextureSpectrum {
	Unimplemented()
	return nil
}

func (t *ImageTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	Unimplemented()
	return NewSpectrum1(0.0)
}

func NewMarbleTextureFloat(oct int, roughness, scale, variation float64, mapping TextureMapping3D) *MarbleTextureFloat {
	return &MarbleTextureFloat{mapping, oct, roughness, scale, variation}
}
func (t *MarbleTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

const (
	MARBLE_NC   = 9
	MARBLE_NSEG = MARBLE_NC - 3
)

var marbleColors [MARBLE_NC][3]float64 = [MARBLE_NC][3]float64{
	{.58, .58, .6}, {.58, .58, .6}, {.58, .58, .6},
	{.5, .5, .5}, {.6, .59, .58}, {.58, .58, .6},
	{.58, .58, .6}, {.2, .2, .33}, {.58, .58, .6}}

func NewMarbleTextureSpectrum(oct int, roughness, scale, variation float64, mapping TextureMapping3D) *MarbleTextureSpectrum {
	return &MarbleTextureSpectrum{mapping, oct, roughness, scale, variation}
}
func (tex *MarbleTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	P = P.Scale(tex.scale)

	marble := P.Y + tex.variation*FBm(P, dpdx.Scale(tex.scale), dpdy.Scale(tex.scale), tex.omega, tex.octaves)
	t := 0.5 + 0.5*math.Sin(marble)

	// Evaluate marble spline at _t_
	first := Floor2Int(t * MARBLE_NSEG)
	t = (t*MARBLE_NSEG - float64(first))
	c0 := NewSpectrum(marbleColors[first])
	c1 := NewSpectrum(marbleColors[first+1])
	c2 := NewSpectrum(marbleColors[first+2])
	c3 := NewSpectrum(marbleColors[first+3])
	// Bezier spline evaluated with de Castilejau's algorithm
	s0 := c0.Scale(1.0 - t).Add(c1.Scale(t))
	s1 := c1.Scale(1.0 - t).Add(c2.Scale(t))
	s2 := c2.Scale(1.0 - t).Add(c3.Scale(t))
	s0 = s0.Scale(1.0 - t).Add(s1.Scale(t))
	s1 = s1.Scale(1.0 - t).Add(s2.Scale(t))
	// Extra scale of 1.5 to increase variation among colors
	return (s0.Scale(1.0 - t).Add(s1.Scale(t))).Scale(1.5)
}

func (tex *MixTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	t1, t2 := tex.tex1.Evaluate(dg), tex.tex2.Evaluate(dg)
	amt := tex.amount.Evaluate(dg)
	return (1.0-amt)*t1 + amt*t2
}

func (tex *MixTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	t1, t2 := tex.tex1.Evaluate(dg), tex.tex2.Evaluate(dg)
	amt := tex.amount.Evaluate(dg)
	return t1.Scale(1.0 - amt).Add(t2.Scale(amt))
}

func (tex *ScaleTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return tex.tex1.Evaluate(dg) * tex.tex2.Evaluate(dg)
}

func (tex *ScaleTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	t1, t2 := tex.tex1.Evaluate(dg), tex.tex2.Evaluate(dg)
	return t1.Mult(t2)
}

func (tex *UVTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (tex *UVTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	s, t, _, _, _, _ := tex.mapping.Map(dg)
	rgb := [3]float64{s - float64(Floor2Int(s)), t - float64(Floor2Int(t)), 0.0}
	return NewSpectrum(rgb)
}

func (tex *WindyTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	windStrength := FBm(P.Scale(0.1), dpdx.Scale(0.1), dpdy.Scale(0.1), 0.5, 3)
	waveHeight := FBm(P, dpdx, dpdy, 0.5, 6)
	return math.Abs(windStrength) * waveHeight
}

func (tex *WindyTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	windStrength := FBm(P.Scale(0.1), dpdx.Scale(0.1), dpdy.Scale(0.1), 0.5, 3)
	waveHeight := FBm(P, dpdx, dpdy, 0.5, 6)
	return NewSpectrum1(math.Abs(windStrength) * waveHeight)
}

func (tex *WrinkledTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	return Turbulence(P, dpdx, dpdy, tex.omega, tex.octaves)
}

func (tex *WrinkledTextureSpectrum) Evaluate(dg *DifferentialGeometry) *Spectrum {
	P, dpdx, dpdy := tex.mapping.Map(dg)
	return NewSpectrum1(Turbulence(P, dpdx, dpdy, tex.omega, tex.octaves))
}

type texInfo struct {
	filename     string
	doTrilinear  bool
	maxAniso     float64
	wrapMode     WrapMode
	scale, gamma float64
}

var (
	textureFloatCache    map[texInfo]*MIPMapFloat
	textureSpectrumCache map[texInfo]*MIPMapSpectrum
)

func getTextureFloat(filename string, doTrilinear bool, maxAniso float64, wrapMode WrapMode, scale, gamma float64) *MIPMapFloat {
	texKey := texInfo{filename, doTrilinear, maxAniso, wrapMode, scale, gamma}

	tex := textureFloatCache[texKey]
	if tex != nil {
		return tex
	}

	texels, width, height := ReadImage(filename)
	if texels != nil {
		texfloat := make([]float64, len(texels), len(texels))
		for i, t := range texels {
			texfloat[i] = convertInFloat(&t, scale, gamma)
		}
		tex = NewMIPMapFloat(width, height, texfloat, doTrilinear, maxAniso, wrapMode)
	} else {
		// Create one-valued _MIPMap_
		texfloat := make([]float64, 1, 1)
		texfloat[0] = math.Pow(scale, gamma)
		tex = NewMIPMapFloat(1, 1, texfloat, doTrilinear, maxAniso, wrapMode)
	}
	textureFloatCache[texKey] = tex
	//PBRT_LOADED_IMAGE_MAP(const_cast<char *>(filename.c_str()), width, height, sizeof(Tmemory), ret);

	return tex
}

func getTextureSpectrum(filename string, doTrilinear bool, maxAniso float64, wrapMode WrapMode, scale, gamma float64) *MIPMapSpectrum {
	texKey := texInfo{filename, doTrilinear, maxAniso, wrapMode, scale, gamma}

	tex := textureSpectrumCache[texKey]
	if tex != nil {
		return tex
	}

	texels, width, height := ReadImage(filename)
	if texels != nil {
		for i, t := range texels {
			texels[i] = *convertInSpectrum(&t, scale, gamma)
		}
		tex = NewMIPMapSpectrum(width, height, texels, doTrilinear, maxAniso, wrapMode)
	} else {
		// Create one-valued _MIPMap_
		texels := make([]Spectrum, 1, 1)
		texels[0] = *NewSpectrum1(math.Pow(scale, gamma))
		tex = NewMIPMapSpectrum(1, 1, texels, doTrilinear, maxAniso, wrapMode)
	}
	textureSpectrumCache[texKey] = tex
	//PBRT_LOADED_IMAGE_MAP(const_cast<char *>(filename.c_str()), width, height, sizeof(Tmemory), ret);

	return tex
}

func convertInSpectrum(from *Spectrum, scale, gamma float64) *Spectrum {
	return PowSpectrum(from.Scale(scale), gamma)
}

func convertInFloat(from *Spectrum, scale, gamma float64) float64 {
	return math.Pow(scale*from.Y(), gamma)
}

func convertOutSpectrum(from *Spectrum) *Spectrum {
	return from
}
func convertOutFloat(from float64) float64 {
	return from
}
