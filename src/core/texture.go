package core

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

	TextureFloat interface {
		Evaluate(dg *DifferentialGeometry) float64
	}

	TextureSpectrum interface {
		Evaluate(dg *DifferentialGeometry) Spectrum
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
	return Noise(p.x, p.y, p.z)
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
	}
	Checkerboard2DTextureSpectrum struct {
		mapping    TextureMapping2D
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
		//mipmap *MIPMap
	}
	ImageTextureSpectrum struct {
		mapping TextureMapping2D
		//mipmap *MIPMap
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

func (t *BilerpTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *BilerpTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateBilerpFloatTexture(tex2world *Transform, tp *TextureParams) *BilerpTextureFloat {
	return nil
}
func CreateBilerpSpectrumTexture(tex2world *Transform, tp *TextureParams) *BilerpTextureSpectrum {
	return nil
}

func (t *Checkerboard2DTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *Checkerboard2DTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateCheckerboardFloatTexture(tex2world *Transform, tp *TextureParams) *Checkerboard2DTextureFloat {
	return nil
}
func CreateCheckerboardSpectrumTexture(tex2world *Transform, tp *TextureParams) *Checkerboard2DTextureSpectrum {
	return nil
}

func (t *ConstantTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return t.value
}

func (t *ConstantTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return t.value
}

func CreateConstantFloatTexture(tex2world *Transform, tp *TextureParams) *ConstantTextureFloat {
	return &ConstantTextureFloat{tp.FindFloat("value", 1.0)}
}
func CreateConstantSpectrumTexture(tex2world *Transform, tp *TextureParams) *ConstantTextureSpectrum {
	return &ConstantTextureSpectrum{tp.FindSpectrum("value", *CreateSpectrum1(1.0))}
}

func (t *DotsTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *DotsTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateDotsFloatTexture(tex2world *Transform, tp *TextureParams) *DotsTextureFloat {
	return nil
}
func CreateDotsSpectrumTexture(tex2world *Transform, tp *TextureParams) *DotsTextureSpectrum {
	return nil
}

func (t *FBmTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *FBmTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateFBmFloatTexture(tex2world *Transform, tp *TextureParams) *FBmTextureFloat {
	return nil
}
func CreateFBmSpectrumTexture(tex2world *Transform, tp *TextureParams) *FBmTextureSpectrum {
	return nil
}

func (t *ImageTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *ImageTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateImageFloatTexture(tex2world *Transform, tp *TextureParams) *ImageTextureFloat {
	return nil
}
func CreateImageSpectrumTexture(tex2world *Transform, tp *TextureParams) *ImageTextureSpectrum {
	return nil
}

func (t *MarbleTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *MarbleTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateMarbleFloatTexture(tex2world *Transform, tp *TextureParams) *MarbleTextureFloat {
	return nil
}
func CreateMarbleSpectrumTexture(tex2world *Transform, tp *TextureParams) *MarbleTextureSpectrum {
	return nil
}

func (t *MixTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *MixTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateMixFloatTexture(tex2world *Transform, tp *TextureParams) *MixTextureFloat {
	return nil
}
func CreateMixSpectrumTexture(tex2world *Transform, tp *TextureParams) *MixTextureSpectrum {
	return nil
}

func (t *ScaleTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *ScaleTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateScaleFloatTexture(tex2world *Transform, tp *TextureParams) *ScaleTextureFloat {
	return nil
}
func CreateScaleSpectrumTexture(tex2world *Transform, tp *TextureParams) *ScaleTextureSpectrum {
	return nil
}

func (t *UVTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *UVTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateUVFloatTexture(tex2world *Transform, tp *TextureParams) *UVTextureFloat {
	return nil
}
func CreateUVSpectrumTexture(tex2world *Transform, tp *TextureParams) *UVTextureSpectrum {
	return nil
}

func (t *WindyTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *WindyTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateWindyFloatTexture(tex2world *Transform, tp *TextureParams) *WindyTextureFloat {
	return nil
}
func CreateWindySpectrumTexture(tex2world *Transform, tp *TextureParams) *WindyTextureSpectrum {
	return nil
}

func (t *WrinkledTextureFloat) Evaluate(dg *DifferentialGeometry) float64 {
	return 0.0
}

func (t *WrinkledTextureSpectrum) Evaluate(dg *DifferentialGeometry) Spectrum {
	return *CreateSpectrum1(0.0)
}

func CreateWrinkledFloatTexture(tex2world *Transform, tp *TextureParams) *WrinkledTextureFloat {
	return nil
}
func CreateWrinkledSpectrumTexture(tex2world *Transform, tp *TextureParams) *WrinkledTextureSpectrum {
	return nil
}
