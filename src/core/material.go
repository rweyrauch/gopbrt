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

const copperSamples int = 56

var copperWavelengths []float64 = []float64{
	298.7570554, 302.4004341, 306.1337728, 309.960445, 313.8839949, 317.9081487, 322.036826,
	326.2741526, 330.6244747, 335.092373, 339.6826795, 344.4004944, 349.2512056, 354.2405086,
	359.374429, 364.6593471, 370.1020239, 375.7096303, 381.4897785, 387.4505563, 393.6005651,
	399.9489613, 406.5055016, 413.2805933, 420.2853492, 427.5316483, 435.0322035, 442.8006357,
	450.8515564, 459.2006593, 467.8648226, 476.8622231, 486.2124627, 495.936712, 506.0578694,
	516.6007417, 527.5922468, 539.0616435, 551.0407911, 563.5644455, 576.6705953, 590.4008476,
	604.8008683, 619.92089, 635.8162974, 652.5483053, 670.1847459, 688.8009889, 708.4810171,
	729.3186941, 751.4192606, 774.9011125, 799.8979226, 826.5611867, 855.0632966, 885.6012714,
}

var copperN []float64 = []float64{
	1.400313, 1.38, 1.358438, 1.34, 1.329063, 1.325, 1.3325, 1.34, 1.334375, 1.325,
	1.317812, 1.31, 1.300313, 1.29, 1.281563, 1.27, 1.249062, 1.225, 1.2, 1.18, 1.174375, 1.175,
	1.1775, 1.18, 1.178125, 1.175, 1.172812, 1.17, 1.165312, 1.16, 1.155312, 1.15, 1.142812, 1.135,
	1.131562, 1.12, 1.092437, 1.04, 0.950375, 0.826, 0.645875, 0.468, 0.35125, 0.272, 0.230813, 0.214,
	0.20925, 0.213, 0.21625, 0.223, 0.2365, 0.25, 0.254188, 0.26, 0.28, 0.3,
}

var copperK []float64 = []float64{
	1.662125, 1.687, 1.703313, 1.72, 1.744563, 1.77, 1.791625, 1.81, 1.822125, 1.834,
	1.85175, 1.872, 1.89425, 1.916, 1.931688, 1.95, 1.972438, 2.015, 2.121562, 2.21, 2.177188, 2.13,
	2.160063, 2.21, 2.249938, 2.289, 2.326, 2.362, 2.397625, 2.433, 2.469187, 2.504, 2.535875, 2.564,
	2.589625, 2.605, 2.595562, 2.583, 2.5765, 2.599, 2.678062, 2.809, 3.01075, 3.24, 3.458187, 3.67,
	3.863125, 4.05, 4.239563, 4.43, 4.619563, 4.817, 5.034125, 5.26, 5.485625, 5.717,
}

var (
	CopperN *Spectrum
	CopperK *Spectrum
)

func init() {
	CopperN = SpectrumFromSampled(copperWavelengths, copperN)
	CopperK = SpectrumFromSampled(copperWavelengths, copperK)
}

type Material interface {
	GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF
	GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF
}

func Bump(texture TextureFloat, dg, dgs *DifferentialGeometry) *DifferentialGeometry {
	// Compute offset positions and evaluate displacement texture
	dgEval := *dgs

	// Shift _dgEval_ _du_ in the $u$ direction
	du := 0.5 * (math.Abs(dgs.dudx) + math.Abs(dgs.dudy))
	if du == 0.0 {
		du = 0.01
	}
	dgEval.p = dgs.p.Add(dgs.dpdu.Scale(du))
	dgEval.u = dgs.u + du
	dgEval.nn = NormalizeNormal(CreateNormalFromVector(CrossVector(dgs.dpdu, dgs.dpdv)).Add(dgs.dndu.Scale(du)))

	uDisplace := float64(texture.Evaluate(&dgEval))

	// Shift _dgEval_ _dv_ in the $v$ direction
	dv := 0.5 * (math.Abs(dgs.dvdx) + math.Abs(dgs.dvdy))
	if dv == 0.0 {
		dv = 0.01
	}
	dgEval.p = dgs.p.Add(dgs.dpdv.Scale(dv))
	dgEval.u = dgs.u
	dgEval.v = dgs.v + dv
	dgEval.nn = NormalizeNormal(CreateNormalFromVector(CrossVector(dgs.dpdu, dgs.dpdv)).Add(dgs.dndv.Scale(dv)))
	vDisplace := float64(texture.Evaluate(&dgEval))
	displace := float64(texture.Evaluate(dgs))

	// Compute bump-mapped differential geometry
	dgBump := *dgs
	dgBump.dpdu = dgs.dpdu.Add(CreateVectorFromNormal(dgs.nn).Scale((uDisplace - displace) / du)).Add(CreateVectorFromNormal(dgs.dndu).Scale(displace))
	dgBump.dpdv = dgs.dpdv.Add(CreateVectorFromNormal(dgs.nn).Scale((vDisplace - displace) / dv)).Add(CreateVectorFromNormal(dgs.dndv).Scale(displace))
	dgBump.nn = CreateNormalFromVector(NormalizeVector(CrossVector(dgBump.dpdu, dgBump.dpdv)))
	if Xor(dgs.shape.ReverseOrientation(), dgs.shape.TransformSwapsHandedness()) {
		dgBump.nn = dgBump.nn.Negate()
	}
	// Orient shading normal to match geometric normal
	dgBump.nn = FaceforwardNormal(dgBump.nn, dg.nn)

	return &dgBump
}

type (
	GlassMaterial struct {
		Kr, Kt         TextureSpectrum
		index, bumpMap TextureFloat
	}

	KdSubsurfaceMaterial struct {
		Kd, Kr                     TextureSpectrum
		meanfreepath, eta, bumpMap TextureFloat
	}

	MatteMaterial struct {
		Kd             TextureSpectrum
		sigma, bumpMap TextureFloat
	}

	MeasuredMaterial struct {
		thetaPhiData *KdTree	// IrregIsotropBRDFSample
		regularHalfangleData    []float64
		nThetaH, nThetaD, nPhiD int
		bumpMap                 TextureFloat
	}

	MetalMaterial struct {
		eta, k             TextureSpectrum
		roughness, bumpMap TextureFloat
	}

	MirrorMaterial struct {
		Kr      TextureSpectrum
		bumpMap TextureFloat
	}

	MixMaterial struct {
		m1, m2 Material
		scale  TextureSpectrum
	}

	PlasticMaterial struct {
		Kd, Ks             TextureSpectrum
		roughness, bumpMap TextureFloat
	}

	ShinyMetalMaterial struct {
		Ks, Kr             TextureSpectrum
		roughness, bumpMap TextureFloat
	}

	SubstrateMaterial struct {
		Kd, Ks          TextureSpectrum
		nu, nv, bumpMap TextureFloat
	}

	SubsurfaceMaterial struct {
		scale                      float64
		Kr, sigma_a, sigma_prime_s TextureSpectrum
		eta, bumpMap               TextureFloat
	}

	TranslucentMaterial struct {
		Kd, Ks, reflect, transmit TextureSpectrum
		roughness, bumpMap        TextureFloat
	}

	UberMaterial struct {
		Kd, Ks, Kr, Kt, opacity TextureSpectrum
		roughness, eta, bumpMap TextureFloat
	}
)

func CreateGlassMaterial(xform *Transform, mp *TextureParams) *GlassMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	Kt := mp.GetSpectrumTexture("Kt", *NewSpectrum1(1.0))
	index := mp.GetFloatTexture("index", 1.5)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &GlassMaterial{Kr, Kt, index, bumpMap}
}

func (m *GlassMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}
	ior := m.index.Evaluate(dgs)
	bsdf := NewBSDF(dgs, dgGeom.nn, ior)

	R := m.Kr.Evaluate(dgs).Clamp(0.0, INFINITY)
	T := m.Kt.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !R.IsBlack() {
		fresnel := &FresnelDielectric{1.0, ior}
		bsdf.Add(NewSpecularReflection(R, fresnel))
	}
	if !T.IsBlack() {
		fresnel := &FresnelDielectric{1.0, ior}
		bsdf.Add(NewSpecularTransmission(T, 1.0, ior, fresnel))
	}
	return bsdf
}
func (m *GlassMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateKdSubsurfaceMaterial(xform *Transform, mp *TextureParams) *KdSubsurfaceMaterial {
	kd := mp.GetSpectrumTexture("Kd", *NewSpectrumRGB(0.5, 0.5, 0.5))
	mfp := mp.GetFloatTexture("meanfreepath", 1.0)
	ior := mp.GetFloatTexture("index", 1.3)
	kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &KdSubsurfaceMaterial{kd, kr, mfp, ior, bumpMap}
}

func (m *KdSubsurfaceMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}
	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)
	R := m.Kr.Evaluate(dgs).Clamp(0.0, INFINITY)
	e := m.eta.Evaluate(dgs)
	if !R.IsBlack() {
		bsdf.Add(NewSpecularReflection(R, &FresnelDielectric{1.0, e}))
	}
	return bsdf
}

func (m *KdSubsurfaceMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	e := m.eta.Evaluate(dgShading)
	mfp := m.meanfreepath.Evaluate(dgShading)
	kd := m.Kd.Evaluate(dgShading).Clamp(0.0, INFINITY)
	sigma_a, sigma_prime_s := SubsurfaceFromDiffuse(kd, mfp, e)
	return NewBSSRDF(sigma_a, sigma_prime_s, e)
}

func CreateMatteMaterial(xform *Transform, mp *TextureParams) *MatteMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.5))
	sigma := mp.GetFloatTexture("sigma", 0.0)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &MatteMaterial{Kd, sigma, bumpMap}
}

func (m *MatteMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}
	bsdf := NewBSDF(dgs, dgGeom.nn, 1)

	// Evaluate textures for _MatteMaterial_ material and allocate BRDF
	kd := m.Kd.Evaluate(dgs).Clamp(0.0, INFINITY)
	sig := Clamp(m.sigma.Evaluate(dgs), 0.0, 90.0)
	if !kd.IsBlack() {
		if sig == 0 {
			bsdf.Add(NewLambertian(kd))
		} else {
			bsdf.Add(NewOrenNayar(kd, sig))
		}
	}
	return bsdf
}

func (m *MatteMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMetalMaterial(xform *Transform, mp *TextureParams) *MetalMaterial {
	eta := mp.GetSpectrumTexture("eta", *CopperN)
	k := mp.GetSpectrumTexture("k", *CopperK)
	roughness := mp.GetFloatTexture("roughness", 0.01)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &MetalMaterial{eta, k, roughness, bumpMap}
}

func (m *MetalMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1)

	rough := m.roughness.Evaluate(dgs)
	md := &Blinn{1.0 / rough}

	frMf := &FresnelConductor{m.eta.Evaluate(dgs), m.k.Evaluate(dgs)}
	bsdf.Add(NewMicrofacet(NewSpectrum1(1.0), frMf, md))

	return bsdf
}

func (m *MetalMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMirrorMaterial(xform *Transform, mp *TextureParams) *MirrorMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(0.9))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &MirrorMaterial{Kr, bumpMap}
}
func (m *MirrorMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1)
	R := m.Kr.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !R.IsBlack() {
		fresnel := &FresnelNoOp{}
		bsdf.Add(NewSpecularReflection(R, fresnel))
	}
	return bsdf
}
func (m *MirrorMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMixMaterial(xform *Transform, mp *TextureParams, m1, m2 Material) *MixMaterial {
	scale := mp.GetSpectrumTexture("amount", *NewSpectrum1(0.5))
	return &MixMaterial{m1, m2, scale}
}

func (m *MixMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	b1 := m.m1.GetBSDF(dgGeom, dgShading, arena)
	b2 := m.m2.GetBSDF(dgGeom, dgShading, arena)
	s1 := m.scale.Evaluate(dgShading).Clamp(0.0, INFINITY)
	s2 := (NewSpectrum1(1.0).Sub(s1)).Clamp(0.0, INFINITY)
	n1, n2 := b1.NumComponents(), b2.NumComponents()
	for i := 0; i < n1; i++ {
		b1.bxdfs[i] = NewScaledBxDF(b1.bxdfs[i], s1)
	}
	for i := 0; i < n2; i++ {
		b1.Add(NewScaledBxDF(b2.bxdfs[i], s2))
	}
	return b1
}
func (m *MixMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreatePlasticMaterial(xform *Transform, mp *TextureParams) *PlasticMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &PlasticMaterial{Kd, Ks, roughness, bumpMap}
}

func (m *PlasticMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1)

	kd := m.Kd.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !kd.IsBlack() {
		bsdf.Add(NewLambertian(kd))
	}
	ks := m.Ks.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !ks.IsBlack() {
		fresnel := &FresnelDielectric{1.5, 1.0}
		rough := m.roughness.Evaluate(dgs)
		spec := NewMicrofacet(ks, fresnel, &Blinn{1.0 / rough})
		bsdf.Add(spec)
	}
	return bsdf
}
func (m *PlasticMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateShinyMetalMaterial(xform *Transform, mp *TextureParams) *ShinyMetalMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(1.0))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &ShinyMetalMaterial{Kr, Ks, roughness, bumpMap}
}

func fresnelApproxEta(Fr *Spectrum) *Spectrum {
	reflectance := Fr.Clamp(0.0, 0.999)
	return NewSpectrum1(1.0).Add(SqrtSpectrum(reflectance)).Divide(NewSpectrum1(1.0).Sub(SqrtSpectrum(reflectance)))
}

func fresnelApproxK(Fr *Spectrum) *Spectrum {
	reflectance := Fr.Clamp(0.0, 0.999)
	return SqrtSpectrum(reflectance.Divide(NewSpectrum1(1.0).Sub(reflectance))).Scale(2.0)
}

func (m *ShinyMetalMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)
	spec := m.Ks.Evaluate(dgs).Clamp(0.0, INFINITY)
	rough := m.roughness.Evaluate(dgs)
	R := m.Kr.Evaluate(dgs).Clamp(0.0, INFINITY)

	md := &Blinn{1.0 / rough}
	k := NewSpectrum1(0.0)
	if !spec.IsBlack() {
		frMf := &FresnelConductor{fresnelApproxEta(spec), k}
		bsdf.Add(NewMicrofacet(NewSpectrum1(1.0), frMf, md))
	}
	if !R.IsBlack() {
		frSr := &FresnelConductor{fresnelApproxEta(R), k}
		bsdf.Add(NewSpecularReflection(NewSpectrum1(1.0), frSr))
	}
	return bsdf
}
func (m *ShinyMetalMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubstrateMaterial(xform *Transform, mp *TextureParams) *SubstrateMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.5))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.5))
	uroughness := mp.GetFloatTexture("uroughness", 0.1)
	vroughness := mp.GetFloatTexture("vroughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &SubstrateMaterial{Kd, Ks, uroughness, vroughness, bumpMap}
}

func (m *SubstrateMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)
	d := m.Kd.Evaluate(dgs).Clamp(0.0, INFINITY)
	s := m.Ks.Evaluate(dgs).Clamp(0.0, INFINITY)
	u := m.nu.Evaluate(dgs)
	v := m.nv.Evaluate(dgs)

	if !d.IsBlack() || !s.IsBlack() {
		bsdf.Add(&FresnelBlend{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)}, d, s, NewAnisotropic(1.0/u, 1.0/v)})
	}
	return bsdf
}
func (m *SubstrateMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubsurfaceMaterial(xform *Transform, mp *TextureParams) *SubsurfaceMaterial {
	sa := NewSpectrumRGB(0.0011, 0.0024, 0.014)
	sps := NewSpectrumRGB(2.55, 3.21, 3.77)
	name := mp.FindString("name", "")
	if len(name) != 0 {
		found, sap, spsp := GetVolumeScatteringProperties(name)
		if !found {
			Warning("Named material \"%s\" not found.  Using defaults.", name)
		} else {
			sa = sap
			sps = spsp
		}
	}
	scale := mp.FindFloat("scale", 1.0)

	sigma_a := mp.GetSpectrumTexture("sigma_a", *sa)
	sigma_prime_s := mp.GetSpectrumTexture("sigma_prime_s", *sps)
	ior := mp.GetFloatTexture("index", 1.3)
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &SubsurfaceMaterial{scale, Kr, sigma_a, sigma_prime_s, ior, bumpMap}
}

func (m *SubsurfaceMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)
	R := m.Kr.Evaluate(dgs).Clamp(0.0, INFINITY)
	e := m.eta.Evaluate(dgs)
	if !R.IsBlack() {
		bsdf.Add(NewSpecularReflection(R, &FresnelDielectric{1.0, e}))
	}
	return bsdf
}

func (m *SubsurfaceMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	e := m.eta.Evaluate(dgShading)
	sa := m.sigma_a.Evaluate(dgShading)
	sps := m.sigma_prime_s.Evaluate(dgShading)
	return NewBSSRDF(sa.Scale(m.scale), sps.Scale(m.scale), e)
}

func CreateTranslucentMaterial(xform *Transform, mp *TextureParams) *TranslucentMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	reflect := mp.GetSpectrumTexture("reflect", *NewSpectrum1(0.5))
	transmit := mp.GetSpectrumTexture("transmit", *NewSpectrum1(0.5))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &TranslucentMaterial{Kd, Ks, reflect, transmit, roughness, bumpMap}
}

func (m *TranslucentMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	ior := 1.5
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}
	bsdf := NewBSDF(dgs, dgGeom.nn, ior)

	r := m.reflect.Evaluate(dgs).Clamp(0.0, INFINITY)
	t := m.transmit.Evaluate(dgs).Clamp(0.0, INFINITY)
	if r.IsBlack() && t.IsBlack() {
		return bsdf
	}

	kd := m.Kd.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !kd.IsBlack() {
		if !r.IsBlack() {
			bsdf.Add(NewLambertian(r.Mult(kd)))
		}
		if !t.IsBlack() {
			bsdf.Add(NewBRDFToBTDF(NewLambertian(t.Mult(kd))))
		}
	}
	ks := m.Ks.Evaluate(dgs).Clamp(0.0, INFINITY)
	if !ks.IsBlack() {
		rough := m.roughness.Evaluate(dgs)
		if !r.IsBlack() {
			fresnel := &FresnelDielectric{ior, 1.0}
			bsdf.Add(NewMicrofacet(r.Mult(ks), fresnel, &Blinn{1.0 / rough}))
		}
		if !t.IsBlack() {
			fresnel := &FresnelDielectric{ior, 1.0}
			bsdf.Add(NewBRDFToBTDF(NewMicrofacet(t.Mult(ks), fresnel, &Blinn{1.0 / rough})))
		}
	}
	return bsdf
}
func (m *TranslucentMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateUberMaterial(xform *Transform, mp *TextureParams) *UberMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(0.0))
	Kt := mp.GetSpectrumTexture("Kt", *NewSpectrum1(0.0))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	eta := mp.GetFloatTexture("index", 1.5)
	opacity := mp.GetSpectrumTexture("opacity", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &UberMaterial{Kd, Ks, Kr, Kt, opacity, roughness, eta, bumpMap}
}

func (m *UberMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
	// Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
	var dgs *DifferentialGeometry
	if m.bumpMap != nil {
		dgs = Bump(m.bumpMap, dgGeom, dgShading)
	} else {
		dgs = dgShading
	}

	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)

	op := m.opacity.Evaluate(dgs).Clamp(0.0, INFINITY)
	if op.Equal(NewSpectrum1(1.0)) {
		fresnel := &FresnelDielectric{1.0, 1.0}
		tr := NewSpecularTransmission(NewSpectrum1(1.0).Sub(op), 1.0, 1.0, fresnel)
		bsdf.Add(tr)
	}

	kd := op.Mult(m.Kd.Evaluate(dgs)).Clamp(0.0, INFINITY)
	if !kd.IsBlack() {
		diff := NewLambertian(kd)
		bsdf.Add(diff)
	}

	e := m.eta.Evaluate(dgs)
	ks := op.Mult(m.Ks.Evaluate(dgs)).Clamp(0.0, INFINITY)
	if !ks.IsBlack() {
		fresnel := &FresnelDielectric{e, 1.0}
		rough := m.roughness.Evaluate(dgs)
		spec := NewMicrofacet(ks, fresnel, &Blinn{1.0 / rough})
		bsdf.Add(spec)
	}

	kr := op.Mult(m.Kr.Evaluate(dgs)).Clamp(0.0, INFINITY)
	if !kr.IsBlack() {
		fresnel := &FresnelDielectric{e, 1.0}
		bsdf.Add(NewSpecularReflection(kr, fresnel))
	}

	kt := op.Mult(m.Kt.Evaluate(dgs)).Clamp(0.0, INFINITY)
	if !kt.IsBlack() {
		fresnel := &FresnelDielectric{e, 1.0}
		bsdf.Add(NewSpecularTransmission(kt, e, 1.0, fresnel))
	}
	return bsdf
}

func (m *UberMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}
