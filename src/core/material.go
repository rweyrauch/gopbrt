package core

import (
	"math"
)

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
		//thetaPhiData *KdTreeIrregIsotropBRDFSample
		regularHalfangleData    []float64
		nThetaH, nThetaD, nPhiD uint32
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
		Kr, sigma_a, sigma_prims_s TextureSpectrum
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
    
    R := m.Kr.Evaluate(dgs)
    R = *R.Clamp(0.0, INFINITY)
    T := m.Kt.Evaluate(dgs)
    T = *T.Clamp(0.0, INFINITY)
    if !R.IsBlack() {
		fresnel := &FresnelDielectric{1.0, ior}
        bsdf.Add(&SpecularReflection{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)}, R, fresnel})
    }        
    if !T.IsBlack() {
		fresnel := &FresnelDielectric{1.0, ior}
        bsdf.Add(&SpecularTransmission{BxDFData{BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)}, T, 1.0, ior, fresnel})
    }    
    return bsdf
}
func (m *GlassMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateKdSubsurfaceMaterial(xform *Transform, mp *TextureParams) *KdSubsurfaceMaterial {
	Unimplemented()
	return nil
}
func (m *KdSubsurfaceMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *KdSubsurfaceMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
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
	kd := m.Kd.Evaluate(dgs)
	kd = *kd.Clamp(0.0, INFINITY)
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

func (m *MatteMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMeasuredMaterial(xform *Transform, mp *TextureParams) *MeasuredMaterial {
	Unimplemented()
	return nil
}
func (m *MeasuredMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MeasuredMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMetalMaterial(xform *Transform, mp *TextureParams) *MetalMaterial {
	Unimplemented()
	return nil
}
func (m *MetalMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MetalMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
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
    R := m.Kr.Evaluate(dgs)
    R = *R.Clamp(0.0, INFINITY)
    if !R.IsBlack() {
		fresnel := &FresnelNoOp{}
        bsdf.Add(&SpecularReflection{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)}, R, fresnel})
    }        
    return bsdf
}
func (m *MirrorMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMixMaterial(xform *Transform, mp *TextureParams, m1, m2 Material) *MixMaterial {
	Unimplemented()
	return nil
}
func (m *MixMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MixMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
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

	kd := m.Kd.Evaluate(dgs)
	kd = *kd.Clamp(0.0, INFINITY)
	if !kd.IsBlack() {
		bsdf.Add(NewLambertian(kd))
	}
	ks := m.Ks.Evaluate(dgs)
	ks = *ks.Clamp(0.0, INFINITY)
	if !ks.IsBlack() {
		fresnel := &FresnelDielectric{1.5, 1.0}
		rough := m.roughness.Evaluate(dgs)
		spec := NewMicrofacet(&ks, fresnel, &Blinn{1.0 / rough})
		bsdf.Add(spec)
	}
	return bsdf
}
func (m *PlasticMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
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
    spec := m.Ks.Evaluate(dgs)
    spec = *spec.Clamp(0.0, INFINITY)
    rough := m.roughness.Evaluate(dgs)
    R := m.Kr.Evaluate(dgs)
    R = *R.Clamp(0.0, INFINITY)

    md := &Blinn{1.0 / rough}
    k := *NewSpectrum1(0.0)
    if !spec.IsBlack() {
        frMf := &FresnelConductor{*fresnelApproxEta(&spec), k}
        bsdf.Add(NewMicrofacet(NewSpectrum1(1.0), frMf, md))
    }
    if !R.IsBlack() {
        frSr := &FresnelConductor{*fresnelApproxEta(&R), k}
        bsdf.Add(&SpecularReflection{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)}, *NewSpectrum1(1.0), frSr})
    }
    return bsdf
}
func (m *ShinyMetalMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubstrateMaterial(xform *Transform, mp *TextureParams) *SubstrateMaterial {
	Unimplemented()
	return nil
}
func (m *SubstrateMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *SubstrateMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubsurfaceMaterial(xform *Transform, mp *TextureParams) *SubsurfaceMaterial {
	Unimplemented()
	return nil
}
func (m *SubsurfaceMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *SubsurfaceMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
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

    r := m.reflect.Evaluate(dgs)
    r = *r.Clamp(0.0, INFINITY)
    t := m.transmit.Evaluate(dgs)
    t = *t.Clamp(0.0, INFINITY)
    if r.IsBlack() && t.IsBlack() { return bsdf }

    kd := m.Kd.Evaluate(dgs)
    kd = *kd.Clamp(0.0, INFINITY)
    if !kd.IsBlack() {
        if !r.IsBlack() { bsdf.Add(NewLambertian(*r.Mult(&kd))) }
        if !t.IsBlack() { bsdf.Add(NewBRDFToBTDF(NewLambertian(*t.Mult(&kd)))) }
    }
    ks := m.Ks.Evaluate(dgs)
    ks = *ks.Clamp(0.0, INFINITY)
    if !ks.IsBlack() {
        rough := m.roughness.Evaluate(dgs)
        if !r.IsBlack() {
            fresnel := &FresnelDielectric{ior, 1.0}
            bsdf.Add(NewMicrofacet(r.Mult(&ks), fresnel, &Blinn{1.0 / rough}))
        }
        if !t.IsBlack() {
            fresnel := &FresnelDielectric{ior, 1.0}
            bsdf.Add(NewBRDFToBTDF(NewMicrofacet(t.Mult(&ks), fresnel, &Blinn{1.0 / rough})))
        }
    }
    return bsdf
}
func (m *TranslucentMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateUberMaterial(xform *Transform, mp *TextureParams) *UberMaterial {
	Unimplemented()
	return nil
}
func (m *UberMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *UberMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}
