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
	return nil
}
func (m *GlassMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *GlassMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateKdSubsurfaceMaterial(xform *Transform, mp *TextureParams) *KdSubsurfaceMaterial {
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
    r := kd.Clamp(0.0, 1.0)
    sig := Clamp(float64(m.sigma.Evaluate(dgs)), 0.0, 90.0)
    if !r.IsBlack() {
        if (sig == 0) {
            bsdf.Add(NewLambertian(*r))
        } else {
            bsdf.Add(NewOrenNayar(*r, sig))
		}
    }
    return bsdf
}

func (m *MatteMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMeasuredMaterial(xform *Transform, mp *TextureParams) *MeasuredMaterial {
	return nil
}
func (m *MeasuredMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MeasuredMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMetalMaterial(xform *Transform, mp *TextureParams) *MetalMaterial {
	return nil
}
func (m *MetalMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MetalMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMirrorMaterial(xform *Transform, mp *TextureParams) *MirrorMaterial {
	return nil
}
func (m *MirrorMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MirrorMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateMixMaterial(xform *Transform, mp *TextureParams, m1, m2 Material) *MixMaterial {
	return nil
}
func (m *MixMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *MixMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreatePlasticMaterial(xform *Transform, mp *TextureParams) *PlasticMaterial {
	return nil
}
func (m *PlasticMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *PlasticMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateShinyMetalMaterial(xform *Transform, mp *TextureParams) *ShinyMetalMaterial {
	return nil
}
func (m *ShinyMetalMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *ShinyMetalMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubstrateMaterial(xform *Transform, mp *TextureParams) *SubstrateMaterial {
	return nil
}
func (m *SubstrateMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *SubstrateMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateSubsurfaceMaterial(xform *Transform, mp *TextureParams) *SubsurfaceMaterial {
	return nil
}
func (m *SubsurfaceMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *SubsurfaceMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateTranslucentMaterial(xform *Transform, mp *TextureParams) *TranslucentMaterial {
	return nil
}
func (m *TranslucentMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *TranslucentMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}

func CreateUberMaterial(xform *Transform, mp *TextureParams) *UberMaterial {
	return nil
}
func (m *UberMaterial) GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF {
	return nil
}
func (m *UberMaterial) GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}
