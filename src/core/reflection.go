package pbrt

import (
    "math"
)

type BxDFType uint32

const (
    BSDF_REFLECTION = 1 << 0
    BSDF_TRANSMISSION = 1 << 1
    BSDF_DIFFUSE = 1 << 2
    BSDF_GLOSSY = 1 << 3
    BSDF_SPECULAR = 1 << 4
    BSDF_ALL_TYPES = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR
    BSDF_ALL_REFLECTION = BSDF_REFLECTION | BSDF_ALL_TYPES
    BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION | BSDF_ALL_TYPES
    BSDF_ALL = BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION
)

type BSDFSample struct {
    uDir [2]float64
    uComponent float64
}

type BSDFSampleOffsets struct {
    nSamples, componentOffset, dirOffset int
}

func CreateBSDFSampleOffsets(count int, sample *Sample) *BSDFSampleOffsets {
    return nil
}

func CreateRandomBSDFSample(rng *RNG) *BSDFSample {
    return &BSDFSample{[2]float64{rng.RandomFloat(), rng.RandomFloat()}, rng.RandomFloat()}
}

func CreateBSDFSample(sample *Sample, offsets *BSDFSampleOffsets, num int) *BSDFSample {
    return nil   
}

type BxDF interface {
    F(wo, wi *Vector) *Spectrum
    Sample_f(wo, wi *Vector, u1, u2 float64) (*Spectrum, float64)
    Rho(wo, nSamples int, samples []float64) *Spectrum
    Rho2(nSamples int, samples1, samples2 []float64) *Spectrum
    Pdf(wi, wo *Vector) float64    
}

type BRDFToBTDF struct {    // BxDF
    brdf BxDF
}
func (b *BRDFToBTDF) F(wo, wi *Vector) *Spectrum {
    return nil
}
func (b *BRDFToBTDF) Sample_f(wo, wi *Vector, u1, u2 float64) (*Spectrum, float64) {
    return nil, 0.0
}
func (b *BRDFToBTDF) Rho(wo, nSamples int, samples []float64) *Spectrum {
    return nil
}
func (b *BRDFToBTDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
    return nil
}
func (b *BRDFToBTDF) Pdf(wi, wo *Vector) float64 {
    return 0.0
}   

type ScaledBxDF struct { // BxDF
    bxdf BxDF
    s Spectrum
}
func (b *ScaledBxDF) F(wo, wi *Vector) *Spectrum {
    return nil
}
func (b *ScaledBxDF) Sample_f(wo, wi *Vector, u1, u2 float64) (*Spectrum, float64) {
    return nil, 0.0
}
func (b *ScaledBxDF) Rho(wo, nSamples int, samples []float64) *Spectrum {
    return nil
}
func (b *ScaledBxDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
    return nil
}
func (b *ScaledBxDF) Pdf(wi, wo *Vector) float64 {
    return 0.0
}   

type Fresnel interface {
    Evaluate(cosi float64) *Spectrum
}

type FresnelConductor struct {
    k, eta Spectrum
}
type FresnelDielectric struct {
    eta_t, eta_i float64
}
type FresnelNoOp struct {
}

type BSDF struct {
    dgShading DifferentialGeometry
    eta float64
	nn, ng Normal
    sn, tn Vector
    nBxDFs int
    bxdfs []BxDF
}
func CreateBSDF(dg *DifferentialGeometry, ngeom *Normal, eta float64) *BSDF {
    return nil
}
func (bsdf *BSDF) Add(bxdf BxDF) {
    
}
func (bsdf *BSDF) NumComponents() int {
    return bsdf.nBxDFs
}
func (bdsf *BSDF) NumComponentsMatching(flags BxDFType) int {
    return 0   
}
func (bsdf *BSDF) WorldToLocal(v *Vector) *Vector {
    return nil
}
func (bsdf *BSDF) LocalToWorld(v *Vector) *Vector {
    return nil
}
func (bsdf *BSDF) f(woW, wiW *Vector, flags BxDFType) *Spectrum {
    return nil
}
func (bsdf *BSDF) rho(rng *RNG, flags BxDFType, sqrtSamples int) *Spectrum {
    return nil
}
func (bsdf *BSDF) rho2(wo *Vector, rng *RNG, flags BxDFType, sqrtSamples int) *Spectrum {
    return nil
}

type SpecularReflection struct { // BxDF
    R *Spectrum
    fresnel Fresnel
}

type SpecularTransmission struct { // BxDF
    T *Spectrum
    etai, etat float64
    fresnel *FresnelDielectric
}

type Lambertian struct { // BxDF
    R *Spectrum
}

type OrenNayar struct { // BxDF
    R *Spectrum
    A, B float64
}

type MicrofacetDistribution interface {
    D(wh *Vector) float64
    Sample_f(wo, wi *Vector, u1, u2 float64) float64
    Pdf(wo, wi *Vector) float64
}

type Microfacet struct { // BxDF
    R *Spectrum
    distribution MicrofacetDistribution
    fresnel Fresnel
}

type Blinn struct { // MicrofacetDistribution
    exponent float64
}

type Anisotropic struct { // MicrofacetDistribution
    ex, ey float64
}

type FresnelBlend struct { // BxDF
    Rd, Rs *Spectrum
    distribution MicrofacetDistribution
}

type RegularHalfangleBRDF struct { // BxDF
    brdf []float64
    nThetaH, nThetaD, nPhiD int
}


type BSSRDF struct {
    e float64
    sig_a, sigp_s *Spectrum	
}
func CreateBSSRDF(sa, sps *Spectrum, et float64) *BSSRDF {
    return &BSSRDF{et, sa, sps}
}

type IrregIsotropicBRDFSample struct {
    p *Point
    v Spectrum  
}

func FrDiel(cosi, cost float64, etai, etat *Spectrum) *Spectrum {
//    var Rparl Spectrum = ((etat * cosi) - (etai * cost)) /
//                     ((etat * cosi) + (etai * cost))
//    var Rperp Spectrum = ((etai * cosi) - (etat * cost)) /
//                     ((etai * cosi) + (etat * cost))
//    return (Rparl*Rparl + Rperp*Rperp) / 2.0  
    return nil
}

func FrCond(cosi float64, eta, k *Spectrum) *Spectrum {
	tmp := eta.Mult(eta).Add(k.Mult(k)).Scale(cosi*cosi) // (eta*eta + k*k) * cosi*cosi
	etaScaled := eta.Scale(2.0 * cosi)
	oneSpec := CreateSpectrum1(1.0)
    cosSpec := CreateSpectrum1(cosi*cosi)
	
	//Rparl2 = (tmp - (2.0 * eta * cosi) + 1) / (tmp + (2.0 * eta * cosi) + 1)
    Rparl2 := tmp.Sub(etaScaled).Add(oneSpec).Divide(tmp.Add(etaScaled).Add(oneSpec)) 
    tmp_f := eta.Mult(eta).Add(k.Mult(k)) //eta*eta + k*k;
    //Rperp2 = (tmp_f - (2.0 * eta * cosi) + cosi*cosi) / (tmp_f + (2.0 * eta * cosi) + cosi*cosi)
    Rperp2 := tmp_f.Sub(etaScaled).Add(cosSpec).Divide(tmp_f.Add(etaScaled).Add(cosSpec))
    return Rparl2.Add(Rperp2).Scale(0.5) 
}

func BRDFRemap(wo, wi *Vector) *Point {
    cosi, coso := CosTheta(wi), CosTheta(wo)
    sini, sino := SinTheta(wi), SinTheta(wo)
    phii, phio := SphericalPhi(wi), SphericalPhi(wo)
    dphi := phii - phio
    if dphi < 0.0 { dphi += 2.0 * math.Pi }
    if dphi > 2.0 * math.Pi { dphi -= 2.0 * math.Pi }
    if dphi > math.Pi { dphi = 2.0 * math.Pi - dphi }
    return &Point{sini * sino, dphi / math.Pi, cosi * coso}    
}

func Fdr(eta float64) float64 {
    if eta >= 1.0 {
        return -1.4399 / (eta*eta) + 0.7099 / eta + 0.6681 + 0.0636 * eta
    } else {
        return -0.4399 + 0.7099 / eta - 0.3319 / (eta * eta) + 0.0636 / (eta*eta*eta)
    }
}

// BSDF Inline Functions
func CosTheta(w *Vector) float64 { return w.z }
func AbsCosTheta(w *Vector) float64 { return math.Abs(w.z) }
func SinTheta2(w *Vector) float64 {
    return math.Max(0.0, 1.0 - CosTheta(w)*CosTheta(w))
}
func SinTheta(w *Vector) float64 {
    return math.Sqrt(SinTheta2(w))
}

func CosPhi(w *Vector) float64 {
    sintheta := SinTheta(w)
    if sintheta == 0.0 { return 1.0 }
    return Clamp(w.x / sintheta, -1.0, 1.0)
}

func SinPhi(w *Vector) float64 {
    sintheta := SinTheta(w)
    if sintheta == 0.0 { return 0.0 }
    return Clamp(w.y / sintheta, -1.0, 1.0)
}

func SameHemisphere(w, wp *Vector) bool {
    return w.z * wp.z > 0.0
}
