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
    
}

func CreateRandomBSDFSample(rng *RNG) *BSDFSample {
    return &BSDFSample{[2]float64{rng.RandomFloat(), rng.RandomFloat}, rng.RandomFloat()}
}

func CreateBSDFSample(sample *Sample, offsets *BSDFSampleOffsets, num int) *BSDFSample {
    
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
    
}
func (b *BRDFToBTDF) Sample_f(wo, wi *Vector, u1, u2 float64) (*Spectrum, float64) {
    
}
func (b *BRDFToBTDF) Rho(wo, nSamples int, samples []float64) *Spectrum {
    
}
func (b *BRDFToBTDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
    
}
func (b *BRDFToBTDF) Pdf(wi, wo *Vector) float64 {
    
}   

type ScaledBxDF struct { // BxDF
    bxdf BxDF
    s Spectrum
}
func (b *ScaledBxDF) F(wo, wi *Vector) *Spectrum {
    
}
func (b *ScaledBxDF) Sample_f(wo, wi *Vector, u1, u2 float64) (*Spectrum, float64) {
    
}
func (b *ScaledBxDF) Rho(wo, nSamples int, samples []float64) *Spectrum {
    
}
func (b *ScaledBxDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
    
}
func (b *ScaledBxDF) Pdf(wi, wo *Vector) float64 {
    
}   

type Fresnel interface {
    Evaluate(cosi float64) *Spectrum
}

type FresnelConductor struct {
    k, eta Spectrum
}
type FresnelDielectric {
    eta_t, eta_i S float64
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
    
}
func (bsdf *BSDF) Add(bxdf BxDF) {
    
}
func (bsdf *BSDF) NumComponents() int {
    return bsdf.nBxDFs
}
func (bdsf *BSDF) NumComponentsMatching(flags BxDFType) int {
    
}
func (bsdf *BSDF) WorldToLocal(v *Vector) *Vector {
    
}
func (bsdf *BSDF) LocalToWorld(v *Vector) *Vector {
    
}
func (bsdf *BSDF) f(woW, wiW *Vector, flags BxDFType) Spectrum {
    
}
func (bsdf *BSDF) rho(rng *RNG, flags BxDFType, sqrtSamples int) Spectrum {
    
}
func (bsdf *BSDF) rho2(wo *Vector, rng *RNG, flags BxDFType, sqrtSamples int) Spectrum {
    
}

type SpecularReflection struct { // BxDF
    R Spectrum
    fresnel Fresnel
}

type SpecularTransmission struct { // BxDF
    T Spectrum
    etai, etat float64
    fresnel *FresnelDielectric
}

type Lambertian struct { // BxDF
    R Spectrum
}

type OrenNayar struct { // BxDF
    R Spectrum
    A, B float64
}

type MicrofacetDistribution interface {
    D(wh *Vector) float64
    Sample_f(wo, wi *Vector, u1, u2 float64) float64
    Pdf(wo, wi *Vector) float64
}

type Microfacet struct { // BxDF
    R Spectrum
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
    Rd, Rs Spectrum
    distribution MicrofacetDistribution
}

type RegularHalfangleBRDF struct { // BxDF
    brdf []float64
    nThetaH, nThetaD, nPhiD int
}


type BSSRDF struct {
    e float64
    sig_a, sigp_s Spectrum	
}
func CreateBSSRDF(sa, sps Spectrum, et float64) *BSSRDF {
    return &BSSRDF{et, sa, sps}
}

type IrregIsotropicBRDFSample struct {
    p *Point
    v Spectrum  
}

func FrDiel(cosi, cost float64, etai, etat Spectrum) Spectrum {
    
}

func FrCond(cosi float64, n, k Spectrum) Spectrum {
    
}

func BRDFRemap(wo, wi *Vector) *Point {
    
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
func AbsCosTheta(w *Vector) { return math.Abs(w.z) }
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
