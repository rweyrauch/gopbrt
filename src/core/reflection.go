package core

import (
	"math"
)

type BxDFType uint32

const (
	BSDF_REFLECTION       BxDFType = 1 << 0
	BSDF_TRANSMISSION              = 1 << 1
	BSDF_DIFFUSE                   = 1 << 2
	BSDF_GLOSSY                    = 1 << 3
	BSDF_SPECULAR                  = 1 << 4
	BSDF_ALL_TYPES                 = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR
	BSDF_ALL_REFLECTION            = BSDF_REFLECTION | BSDF_ALL_TYPES
	BSDF_ALL_TRANSMISSION          = BSDF_TRANSMISSION | BSDF_ALL_TYPES
	BSDF_ALL                       = BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION
)
const (
	MAX_BxDFS = 8
)

type BSDFSample struct {
	uDir       [2]float64
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
	Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64)
	Rho(wo *Vector, nSamples int, samples []float64) *Spectrum
	Rho2(nSamples int, samples1, samples2 []float64) *Spectrum
	Pdf(wi, wo *Vector) float64
	Type() BxDFType
}
type BxDFData struct {
	bxdftype BxDFType
}

func BxDFSample_f(bxdf BxDF, wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	// Cosine-sample the hemisphere, flipping the direction if necessary
	wi = CosineSampleHemisphere(u1, u2)
	if wo.z < 0.0 {
		wi.z *= -1.0
	}
	pdf = BxDFPdf(wo, wi)
	return wi, bxdf.F(wo, wi), pdf
}

func BxDFrho(bxdf BxDF, w *Vector, nSamples int, samples []float64) *Spectrum {
	r := NewSpectrum1(0.0)
	for i := 0; i < nSamples; i++ {
		// Estimate one term of $\rho_\roman{hd}$
		wi, f, pdf := bxdf.Sample_f(w, samples[2*i], samples[2*i+1])
		if pdf > 0.0 {
			r = r.Add(f.Scale(float32(AbsCosTheta(wi) / pdf)))
		}
	}
	return r.Scale(1.0 / float32(nSamples))
}

func BxDFrho2(bxdf BxDF, nSamples int, samples1, samples2 []float64) *Spectrum {
	r := NewSpectrum1(0.0)
	for i := 0; i < nSamples; i++ {
		wo := UniformSampleHemisphere(samples1[2*i], samples1[2*i+1])
		pdf_o := 1.0 / 2.0 * math.Pi
		wi, f, pdf_i := bxdf.Sample_f(wo, samples2[2*i], samples2[2*i+1])
		if pdf_i > 0.0 {
			r = r.Add(f.Scale(float32(AbsCosTheta(wi) * AbsCosTheta(wo) / (pdf_o * pdf_i))))
		}
	}
	return r.Scale(1.0 / float32(nSamples) * float32(math.Pi))
}

func BxDFPdf(wi, wo *Vector) float64 {
	if SameHemisphere(wo, wi) {
		return AbsCosTheta(wi) / math.Pi
	} else {
		return 0.0
	}
}

func matchesFlags(bxdf BxDF, flags BxDFType) bool {
	return (bxdf.Type() & flags) == bxdf.Type()
}

type BRDFToBTDF struct { // BxDF
	BxDFData
	brdf BxDF
}

func (b *BRDFToBTDF) F(wo, wi *Vector) *Spectrum {
	return b.brdf.F(wo, b.otherHemisphere(wi))
}

func (b *BRDFToBTDF) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	wi, f, pdf = b.brdf.Sample_f(wo, u1, u2)
	wi = b.otherHemisphere(wi)
	return wi, f, pdf
}

func (b *BRDFToBTDF) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return nil
}
func (b *BRDFToBTDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return nil
}
func (b *BRDFToBTDF) Pdf(wi, wo *Vector) float64 {
	return 0.0
}
func (b *BRDFToBTDF) Type() BxDFType { return b.bxdftype }

func (b *BRDFToBTDF) otherHemisphere(w *Vector) *Vector {
	return CreateVector(w.x, w.y, -w.z)
}

type ScaledBxDF struct { // BxDF
	BxDFData
	bxdf BxDF
	s    Spectrum
}

func (b *ScaledBxDF) F(wo, wi *Vector) *Spectrum {
	return b.bxdf.F(wo, wi).Mult(&b.s)
}
func (b *ScaledBxDF) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	wi, f, pdf = b.bxdf.Sample_f(wo, u1, u2)
	return wi, f.Mult(&b.s), pdf
}
func (b *ScaledBxDF) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return nil
}
func (b *ScaledBxDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return nil
}
func (b *ScaledBxDF) Pdf(wi, wo *Vector) float64 {
	return 0.0
}
func (b *ScaledBxDF) Type() BxDFType { return b.bxdftype }

type Fresnel interface {
	Evaluate(cosi float64) *Spectrum
}

type FresnelConductor struct {
	k, eta Spectrum
}

func (fresnel *FresnelConductor) Evaluate(cosi float64) *Spectrum {
	return FrCond(math.Abs(cosi), &fresnel.eta, &fresnel.k)
}

type FresnelDielectric struct {
	eta_t, eta_i float64
}

func (fresnel *FresnelDielectric) Evaluate(cosi float64) *Spectrum {
	// Compute Fresnel reflectance for dielectric
	cosi = Clamp(cosi, -1.0, 1.0)

	// Compute indices of refraction for dielectric
	entering := cosi > 0.0
	ei, et := fresnel.eta_i, fresnel.eta_t
	if !entering {
		ei, et = et, ei
	}
	// Compute _sint_ using Snell's law
	sint := ei / et * math.Sqrt(math.Max(0.0, 1.0-cosi*cosi))
	if sint >= 1.0 {
		// Handle total internal reflection
		return NewSpectrum1(1.0)
	} else {
		cost := math.Sqrt(math.Max(0.0, 1.0-sint*sint))
		return FrDiel(math.Abs(cosi), cost, NewSpectrum1(float32(ei)), NewSpectrum1(float32(et)))
	}
}

type FresnelNoOp struct {
}

type BSDF struct {
	dgShading DifferentialGeometry
	eta       float64
	nn, ng    Normal
	sn, tn    Vector
	nBxDFs    int
	bxdfs     [MAX_BxDFS]BxDF
}

func NewBSDF(dg *DifferentialGeometry, ngeom *Normal, eta float64) *BSDF {
	bsdf := new(BSDF)
	bsdf.dgShading = *dg
	bsdf.eta = eta
	bsdf.ng = *ngeom
	bsdf.nn = *bsdf.dgShading.nn
	bsdf.sn = *NormalizeVector(bsdf.dgShading.dpdu)
	bsdf.tn = *CrossNormalVector(&bsdf.nn, &bsdf.sn)
	bsdf.nBxDFs = 0
	return bsdf
}

func (bsdf *BSDF) Sample_f(wo *Vector, bsdfSample *BSDFSample, flags BxDFType) (f *Spectrum, wi *Vector, pdf float64, sampledType BxDFType) {
	return nil, nil, 0.0, BSDF_REFLECTION
}

func (bsdf *BSDF) Pdf(woW, wiW *Vector,flags BxDFType) float64 {
    if (bsdf.nBxDFs == 0) { return 0.0 }
    //PBRT_STARTED_BSDF_PDF()
    wo, wi := bsdf.WorldToLocal(woW), bsdf.WorldToLocal(wiW)
    var pdf float64
    matchingComps := 0
    for i := 0; i < bsdf.nBxDFs; i++ {
        if matchesFlags(bsdf.bxdfs[i], flags) {
            matchingComps++
            pdf += bsdf.bxdfs[i].Pdf(wo, wi)
        }
	}        
    v := 0.0
    if matchingComps > 0 {
    	v = pdf / float64(matchingComps)
   	}
    //PBRT_FINISHED_BSDF_PDF()
    return v  	
}

func (bsdf *BSDF) Add(bxdf BxDF) {
	Assert(bsdf.nBxDFs < MAX_BxDFS)
	bsdf.bxdfs[bsdf.nBxDFs] = bxdf
	bsdf.nBxDFs++
}

func (bsdf *BSDF) NumComponents() int {
	return bsdf.nBxDFs
}

func (bsdf *BSDF) NumComponentsMatching(flags BxDFType) int {
	num := 0
	for i := 0; i < bsdf.nBxDFs; i++ {
		if matchesFlags(bsdf.bxdfs[i], flags) {
			num++
		}
	}
	return num
}

func (bsdf *BSDF) WorldToLocal(v *Vector) *Vector {
	return CreateVector(DotVector(v, &bsdf.sn), DotVector(v, &bsdf.tn), DotVectorNormal(v, &bsdf.nn))
}

func (bsdf *BSDF) LocalToWorld(v *Vector) *Vector {
	return CreateVector(bsdf.sn.x*v.x+bsdf.tn.x*v.y+bsdf.nn.x*v.z,
		bsdf.sn.y*v.x+bsdf.tn.y*v.y+bsdf.nn.y*v.z,
		bsdf.sn.z*v.x+bsdf.tn.z*v.y+bsdf.nn.z*v.z)
}

func (bsdf *BSDF) f(woW, wiW *Vector, flags BxDFType) *Spectrum {
	//PBRT_STARTED_BSDF_EVAL();
	wi, wo := bsdf.WorldToLocal(wiW), bsdf.WorldToLocal(woW)
	if DotVectorNormal(wiW, &bsdf.ng)*DotVectorNormal(woW, &bsdf.ng) > 0 { // ignore BTDFs
		flags = BxDFType(flags ^ BSDF_TRANSMISSION)
	} else { // ignore BRDFs
		flags = BxDFType(flags ^ BSDF_REFLECTION)
	}
	ff := NewSpectrum1(0.0)
	for i := 0; i < bsdf.nBxDFs; i++ {
		if matchesFlags(bsdf.bxdfs[i], flags) {
			ff = ff.Add(bsdf.bxdfs[i].F(wo, wi))
		}
	}
	//PBRT_FINISHED_BSDF_EVAL();
	return ff
}

func (bsdf *BSDF) rho(rng *RNG, flags BxDFType, sqrtSamples int) *Spectrum {
	nSamples := sqrtSamples * sqrtSamples
	s1 := make([]float64, 2*nSamples, 2*nSamples)
	StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng, true)
	s2 := make([]float64, 2*nSamples, 2*nSamples)
	StratifiedSample2D(s2, sqrtSamples, sqrtSamples, rng, true)

	ret := NewSpectrum1(0.0)
	for i := 0; i < bsdf.nBxDFs; i++ {
		if matchesFlags(bsdf.bxdfs[i], flags) {
			ret = ret.Add(bsdf.bxdfs[i].Rho2(nSamples, s1, s2))
		}
	}
	return ret
}
func (bsdf *BSDF) rho2(wo *Vector, rng *RNG, flags BxDFType, sqrtSamples int) *Spectrum {
	nSamples := sqrtSamples * sqrtSamples
	s1 := make([]float64, 2*nSamples, 2*nSamples)
	StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng, true)
	ret := NewSpectrum1(0.0)
	for i := 0; i < bsdf.nBxDFs; i++ {
		if matchesFlags(bsdf.bxdfs[i], flags) {
			ret = ret.Add(bsdf.bxdfs[i].Rho(wo, nSamples, s1))
		}
	}
	return ret
}

type SpecularReflection struct { // BxDF
	BxDFData
	R       Spectrum
	fresnel Fresnel
}

type SpecularTransmission struct { // BxDF
	BxDFData
	T          Spectrum
	etai, etat float64
	fresnel    *FresnelDielectric
}

type Lambertian struct { // BxDF
	BxDFData
	R Spectrum
}

type OrenNayar struct { // BxDF
	BxDFData
	R    Spectrum
	A, B float64
}

type MicrofacetDistribution interface {
	D(wh *Vector) float64
	Sample_f(wo, wi *Vector, u1, u2 float64) float64
	Pdf(wo, wi *Vector) float64
}

type Microfacet struct { // BxDF
	BxDFData
	R            Spectrum
	distribution MicrofacetDistribution
	fresnel      Fresnel
}

type Blinn struct { // MicrofacetDistribution
	exponent float64
}

type Anisotropic struct { // MicrofacetDistribution
	ex, ey float64
}

type FresnelBlend struct { // BxDF
	BxDFData
	Rd, Rs       Spectrum
	distribution MicrofacetDistribution
}

type RegularHalfangleBRDF struct { // BxDF
	BxDFData
	brdf                    []float64
	nThetaH, nThetaD, nPhiD int
}

type BSSRDF struct {
	e             float64
	sig_a, sigp_s Spectrum
}

func CreateBSSRDF(sa, sps Spectrum, et float64) *BSSRDF {
	return &BSSRDF{et, sa, sps}
}

type IrregIsotropicBRDFSample struct {
	p *Point
	v Spectrum
}

func FrDiel(cosi, cost float64, etai, etat *Spectrum) *Spectrum {
	cosi_ := float32(cosi)
	cost_ := float32(cost)
	Rparl := (etat.Scale(cosi_).Sub(etai.Scale(cost_))).Divide((etat.Scale(cosi_).Add(etai.Scale(cost_))))
	Rperp := (etai.Scale(cosi_).Sub(etat.Scale(cost_))).Divide((etai.Scale(cosi_).Add(etat.Scale(cost_))))
	return ((Rparl.Mult(Rparl)).Add(Rperp.Mult(Rperp))).Scale(1.0 / 2.0)
}

func FrCond(cosi float64, eta, k *Spectrum) *Spectrum {
	tmp := eta.Mult(eta).Add(k.Mult(k)).Scale(float32(cosi * cosi)) // (eta*eta + k*k) * cosi*cosi
	etaScaled := eta.Scale(2.0 * float32(cosi))
	oneSpec := NewSpectrum1(1.0)
	cosSpec := NewSpectrum1(float32(cosi * cosi))

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
	if dphi < 0.0 {
		dphi += 2.0 * math.Pi
	}
	if dphi > 2.0*math.Pi {
		dphi -= 2.0 * math.Pi
	}
	if dphi > math.Pi {
		dphi = 2.0*math.Pi - dphi
	}
	return &Point{sini * sino, dphi / math.Pi, cosi * coso}
}

func Fdr(eta float64) float64 {
	if eta >= 1.0 {
		return -1.4399/(eta*eta) + 0.7099/eta + 0.6681 + 0.0636*eta
	} else {
		return -0.4399 + 0.7099/eta - 0.3319/(eta*eta) + 0.0636/(eta*eta*eta)
	}
}

// BSDF Inline Functions
func CosTheta(w *Vector) float64    { return w.z }
func AbsCosTheta(w *Vector) float64 { return math.Abs(w.z) }
func SinTheta2(w *Vector) float64 {
	return math.Max(0.0, 1.0-CosTheta(w)*CosTheta(w))
}
func SinTheta(w *Vector) float64 {
	return math.Sqrt(SinTheta2(w))
}

func CosPhi(w *Vector) float64 {
	sintheta := SinTheta(w)
	if sintheta == 0.0 {
		return 1.0
	}
	return Clamp(w.x/sintheta, -1.0, 1.0)
}

func SinPhi(w *Vector) float64 {
	sintheta := SinTheta(w)
	if sintheta == 0.0 {
		return 0.0
	}
	return Clamp(w.y/sintheta, -1.0, 1.0)
}

func SameHemisphere(w, wp *Vector) bool {
	return w.z*wp.z > 0.0
}

func NewLambertian(reflectance Spectrum) *Lambertian {
	b := new(Lambertian)
	b.bxdftype = BSDF_REFLECTION | BSDF_DIFFUSE
	b.R = reflectance
	return b
}
func (b *Lambertian) F(wo, wi *Vector) *Spectrum {
	return b.R.Scale(float32(1.0 / math.Pi))
}

func (b *Lambertian) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	return BxDFSample_f(b, wo, u1, u2)
}

func (b *Lambertian) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)
}

func (b *Lambertian) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)
}

func (b *Lambertian) Pdf(wi, wo *Vector) float64 {
	return BxDFPdf(wi, wo)
}
func (b *Lambertian) Type() BxDFType { return b.bxdftype }

func NewOrenNayar(reflectance Spectrum, sig float64) *OrenNayar {
	on := new(OrenNayar)
	on.bxdftype = BSDF_REFLECTION | BSDF_DIFFUSE
	on.R = reflectance
	sigma := Radians(sig)
	sigma2 := sigma * sigma
	on.A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)))
	on.B = 0.45 * sigma2 / (sigma2 + 0.09)
	return on
}

func (b *OrenNayar) F(wo, wi *Vector) *Spectrum {
	sinthetai := SinTheta(wi)
	sinthetao := SinTheta(wo)
	// Compute cosine term of Oren-Nayar model
	maxcos := 0.0
	if sinthetai > 1.0e-4 && sinthetao > 1.0e-4 {
		sinphii, cosphii := SinPhi(wi), CosPhi(wi)
		sinphio, cosphio := SinPhi(wo), CosPhi(wo)
		dcos := cosphii*cosphio + sinphii*sinphio
		maxcos = math.Max(0.0, dcos)
	}

	// Compute sine and tangent terms of Oren-Nayar model
	var sinalpha, tanbeta float64
	if AbsCosTheta(wi) > AbsCosTheta(wo) {
		sinalpha = sinthetao
		tanbeta = sinthetai / AbsCosTheta(wi)
	} else {
		sinalpha = sinthetai
		tanbeta = sinthetao / AbsCosTheta(wo)
	}
	return b.R.Scale(float32((b.A + b.B*maxcos*sinalpha*tanbeta) / math.Pi))
}

func (b *OrenNayar) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	return BxDFSample_f(b, wo, u1, u2)
}

func (b *OrenNayar) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)
}

func (b *OrenNayar) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)
}

func (b *OrenNayar) Pdf(wi, wo *Vector) float64 {
	return BxDFPdf(wi, wo)
}

func (b *OrenNayar) Type() BxDFType { return b.bxdftype }
