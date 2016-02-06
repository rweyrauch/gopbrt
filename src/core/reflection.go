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

type (
	Fresnel interface {
		Evaluate(cosi float64) *Spectrum
	}

	FresnelConductor struct {
		k, eta *Spectrum
	}

	FresnelNoOp struct {
	}

	FresnelDielectric struct {
		eta_i, eta_t float64
	}

	BSDFSample struct {
		uDir       [2]float64
		uComponent float64
	}

	BSDFSampleOffsets struct {
		nSamples, componentOffset, dirOffset int
	}

	IrregIsotropicBRDFSample struct {
		p Point
		v Spectrum
	}

    IrregIsoProc struct {
    	v Spectrum
    	sumWeights float64
    	nFound int
	}

	MicrofacetDistribution interface {
		D(wh *Vector) float64
		Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, pdf float64)
		Pdf(wo, wi *Vector) float64
	}

	Blinn struct { // MicrofacetDistribution
		exponent float64
	}

	Anisotropic struct { // MicrofacetDistribution
		ex, ey float64
	}

	BxDF interface {
		F(wo, wi *Vector) *Spectrum
		Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64)
		Rho(wo *Vector, nSamples int, samples []float64) *Spectrum
		Rho2(nSamples int, samples1, samples2 []float64) *Spectrum
		Pdf(wi, wo *Vector) float64
		Type() BxDFType
	}
	BxDFData struct {
		bxdftype BxDFType
	}

	BRDFToBTDF struct { // BxDF
		BxDFData
		brdf BxDF
	}

	ScaledBxDF struct { // BxDF
		BxDFData
		bxdf BxDF
		s    *Spectrum
	}

	SpecularReflection struct { // BxDF
		BxDFData
		R       *Spectrum
		fresnel Fresnel
	}

	SpecularTransmission struct { // BxDF
		BxDFData
		T          *Spectrum
		etai, etat float64
		fresnel    *FresnelDielectric
	}

	Lambertian struct { // BxDF
		BxDFData
		R *Spectrum
	}

	OrenNayar struct { // BxDF
		BxDFData
		R    *Spectrum
		A, B float64
	}

	Microfacet struct { // BxDF
		BxDFData
		R            *Spectrum
		distribution MicrofacetDistribution
		fresnel      Fresnel
	}

	FresnelBlend struct { // BxDF
		BxDFData
		Rd, Rs       *Spectrum
		distribution MicrofacetDistribution
	}

	IrregIsotropicBRDF struct { // BxDF
		BxDFData
		isoBRDFData *KdTree
	}
	
	RegularHalfangleBRDF struct { // BxDF
		BxDFData
		brdf                    []float64
		nThetaH, nThetaD, nPhiD int
	}

	BSDF struct {
		dgShading DifferentialGeometry
		eta       float64
		nn, ng    Normal
		sn, tn    Vector
		nBxDFs    int
		bxdfs     [MAX_BxDFS]BxDF
	}

	BSSRDF struct {
		eta                    float64
		sigma_a, sigma_prime_s *Spectrum
	}
)

func CreateBSDFSampleOffsets(count int, sample *Sample) *BSDFSampleOffsets {
	return &BSDFSampleOffsets{count, sample.Add1D(count), sample.Add2D(count)}
}

func CreateRandomBSDFSample(rng *RNG) *BSDFSample {
	return &BSDFSample{[2]float64{rng.RandomFloat(), rng.RandomFloat()}, rng.RandomFloat()}
}

func CreateBSDFSample(sample *Sample, offsets *BSDFSampleOffsets, n int) *BSDFSample {
	//Assert(n < sample.n2D[offsets.dirOffset])
	//Assert(n < sample.n1D[offsets.componentOffset])
	bsdfSample := new(BSDFSample)
	bsdfSample.uDir[0] = sample.twoD[offsets.dirOffset][2*n]
	bsdfSample.uDir[1] = sample.twoD[offsets.dirOffset][2*n+1]
	bsdfSample.uComponent = sample.oneD[offsets.componentOffset][n]
	//Assert(bsdfSample.uDir[0] >= 0.0 && bsdfSample.uDir[0] < 1.0)
	//Assert(bsdfSample.uDir[1] >= 0.0 && bsdfSample.uDir[1] < 1.0)
	//Assert(bsdfSample.uComponent >= 0.0 && bsdfSample.uComponent < 1.0)
	return bsdfSample
}

func BxDFSample_f(bxdf BxDF, wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	// Cosine-sample the hemisphere, flipping the direction if necessary
	wi = CosineSampleHemisphere(u1, u2)
	if wo.Z < 0.0 {
		wi.Z *= -1.0
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
			r = r.Add(f.Scale(AbsCosTheta(wi) / pdf))
		}
	}
	return r.Scale(1.0 / float64(nSamples))
}

func BxDFrho2(bxdf BxDF, nSamples int, samples1, samples2 []float64) *Spectrum {
	r := NewSpectrum1(0.0)
	for i := 0; i < nSamples; i++ {
		wo := UniformSampleHemisphere(samples1[2*i], samples1[2*i+1])
		pdf_o := 1.0 / 2.0 * math.Pi
		wi, f, pdf_i := bxdf.Sample_f(wo, samples2[2*i], samples2[2*i+1])
		if pdf_i > 0.0 {
			r = r.Add(f.Scale(AbsCosTheta(wi) * AbsCosTheta(wo) / (pdf_o * pdf_i)))
		}
	}
	return r.Scale(1.0 / float64(nSamples) * math.Pi)
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

func NewBRDFToBTDF(brdf BxDF) *BRDFToBTDF {
	return &BRDFToBTDF{BxDFData{brdf.Type() ^ (BSDF_REFLECTION | BSDF_TRANSMISSION)}, brdf}
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
	return b.brdf.Rho(b.otherHemisphere(wo), nSamples, samples)
}
func (b *BRDFToBTDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return b.brdf.Rho2(nSamples, samples1, samples2)
}
func (b *BRDFToBTDF) Pdf(wi, wo *Vector) float64 {
	return b.brdf.Pdf(wo, b.otherHemisphere(wi))
}
func (b *BRDFToBTDF) Type() BxDFType { return b.bxdftype }

func (b *BRDFToBTDF) otherHemisphere(w *Vector) *Vector {
	return CreateVector(w.X, w.Y, -w.Z)
}

func NewScaledBxDF(bxdf BxDF, s *Spectrum) *ScaledBxDF {
	return &ScaledBxDF{BxDFData{bxdf.Type()}, bxdf, s}
}

func (b *ScaledBxDF) F(wo, wi *Vector) *Spectrum {
	return b.bxdf.F(wo, wi).Mult(b.s)
}
func (b *ScaledBxDF) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	wi, f, pdf = b.bxdf.Sample_f(wo, u1, u2)
	return wi, f.Mult(b.s), pdf
}
func (b *ScaledBxDF) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return b.bxdf.Rho(wo, nSamples, samples).Mult(b.s)
}
func (b *ScaledBxDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return b.bxdf.Rho2(nSamples, samples1, samples2).Mult(b.s)
}
func (b *ScaledBxDF) Pdf(wi, wo *Vector) float64 {
	return BxDFPdf(wi, wo)
}
func (b *ScaledBxDF) Type() BxDFType { return b.bxdftype }

func (fresnel *FresnelConductor) Evaluate(cosi float64) *Spectrum {
	return FrCond(math.Abs(cosi), fresnel.eta, fresnel.k)
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
		return FrDiel(math.Abs(cosi), cost, NewSpectrum1(ei), NewSpectrum1(et))
	}
}

func (fresnel *FresnelNoOp) Evaluate(cosi float64) *Spectrum {
	return NewSpectrum1(1.0)
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

func (bsdf *BSDF) Sample_f(woW *Vector, bsdfSample *BSDFSample, flags BxDFType) (f *Spectrum, wi *Vector, pdf float64, sampledType BxDFType) {
	//PBRT_STARTED_BSDF_SAMPLE();
	// Choose which _BxDF_ to sample
	matchingComps := bsdf.NumComponentsMatching(flags)
	if matchingComps == 0 {
		pdf = 0.0
		sampledType = BxDFType(0)
		//PBRT_FINISHED_BSDF_SAMPLE()
		wi = nil
		return NewSpectrum1(0.0), wi, pdf, sampledType
	}
	which := Mini(Floor2Int(bsdfSample.uComponent*float64(matchingComps)), matchingComps-1)
	var bxdf BxDF = nil
	count := which
	for i := 0; i < bsdf.nBxDFs; i++ {
		if matchesFlags(bsdf.bxdfs[i], flags) {
			if count == 0 {
				bxdf = bsdf.bxdfs[i]
				break
			}
			count--
		}
	}
	Assert(bxdf != nil)

	// Sample chosen _BxDF_
	wo := bsdf.WorldToLocal(woW)
	pdf = 0.0
	wi, f, pdf = bxdf.Sample_f(wo, bsdfSample.uDir[0], bsdfSample.uDir[1])
	if pdf == 0.0 {
		sampledType = BxDFType(0)
		//PBRT_FINISHED_BSDF_SAMPLE()
		return f, wi, pdf, sampledType
	}
	sampledType = bxdf.Type()
	wiW := bsdf.LocalToWorld(wi)

	// Compute overall PDF with all matching _BxDF_s
	if !((bxdf.Type()&BSDF_SPECULAR) == 0 && matchingComps > 1) {
		for i := 0; i < bsdf.nBxDFs; i++ {
			if bsdf.bxdfs[i] != bxdf && matchesFlags(bsdf.bxdfs[i], flags) {
				pdf += bsdf.bxdfs[i].Pdf(wo, wi)
			}
		}
	}
	if matchingComps > 1 {
		pdf /= float64(matchingComps)
	}

	// Compute value of BSDF for sampled direction
	if (bxdf.Type() & BSDF_SPECULAR) != 0 {
		f = NewSpectrum1(0.0)
		if DotVectorNormal(wiW, &bsdf.ng)*DotVectorNormal(woW, &bsdf.ng) > 0 { // ignore BTDFs
			flags = BxDFType(flags ^ BSDF_TRANSMISSION)
		} else { // ignore BRDFs
			flags = BxDFType(flags ^ BSDF_REFLECTION)
		}
		for i := 0; i < bsdf.nBxDFs; i++ {
			if matchesFlags(bsdf.bxdfs[i], flags) {
				f = f.Add(bsdf.bxdfs[i].F(wo, wi))
			}
		}
	}
	//PBRT_FINISHED_BSDF_SAMPLE()
	return f, wiW, pdf, sampledType
}

func (bsdf *BSDF) Pdf(woW, wiW *Vector, flags BxDFType) float64 {
	if bsdf.nBxDFs == 0 {
		return 0.0
	}
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
	Assert(bxdf.Type() != 0)
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
	return CreateVector(bsdf.sn.X*v.X+bsdf.tn.X*v.Y+bsdf.nn.X*v.Z,
		bsdf.sn.Y*v.X+bsdf.tn.Y*v.Y+bsdf.nn.Y*v.Z,
		bsdf.sn.Z*v.X+bsdf.tn.Z*v.Y+bsdf.nn.Z*v.Z)
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

func NewBSSRDF(sa, sps *Spectrum, et float64) *BSSRDF {
	return &BSSRDF{et, sa, sps}
}

func FrDiel(cosi, cost float64, etai, etat *Spectrum) *Spectrum {
	Rparl := (etat.Scale(cosi).Sub(etai.Scale(cost))).Divide((etat.Scale(cosi).Add(etai.Scale(cost))))
	Rperp := (etai.Scale(cosi).Sub(etat.Scale(cost))).Divide((etai.Scale(cosi).Add(etat.Scale(cost))))
	return ((Rparl.Mult(Rparl)).Add(Rperp.Mult(Rperp))).Scale(1.0 / 2.0)
}

func FrCond(cosi float64, eta, k *Spectrum) *Spectrum {
	tmp := eta.Mult(eta).Add(k.Mult(k)).Scale(cosi * cosi) // (eta*eta + k*k) * cosi*cosi
	etaScaled := eta.Scale(2.0 * cosi)
	oneSpec := NewSpectrum1(1.0)
	cosSpec := NewSpectrum1(cosi * cosi)

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
func CosTheta(w *Vector) float64    { return w.Z }
func AbsCosTheta(w *Vector) float64 { return math.Abs(w.Z) }
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
	return Clamp(w.X/sintheta, -1.0, 1.0)
}

func SinPhi(w *Vector) float64 {
	sintheta := SinTheta(w)
	if sintheta == 0.0 {
		return 0.0
	}
	return Clamp(w.Y/sintheta, -1.0, 1.0)
}

func SameHemisphere(w, wp *Vector) bool {
	return w.Z*wp.Z > 0.0
}

func NewSpecularReflection(r *Spectrum, fresnel Fresnel) *SpecularReflection {
	return &SpecularReflection{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)}, r, fresnel}
}

func (b *SpecularReflection) F(wo, wi *Vector) *Spectrum {
	return NewSpectrum1(0.0)
}
func (b *SpecularReflection) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	// Compute perfect specular reflection direction
	wi = CreateVector(-wo.X, -wo.Y, wo.Z)
	pdf = 1.0
	f = (b.fresnel.Evaluate(CosTheta(wo)).Mult(b.R)).InvScale(AbsCosTheta(wi))
	return wi, f, pdf
}
func (b *SpecularReflection) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)
}
func (b *SpecularReflection) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)
}
func (b *SpecularReflection) Pdf(wi, wo *Vector) float64 {
	return 0.0
}
func (b *SpecularReflection) Type() BxDFType { return b.bxdftype }

func NewSpecularTransmission(t *Spectrum, etai, etat float64, fresnel *FresnelDielectric) *SpecularTransmission {
	return &SpecularTransmission{BxDFData{BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)}, t, etai, etat, fresnel}
}

func (b *SpecularTransmission) F(wo, wi *Vector) *Spectrum {
	return NewSpectrum1(0.0)
}

func (b *SpecularTransmission) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	// Figure out which $\eta$ is incident and which is transmitted
	entering := CosTheta(wo) > 0.0
	ei, et := b.etai, b.etat
	if !entering {
		ei, et = et, ei
	}
	// Compute transmitted ray direction
	sini2 := SinTheta2(wo)
	eta := ei / et
	sint2 := eta * eta * sini2

	// Handle total internal reflection for transmission
	if sint2 >= 1.0 {
		return nil, NewSpectrum1(0.0), 0.0
	}
	cost := math.Sqrt(math.Max(0.0, 1.0-sint2))
	if entering {
		cost = -cost
	}
	sintOverSini := eta
	wi = CreateVector(sintOverSini*-wo.X, sintOverSini*-wo.Y, cost)
	pdf = 1.0
	F := b.fresnel.Evaluate(CosTheta(wo))
	f = (NewSpectrum1(1.0).Sub(F)).Mult(b.T.InvScale(AbsCosTheta(wi)))
	return wi, f, pdf
}

func (b *SpecularTransmission) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)
}
func (b *SpecularTransmission) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)
}
func (b *SpecularTransmission) Pdf(wi, wo *Vector) float64 {
	return 0.0
}
func (b *SpecularTransmission) Type() BxDFType { return b.bxdftype }

func NewLambertian(reflectance *Spectrum) *Lambertian {
	b := new(Lambertian)
	b.bxdftype = BSDF_REFLECTION | BSDF_DIFFUSE
	b.R = reflectance
	return b
}
func (b *Lambertian) F(wo, wi *Vector) *Spectrum {
	return b.R.Scale(1.0 / math.Pi)
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

func NewOrenNayar(reflectance *Spectrum, sig float64) *OrenNayar {
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
	return b.R.Scale((b.A + b.B*maxcos*sinalpha*tanbeta) / math.Pi)
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

func NewMicrofacet(reflectance *Spectrum, f Fresnel, d MicrofacetDistribution) *Microfacet {
	return &Microfacet{BxDFData{BSDF_REFLECTION | BSDF_GLOSSY}, reflectance, d, f}
}
func (b *Microfacet) F(wo, wi *Vector) *Spectrum {
	cosThetaO := AbsCosTheta(wo)
	cosThetaI := AbsCosTheta(wi)
	if cosThetaI == 0.0 || cosThetaO == 0.0 {
		return NewSpectrum1(0.0)
	}
	wh := wi.Add(wo)
	if wh.X == 0.0 && wh.Y == 0.0 && wh.Z == 0.0 {
		return NewSpectrum1(0.0)
	}
	wh = NormalizeVector(wh)
	cosThetaH := DotVector(wi, wh)
	fr := b.fresnel.Evaluate(cosThetaH)
	return b.R.Scale(b.distribution.D(wh) * b.G(wo, wi, wh)).Mult(fr).Scale((1.0 / (4.0 * cosThetaI * cosThetaO)))
}

func (b *Microfacet) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	wi, pdf = b.distribution.Sample_f(wo, u1, u2)
	if !SameHemisphere(wo, wi) {
		return wi, NewSpectrum1(0.0), pdf
	}
	return wi, b.F(wo, wi), pdf
}

func (b *Microfacet) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)
}

func (b *Microfacet) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)
}

func (b *Microfacet) Pdf(wi, wo *Vector) float64 {
	if !SameHemisphere(wo, wi) {
		return 0.0
	}
	return b.distribution.Pdf(wo, wi)
}

func (b *Microfacet) Type() BxDFType { return b.bxdftype }

func (b *Microfacet) G(wo, wi, wh *Vector) float64 {
	NdotWh := AbsCosTheta(wh)
	NdotWo := AbsCosTheta(wo)
	NdotWi := AbsCosTheta(wi)
	WOdotWh := AbsDotVector(wo, wh)
	return math.Min(1.0, math.Min((2.0*NdotWh*NdotWo/WOdotWh), (2.0*NdotWh*NdotWi/WOdotWh)))
}

func isoBRDFProc(p *Point, nodeData NodeData, dist2 float64, maxDistSquared *float64) {
	Unimplemented()
}

func (b *IrregIsotropicBRDF) F(wo, wi *Vector) *Spectrum {
	Unimplemented()
	return NewSpectrum1(0.0)
/*	
    m := BRDFRemap(wo, wi)
    lastMaxDist2 := 0.001
    for {
        // Try to find enough BRDF samples around _m_ within search radius
        var proc IrregIsoProc
        maxDist2 := lastMaxDist2
        b.isoBRDFData.Lookup(m, isoBRDFProc, &maxDist2)
        if proc.nFound > 2 || lastMaxDist2 > 1.5 {
            return proc.v.Clamp(0.0, Infinity).InvScale(proc.sumWeights)
        }    
        lastMaxDist2 *= 2.0
    }
*/    	
}

func (b *IrregIsotropicBRDF) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	return BxDFSample_f(b, wo, u1, u2)	
}
func (b *IrregIsotropicBRDF) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)	
}
func (b *IrregIsotropicBRDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)	
}
func (b *IrregIsotropicBRDF) Pdf(wi, wo *Vector) float64 {
	return BxDFPdf(wi, wo)	
}

func (b *IrregIsotropicBRDF) Type() BxDFType { return b.bxdftype }


func (b *RegularHalfangleBRDF) F(WO, WI *Vector) *Spectrum {
    // Compute $\wh$ and transform $\wi$ to halfangle coordinate system
    wo, wi := WO, WI
    wh := wo.Add(wi)
    if wh.Z < 0.0 {
        wo = wo.Negate()
        wi = wi.Negate()
        wh = wh.Negate()
    }
    if wh.X == 0.0 && wh.Y == 0.0 && wh.Z == 0.0 { return NewSpectrum1(0.0) }
    wh = NormalizeVector(wh)
    whTheta := SphericalTheta(wh)
    whCosPhi, whSinPhi := CosPhi(wh), SinPhi(wh)
    whCosTheta, whSinTheta := CosTheta(wh), SinTheta(wh)
    whx := CreateVector(whCosPhi * whCosTheta, whSinPhi * whCosTheta, -whSinTheta)
    why := CreateVector(-whSinPhi, whCosPhi, 0)
    wd := CreateVector(DotVector(wi, whx), DotVector(wi, why), DotVector(wi, wh))

    // Compute _index_ into measured BRDF tables
    wdTheta, wdPhi := SphericalTheta(wd), SphericalPhi(wd)
    if wdPhi > math.Pi { wdPhi -= math.Pi }

    // Compute indices _whThetaIndex_, _wdThetaIndex_, _wdPhiIndex_
 	REMAP := func(V, MAX float64, COUNT int) int {
        return Clampi(int(V / MAX * float64(COUNT)), 0, COUNT-1)
    }
    whThetaIndex := REMAP(math.Sqrt(math.Max(0.0, whTheta / (math.Pi / 2.0))), 1.0, b.nThetaH)
    wdThetaIndex := REMAP(wdTheta, math.Pi / 2.0, b.nThetaD)
    wdPhiIndex := REMAP(wdPhi, math.Pi, b.nPhiD)

    index := wdPhiIndex + b.nPhiD * (wdThetaIndex + whThetaIndex * b.nThetaD)
    return NewSpectrumRGB(b.brdf[3*index], b.brdf[3*index+1], b.brdf[3*index+2])	
}
func (b *RegularHalfangleBRDF) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	return BxDFSample_f(b, wo, u1, u2)	
}
func (b *RegularHalfangleBRDF) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(b, wo, nSamples, samples)	
}
func (b *RegularHalfangleBRDF) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(b, nSamples, samples1, samples2)	
}
func (b *RegularHalfangleBRDF) Pdf(wi, wo *Vector) float64 {
	return BxDFPdf(wi, wo)	
}

func (b *RegularHalfangleBRDF) Type() BxDFType { return b.bxdftype }


func NewFresnelBlend(d, s *Spectrum, dist MicrofacetDistribution) *FresnelBlend {
	fb := &FresnelBlend{BxDFData{BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)}, d, s, dist}
	return fb
}

func (blend *FresnelBlend) schlickFresnel(costheta float64) *Spectrum {
	return blend.Rs.Add((NewSpectrum1(1.0).Sub(blend.Rs)).Scale(math.Pow(1.0-costheta, 5.0)))
}

func (blend *FresnelBlend) F(wo, wi *Vector) *Spectrum {

	diffuse := blend.Rd.Scale(28.0 / (23.0 * math.Pi)).Mult((NewSpectrum1(1.0).Sub(blend.Rs)).Scale((1.0 - math.Pow(1.0-0.5*AbsCosTheta(wi), 5)) * (1.0 - math.Pow(1.0-0.5*AbsCosTheta(wo), 5))))
	wh := wi.Add(wo)
	if wh.X == 0.0 && wh.Y == 0.0 && wh.Z == 0.0 {
		return NewSpectrum1(0.0)
	}
	wh = NormalizeVector(wh)
	specular := blend.schlickFresnel(DotVector(wi, wh)).Scale(blend.distribution.D(wh) / (4.0 * AbsDotVector(wi, wh) * math.Max(AbsCosTheta(wi), AbsCosTheta(wo))))
	return diffuse.Add(specular)
}

func (blend *FresnelBlend) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, f *Spectrum, pdf float64) {
	if u1 < 0.5 {
		u1 = 2.0 * u1
		// Cosine-sample the hemisphere, flipping the direction if necessary
		wi = CosineSampleHemisphere(u1, u2)
		if wo.Z < 0.0 {
			wi.Z *= -1.0
		}
	} else {
		u1 = 2.0 * (u1 - 0.5)
		wi, pdf = blend.distribution.Sample_f(wo, u1, u2)
		if !SameHemisphere(wo, wi) {
			return wi, NewSpectrum1(0.0), pdf
		}
	}
	pdf = blend.Pdf(wo, wi)
	return wi, blend.F(wo, wi), pdf
}

func (blend *FresnelBlend) Rho(wo *Vector, nSamples int, samples []float64) *Spectrum {
	return BxDFrho(blend, wo, nSamples, samples)
}
func (blend *FresnelBlend) Rho2(nSamples int, samples1, samples2 []float64) *Spectrum {
	return BxDFrho2(blend, nSamples, samples1, samples2)
}
func (blend *FresnelBlend) Pdf(wi, wo *Vector) float64 {
	if !SameHemisphere(wo, wi) {
		return 0.0
	}
	return 0.5 * (AbsCosTheta(wi)*INV_PI + blend.distribution.Pdf(wo, wi))
}
func (f *FresnelBlend) Type() BxDFType {
	return f.bxdftype
}

func (b *Blinn) D(wh *Vector) float64 {
	costhetah := AbsCosTheta(wh)
	return (b.exponent + 2) / (2.0 * math.Pi) * math.Pow(costhetah, b.exponent)
}

func (b *Blinn) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, pdf float64) {
	// Compute sampled half-angle vector $\wh$ for Blinn distribution
	costheta := math.Pow(u1, 1.0/(b.exponent+1))
	sintheta := math.Sqrt(math.Max(0.0, 1.0-costheta*costheta))
	phi := u2 * 2.0 * math.Pi
	wh := SphericalDirection(sintheta, costheta, phi)
	if !SameHemisphere(wo, wh) {
		wh = wh.Negate()
	}

	// Compute incident direction by reflecting about $\wh$
	wi = wo.Negate().Add(wh.Scale(2.0 * DotVector(wo, wh)))

	// Compute PDF for $\wi$ from Blinn distribution
	blinn_pdf := ((b.exponent + 1.0) * math.Pow(costheta, b.exponent)) / (2.0 * math.Pi * 4.0 * DotVector(wo, wh))
	if DotVector(wo, wh) <= 0.0 {
		blinn_pdf = 0.0
	}
	pdf = blinn_pdf
	return wi, pdf
}

func (b *Blinn) Pdf(wo, wi *Vector) float64 {
	wh := NormalizeVector(wo.Add(wi))
	costheta := AbsCosTheta(wh)
	// Compute PDF for $\wi$ from Blinn distribution
	blinn_pdf := ((b.exponent + 1.0) * math.Pow(costheta, b.exponent)) / (2.0 * math.Pi * 4.0 * DotVector(wo, wh))
	if DotVector(wo, wh) <= 0.0 {
		blinn_pdf = 0.0
	}
	return blinn_pdf
}

func NewAnisotropic(x, y float64) *Anisotropic {
	if x > 10000.0 || math.IsNaN(x) {
		x = 10000.0
	}
	if y > 10000.0 || math.IsNaN(y) {
		y = 10000.0
	}
	return &Anisotropic{x, y}
}

func (a *Anisotropic) D(wh *Vector) float64 {
	costhetah := AbsCosTheta(wh)
	d := 1.0 - costhetah*costhetah
	if d == 0.0 {
		return 0.0
	}
	e := (a.ex*wh.X*wh.X + a.ey*wh.Y*wh.Y) / d
	return math.Sqrt((a.ex+2.0)*(a.ey+2.0)) * (1.0 / (2.0 * math.Pi)) * math.Pow(costhetah, e)
}

func (a *Anisotropic) Sample_f(wo *Vector, u1, u2 float64) (wi *Vector, pdf float64) {
	// Sample from first quadrant and remap to hemisphere to sample $\wh$
	var phi, costheta float64
	if u1 < 0.25 {
		phi, costheta = a.sampleFirstQuadrant(4.0*u1, u2)
	} else if u1 < 0.5 {
		u1 = 4.0 * (0.5 - u1)
		phi, costheta = a.sampleFirstQuadrant(u1, u2)
		phi = math.Pi - phi
	} else if u1 < 0.75 {
		u1 = 4.0 * (u1 - 0.5)
		phi, costheta = a.sampleFirstQuadrant(u1, u2)
		phi += math.Pi
	} else {
		u1 = 4.0 * (1.0 - u1)
		phi, costheta = a.sampleFirstQuadrant(u1, u2)
		phi = 2.0*math.Pi - phi
	}
	sintheta := math.Sqrt(math.Max(0.0, 1.0-costheta*costheta))
	wh := SphericalDirection(sintheta, costheta, phi)
	if !SameHemisphere(wo, wh) {
		wh = wh.Negate()
	}

	// Compute incident direction by reflecting about $\wh$
	dotwowh := DotVector(wo, wh)
	wi = wh.Scale(2.0 * dotwowh).Sub(wo)

	// Compute PDF for $\wi$ from anisotropic distribution
	costhetah := AbsCosTheta(wh)
	ds := 1.0 - costhetah*costhetah
	anisotropic_pdf := 0.0
	if ds > 0.0 && DotVector(wo, wh) > 0.0 {
		e := (a.ex*wh.X*wh.X + a.ey*wh.Y*wh.Y) / ds
		d := math.Sqrt((a.ex+1.0)*(a.ey+1.0)) * INV_TWOPI * math.Pow(costhetah, e)
		anisotropic_pdf = d / (4.0 * DotVector(wo, wh))
	}
	pdf = anisotropic_pdf
	return wi, pdf
}

func (a *Anisotropic) Pdf(wo, wi *Vector) float64 {
	wh := NormalizeVector(wo.Add(wi))
	// Compute PDF for $\wi$ from anisotropic distribution
	costhetah := AbsCosTheta(wh)
	ds := 1.0 - costhetah*costhetah
	anisotropic_pdf := 0.0
	if ds > 0.0 && DotVector(wo, wh) > 0.0 {
		e := (a.ex*wh.X*wh.X + a.ey*wh.Y*wh.Y) / ds
		d := math.Sqrt((a.ex+1.0)*(a.ey+1.0)) * INV_TWOPI * math.Pow(costhetah, e)
		anisotropic_pdf = d / (4.0 * DotVector(wo, wh))
	}
	return anisotropic_pdf
}

func (a *Anisotropic) sampleFirstQuadrant(u1, u2 float64) (phi, costheta float64) {
	if a.ex == a.ey {
		phi = math.Pi * u1 * 0.5
	} else {
		phi = math.Atan(math.Sqrt((a.ex+1.0)/(a.ey+1.0)) * math.Tan(math.Pi*u1*0.5))
	}
	cosphi, sinphi := math.Cos(phi), math.Sin(phi)
	costheta = math.Pow(u2, 1.0/(a.ex*cosphi*cosphi+a.ey*sinphi*sinphi+1))
	return phi, costheta
}
