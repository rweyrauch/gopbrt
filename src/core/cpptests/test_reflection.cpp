#include "api.h"
#include "pbrt.h"
#include "reflection.h"
#include "montecarlo.h"
#include <stdio.h>

struct FrDielInputs {
	float cosi, cost;
	float etai[3];
	float etat[3];
};

const FrDielInputs frDielInputVector[] =
{
	{ 1.0, 1.0, { 0.5, 0.5, 0.5 }, { 0.25, 0.25, 0.25 } },
	{ 0.5, 0.717, { 1.0, 0.5, 0.25 }, { 0.1, 0.2, 0.3 } },
	{ -0.2, 0.0, { 0.5, 0.85, 0.65 }, { 0.75, 0.05, 0.5 } },
};

struct FrCondInputs {
	float cosi;
	float n[3];
	float k[3];
};

const FrCondInputs frCondInputVector[] =
{
	{ 1.0, { 0.5, 0.5, 0.5 }, { 0.25, 0.25, 0.25 } },
	{ 0.717, { 1.0, 0.5, 0.25 }, { 0.1, 0.2, 0.3 } },
	{ -0.2, { 0.5, 0.85, 0.65 }, { 0.75, 0.05, 0.5 } },
};

int main(int argc, char *argv[]) {
	Options options;
	
	pbrtInit(options);
	
	Spectrum etai();
	Spectrum etat();
	Spectrum n();
	Spectrum k();
	
	// Fresnel utility functions.
	for (size_t i = 0; i < sizeof(frDielInputVector)/sizeof(frDielInputVector[0]); i++) {
		Spectrum frd = FrDiel(frDielInputVector[i].cosi, frDielInputVector[i].cost, RGBSpectrum::FromRGB(frDielInputVector[i].etai), RGBSpectrum::FromRGB(frDielInputVector[i].etat));
		fprintf(stdout, "FrDiel In: %2.10f, %2.10f [3]float64{ %2.10f, %2.10f, %2.10f }, [3]float64{ %2.10f, %2.10f, %2.10f },\n", 
			frDielInputVector[i].cosi, frDielInputVector[i].cost, frDielInputVector[i].etai[0], frDielInputVector[i].etai[1], frDielInputVector[i].etai[2], 
			frDielInputVector[i].etat[0], frDielInputVector[i].etat[1], frDielInputVector[i].etat[2]);
		fprintf(stdout, "Output: "); frd.Write(stdout); fprintf(stdout, "\n");
    }     
	for (size_t i = 0; i < sizeof(frCondInputVector)/sizeof(frCondInputVector[0]); i++) {
        Spectrum frc = FrCond(frCondInputVector[i].cosi, RGBSpectrum::FromRGB(frCondInputVector[i].n), RGBSpectrum::FromRGB(frCondInputVector[i].k));
        fprintf(stdout, "FrCond In: %2.10f [3]float64{ %2.10f, %2.10f, %2.10f }, [3]float64{ %2.10f, %2.10f, %2.10f },\n",
			frCondInputVector[i].cosi, frCondInputVector[i].n[0], frCondInputVector[i].n[1], frCondInputVector[i].n[2],
			frCondInputVector[i].k[0], frCondInputVector[i].k[1], frCondInputVector[i].k[2]);
		fprintf(stdout, "Output: "); frc.Write(stdout); fprintf(stdout, "\n");
	}
	
	float pdf;
	float u1 = 0.5;
	float u2 = 0.25;
	Vector wo = UniformSampleHemisphere(u1, u2);
	Vector wi = UniformSampleHemisphere(u1, u2);
	
	Spectrum r = RGBSpectrum(0.333);
	float ei = 0.5;
	float et = 0.25;
	Fresnel *f = new FresnelDielectric(ei, et);
	SpecularReflection spec(r, f);
	
	Spectrum specF = spec.f(wo, wi);	
	float specPdf = spec.Pdf(wo, wi);	
	Spectrum sampF = spec.Sample_f(wo, &wi, u1, u2, &pdf);
	
	fprintf(stdout, "SpecRefl In: r: "); r.Write(stdout); fprintf(stdout, ", fd: %f, %f\n", ei, et);
	fprintf(stdout, "wo: %f,%f,%f  wi: %f,%f,%f\n", wi.x, wi.y, wi.z, wo.x, wo.y, wo.z);
	
	fprintf(stdout, "SpecRefl Out: f: "); specF.Write(stdout); fprintf(stdout, " pdf: %f\n", specPdf);
	fprintf(stdout, "SpecRefl Out: sampF: "); sampF.Write(stdout); fprintf(stdout, " pdf: %f\n", pdf);


	wo = UniformSampleHemisphere(u1, u2);
	wi = UniformSampleHemisphere(u1, u2);
	ei = 1.0;
	et = 1.0;
	
	Spectrum t = RGBSpectrum(0.22);
	SpecularTransmission trans(t, ei, et);
	
	Spectrum transF = trans.f(wo, wi);
	float transPdf = trans.Pdf(wo, wi);
	Spectrum transSampF = trans.Sample_f(wo, &wi, u1, u2, &pdf);

	fprintf(stdout, "SpecTrans In: t: "); t.Write(stdout); fprintf(stdout, ", fd: %f, %f\n", ei, et);
	fprintf(stdout, "wo: %f,%f,%f  wi: %f,%f,%f\n", wi.x, wi.y, wi.z, wo.x, wo.y, wo.z);
	
	fprintf(stdout, "SpecTrans Out: f: "); transF.Write(stdout); fprintf(stdout, " pdf: %f\n", transPdf);
	fprintf(stdout, "SpecTrans Out: sampF: "); transSampF.Write(stdout); fprintf(stdout, " pdf: %f\n", pdf);
	
	pbrtCleanup();
}
