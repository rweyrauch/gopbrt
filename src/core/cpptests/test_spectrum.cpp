#include "api.h"
#include "pbrt.h"
#include "spectrum.h"
#include <stdio.h>

const float rgbInputVector[][3] = {
	{ 1.0, 1.0, 1.0 },
	{ 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0 },
	{ 0.5, 0.5, 0.0 },
	{ 0.0, 0.75, 1.0 },
};

const float blackbodyInputVector[][2] = {	// {temp, scale}
	{ 20000.0, 0.8 },
	{ 5500.0, 0.5 },
	{ 12005.0, 1.0 },
	{ 0.0, 1.0 },
	{ 0.0, 0.0 },
	{ 400.0, 0.5 },
};

int main(int argc, char *argv[]) {
	Options options;
	
	pbrtInit(options);
	
	RGBSpectrum rgb;
	float xyz[3];
	
	const size_t numInputs = sizeof(rgbInputVector) / sizeof(float[3]);
	for (size_t i = 0; i < numInputs; i++) {
		rgb = RGBSpectrum::FromRGB(rgbInputVector[i]);
		fprintf(stdout, "RGBin: %2.10f,%2.10f,%2.10f RGBout: ", rgbInputVector[i][0], rgbInputVector[i][1], rgbInputVector[i][2]);
		rgb.Write(stdout);
		rgb.ToXYZ(xyz);
		fprintf(stdout, " XYZout: %2.10f,%2.10f,%2.10f  Y: %2.10f\n", xyz[0], xyz[1], xyz[2], rgb.y());
	}
	
	const size_t numBbInputs = sizeof(blackbodyInputVector) / sizeof(float[2]);
	
    float *v = new float[nCIESamples];	
	for (size_t i = 0; i < numBbInputs; i++) {
        Blackbody(CIE_lambda, nCIESamples, blackbodyInputVector[i][0], v);
        rgb = blackbodyInputVector[i][1] * Spectrum::FromSampled(CIE_lambda, v, nCIESamples);	
        fprintf(stdout, "BlackBodyin: T:%2.10f Sc:%2.10f\n", blackbodyInputVector[i][0], blackbodyInputVector[i][1]);
        fprintf(stdout, "BBout: ");
        rgb.Write(stdout);
        fprintf(stdout, "\n"); 	
	}
	pbrtCleanup();
}
