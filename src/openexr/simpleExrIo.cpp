/*
	gopbrt

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
#include <vector>
#include <cstdio>
#include "simpleExrIo.h"
#include <ImfRgba.h>
#include <ImfRgbaFile.h>

Pixel *ReadImageEXR(const char* filename, int *width, int *height) {
    using namespace Imf;
    using namespace Imath;
    try {
        RgbaInputFile file(filename);
        Box2i dw = file.dataWindow();
        *width = dw.max.x - dw.min.x + 1;
        *height = dw.max.y - dw.min.y + 1;
        std::vector<Rgba> pixels(*width * *height);
        file.setFrameBuffer(&pixels[0] - dw.min.x - dw.min.y * *width, 1,
                            *width);
        file.readPixels(dw.min.y, dw.max.y);

        Pixel *ret = new Pixel[*width * *height];
        for (int i = 0; i < *width * *height; ++i) {
            ret[i].r = pixels[i].r;
            ret[i].g = pixels[i].g;
            ret[i].b = pixels[i].b;
            ret[i].a = pixels[i].a;
        }
        fprintf(stdout, "Read EXR image %s (%d x %d)", filename, *width, *height);
        return ret;
    } catch (const std::exception &e) {
        fprintf(stderr, "Unable to read image file \"%s\": %s", filename, e.what());
    }

	return NULL;
}

void WriteImageEXR(const char* filename, int width, int height, const Pixel* pixels)
{
    using namespace Imf;
    using namespace Imath;

    Rgba *hrgba = new Rgba[width * height];
    for (int i = 0; i < width * height; ++i)
        hrgba[i] = Rgba(pixels[i].r, pixels[i].g, pixels[i].b);

    int totalXRes = width;
    int totalYRes = height;
    int xOffset = 0;
    int yOffset = 0;

    Box2i displayWindow(V2i(0, 0), V2i(totalXRes - 1, totalYRes - 1));
    Box2i dataWindow(V2i(xOffset, yOffset),
                     V2i(xOffset + width - 1, yOffset + height - 1));

    fprintf(stdout, "Write EXR image: %s (%d x %d)", filename, width, height);
    try {
        RgbaOutputFile file(filename, displayWindow, dataWindow,
                            WRITE_RGBA);
        file.setFrameBuffer(hrgba - xOffset - yOffset * width, 1, width);
        file.writePixels(height);
    } catch (const std::exception &exc) {
    	fprintf(stderr, "Error writing \"%s\": %s", filename, exc.what());
    }

    delete[] hrgba;
}

