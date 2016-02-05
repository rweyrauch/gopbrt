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
package tinyexr

/*
#include <stdlib.h>
#include <string.h>
#include "tinyexr.h"
*/
import "C"
import "unsafe"
import "fmt"
import "reflect"

const (
	TINYEXR_PIXELTYPE_UINT  = 0
	TINYEXR_PIXELTYPE_HALF  = 1
	TINYEXR_PIXELTYPE_FLOAT = 2
)

const (
	TINYEXR_COMPRESSIONTYPE_NONE = 0
	TINYEXR_COMPRESSIONTYPE_RLE = 1  // not supported yet
	TINYEXR_COMPRESSIONTYPE_ZIPS = 2
	TINYEXR_COMPRESSIONTYPE_ZIP = 3
	TINYEXR_COMPRESSIONTYPE_PIZ = 4
)

func initEXRImage(exrImage *C.EXRImage) {
	C.InitEXRImage(exrImage)
}
func freeEXRImage(exrImage *C.EXRImage) {
	C.FreeEXRImage(exrImage)
}
func parseMultiChannelEXRHeaderFromFile(image *C.EXRImage, filename string) (int, string) {
	cfilename := C.CString(filename)
	var cerr *C.char
	res := C.ParseMultiChannelEXRHeaderFromFile(image, cfilename, &cerr)
	C.free(unsafe.Pointer(cfilename))
	return int(res), C.GoString(cerr)
}

func loadMultiChannelEXRFromFile(image *C.EXRImage, filename string) (int, string) {
	cfilename := C.CString(filename)
	var cerr *C.char
	res := C.LoadMultiChannelEXRFromFile(image, cfilename, &cerr)
	C.free(unsafe.Pointer(cfilename))
	return int(res), C.GoString(cerr)
}

func saveMultiChannelEXRToFile(image *C.EXRImage, filename string) (int, string) {
	cfilename := C.CString(filename)
	var cerr *C.char
	res := C.SaveMultiChannelEXRToFile(image, cfilename, &cerr)
	C.free(unsafe.Pointer(cfilename))
	return int(res), C.GoString(cerr)
}

func ReadImageEXR(filename string) (width, height, channels int, planarRGBA []float32) {
	var exrImage C.EXRImage

	fmt.Printf("Reading EXR image from %s\n", filename)
	initEXRImage(&exrImage)
	defer freeEXRImage(&exrImage)

	res, errmsg := parseMultiChannelEXRHeaderFromFile(&exrImage, filename)
	if res != 0 {
		fmt.Printf("Failed to parse EXR header from file %s.  Error: %s\n", filename, errmsg)
		return 0, 0, 0, nil
	}

	// request that image load each channel as floats
	rpt := (*[4]C.int)(unsafe.Pointer(exrImage.requested_pixel_types))
	for c := 0; c < int(exrImage.num_channels); c++ {
		rpt[c] = TINYEXR_PIXELTYPE_FLOAT
	}

	res, errmsg = loadMultiChannelEXRFromFile(&exrImage, filename)
	if res != 0 {
		fmt.Printf("Failed to load EXR from file %s.  Error: %s\n", filename, errmsg)
		return 0, 0, 0, nil
	}

	width = int(exrImage.width)
	height = int(exrImage.height)
	channels = int(exrImage.num_channels)
	totalValues := width * height * channels

	planarRGBA = make([]float32, totalValues, totalValues)
	sliceHeader := (*reflect.SliceHeader)(unsafe.Pointer(&planarRGBA))
	C.memcpy(unsafe.Pointer(sliceHeader.Data), unsafe.Pointer(exrImage.images), C.size_t(4*totalValues))

	return width, height, channels, planarRGBA
}

func WriteImageEXR(filename string, width, height, channels int, planarRGBA []float32) {
	var exrImage C.EXRImage

	if channels != 3 {
		panic("Only 3-channel export implemented.")
	}
	
	initEXRImage(&exrImage)
	defer freeEXRImage(&exrImage)

	exrImage.width = C.int(width)
	exrImage.height = C.int(height)
	exrImage.num_channels = C.int(channels)
	exrImage.compression = TINYEXR_COMPRESSIONTYPE_ZIP
	pixel_types := make([]C.int, channels, channels)
	requested_pixel_types := make([]C.int, channels, channels)
	for i := 0; i < channels; i++ {
		pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT
		requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF
	}
	exrImage.pixel_types = (*C.int)(unsafe.Pointer(&pixel_types[0]))
	exrImage.requested_pixel_types = (*C.int)(unsafe.Pointer(&requested_pixel_types[0]))
	
	imageR := make([]float32, width*height, width*height)
	imageG := make([]float32, width*height, width*height)
	imageB := make([]float32, width*height, width*height)
	for i, _ := range imageR {
		imageR[i] = planarRGBA[i*channels+0]
		imageB[i] = planarRGBA[i*channels+1]
		imageG[i] = planarRGBA[i*channels+2]
	}
	imagePtrs := make([]*float32, channels, channels)
	imagePtrs[0] = &imageR[0]
	imagePtrs[1] = &imageG[0]
	imagePtrs[2] = &imageB[0]
	exrImage.images = (**C.uchar)(unsafe.Pointer(&imagePtrs[0]))
	
	res, errmsg := saveMultiChannelEXRToFile(&exrImage, filename)
	if res != 0 {
		fmt.Printf("Failed to save EXR to file %s.  Error: %s\n", filename, errmsg)
	}
}
