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
package openexr

// #cgo CPPFLAGS: -I/usr/include/OpenEXR
// #cgo LDFLAGS: -lIlmImf -lImath -lHalf
// #include "simpleExrIo.h"
// #include <stdlib.h>
import "C"
import "unsafe"

func ReadImageEXR(filename string) (width, height, channels int, packedPixels []float32) {
	
	var pixels *C.Pixel
	var cpixel C.Pixel
	var cwidth, cheight C.int
	cfilename := C.CString(filename)
	
	pixels = C.ReadImageEXR(cfilename, &cwidth, &cheight)
	
	width = int(cwidth)
	height = int(cheight)
	channels = 4
	packedPixels = make([]float32, width*height*channels, width*height*channels)
	pixelPtr := unsafe.Pointer(pixels)
	for i := 0; i < width*height; i++ {
		curPixel := (*C.Pixel)(pixelPtr)
		packedPixels[i*4] = float32(curPixel.r)
		packedPixels[i*4+1] = float32(curPixel.g)
		packedPixels[i*4+2] = float32(curPixel.b)
		packedPixels[i*4+3] = float32(curPixel.a)
		pixelPtr = unsafe.Pointer(uintptr(pixelPtr) + unsafe.Sizeof(cpixel))
	}
	return width, height, channels, packedPixels
}

func WriteImageEXR(filename string, width, height, channels int, packedPixels []float32) {	
	cfilename := C.CString(filename)
	cwidth := C.int(width)
	cheight := C.int(height)	
	var cpixel C.Pixel
	
	cpixelsptr := C.malloc(C.size_t(cwidth*cheight)*C.size_t(unsafe.Sizeof(cpixel)))
	cpixels := (*C.Pixel)(cpixelsptr)
	
	for i := 0; i < width*height; i++ {
		curPixel := (*C.Pixel)(cpixelsptr)
		curPixel.r = C.float(packedPixels[i*channels])
		curPixel.g = C.float(packedPixels[i*channels+1])
		curPixel.b = C.float(packedPixels[i*channels+2])
		curPixel.a = 1.0
		cpixelsptr = unsafe.Pointer(uintptr(cpixelsptr) + unsafe.Sizeof(cpixel))
	}
	
	C.WriteImageEXR(cfilename, cwidth, cheight, cpixels)
	
	C.free(unsafe.Pointer(cpixels))
}
	