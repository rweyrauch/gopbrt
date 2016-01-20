package tinyexr 

/*
#include "tinyexr.h"
#include <stdlib.h>
*/
import "C"
import "unsafe"
import "fmt"

const (
 TINYEXR_PIXELTYPE_UINT = 0
 TINYEXR_PIXELTYPE_HALF = 1
 TINYEXR_PIXELTYPE_FLOAT = 2
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
	rpt :=(*[4]C.int)(unsafe.Pointer(exrImage.requested_pixel_types))
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


	return width, height, channels, planarRGBA
}
