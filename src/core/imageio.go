package core

import (
	"image"
	"image/png"
	"image/color"
	"strings"
	"path/filepath"
	"os"
	"math"
)

func ReadImage(name string) (pixels []Spectrum, xSize, ySize int) {
	ext := strings.ToLower(filepath.Ext(name))
	if strings.Compare(ext, ".exr") == 0 {
		Warning("EXR file load not implemented.")				
	} else if strings.Compare(ext, ".pfm") == 0 {
		Warning("PFM file load not implemented.")						
	} else if strings.Compare(ext, ".tga") == 0 {
		Warning("TGA file load not implemented.")						
	}
	
    Error("Unable to load image stored in format \"%s\" for filename \"%s\". Returning a constant grey image instead.", ext, name)
    ret := make([]Spectrum, 1, 1)
    ret[0] = *CreateSpectrum1(0.5) 
    return ret, 1, 1
}

func WriteImage(name string, pixels, alpha []float64, XRes, YRes, totalXRes, totalYRes, xOffset, yOffset int) {
	ext := strings.ToLower(filepath.Ext(name))
	if strings.Compare(ext, ".exr") == 0 {
		Warning("EXR file save not implemented.")	
	} else if strings.Compare(ext, ".pfm") == 0 {
		Warning("PFM file save not implemented.")			
	} else if strings.Compare(ext, ".tga") == 0 {
		Warning("TGA file save not implemented.")			
	} else if strings.Compare(ext, ".png") == 0 {
		writeImagePng(name, pixels, XRes, YRes)
	}
	Error("Can't determine image file type from suffix of filename \"%s\"", name)
}


func writeImagePng(filename string, pixels []float64, xres, yres int) {
	outImage := image.NewNRGBA(image.Rect(0, 0, xres, yres))
	
	to_byte := func (v float64) uint8 {
		// apply gamma and convert to 0..255
		return uint8(Clamp(255.0 * math.Pow(v, 1.0/2.2), 0.0, 255.0))
	}
	
	for y := 0; y < yres; y++ {
		for x := 0; x < xres; x++ {
            var fcolor color.NRGBA
            fcolor.R = to_byte(pixels[3*(y*xres+x)+2])
            fcolor.G = to_byte(pixels[3*(y*xres+x)+1])
            fcolor.B = to_byte(pixels[3*(y*xres+x)+0])
            fcolor.A = 0xff           
            outImage.Set(x, y, fcolor)			
		}
	}
    f, err := os.Create(filename)
    defer f.Close()

    if err != nil {
        Error("Error writing PNG \"%s\"", filename)
    }
    png.Encode(f, outImage)
	
}