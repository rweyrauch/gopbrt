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
	"github.com/ftrvxmtrx/tga"
	"github.com/rweyrauch/gopbrt/src/openexr"
	"image"
	"image/color"
	"image/png"
	"image/jpeg"
	"math"
	"os"
	"path/filepath"
	"strings"
)

func ReadImage(name string) (pixels []Spectrum, xSize, ySize int) {
	ext := strings.ToLower(filepath.Ext(name))
	if strings.Compare(ext, ".exr") == 0 {
		return readImageExr(name)
	} else if strings.Compare(ext, ".pfm") == 0 {
		Warning("PFM file load not implemented.")
	} else if strings.Compare(ext, ".tga") == 0 {
		return readImageTga(name) 
	} else if strings.Compare(ext, ".png") == 0 {
		return readImagePng(name)
	} else if strings.Compare(ext, ".jpg") == 0 {
		return readImageJpeg(name)
	}

	Error("Unable to load image stored in format \"%s\" for filename \"%s\". Returning a constant grey image instead.", ext, name)
	ret := make([]Spectrum, 1, 1)
	ret[0] = *NewSpectrum1(0.5)
	return ret, 1, 1
}

func WriteImage(name string, pixels, alpha []float32, XRes, YRes, totalXRes, totalYRes, xOffset, yOffset int) {
	ext := strings.ToLower(filepath.Ext(name))
	Debug("Writing image: %s", name)
	if strings.Compare(ext, ".exr") == 0 {
		writeImageExr(name, pixels, XRes, YRes)
	} else if strings.Compare(ext, ".pfm") == 0 {
		Warning("PFM file save not implemented.")
	} else if strings.Compare(ext, ".tga") == 0 {
		writeImageTga(name, pixels, XRes, YRes)
	} else if strings.Compare(ext, ".png") == 0 {
		writeImagePng(name, pixels, XRes, YRes)
	} else if strings.Compare(ext, ".jpg") == 0 {
		writeImageJpeg(name, pixels, XRes, YRes)
	} else {
		Error("Can't determine image file type from suffix of filename \"%s\"", name)
	}
}

func writeImagePng(filename string, pixels []float32, xres, yres int) {
	outImage := image.NewNRGBA(image.Rect(0, 0, xres, yres))

	to_byte := func(v float32) uint8 {
		// apply gamma and convert to 0..255
		return uint8(Clamp(255.0*math.Pow(float64(v), 1.0/2.2), 0.0, 255.0))
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
	} else {
		png.Encode(f, outImage)
	}
}

func readImagePng(filename string) (image []Spectrum, width, height int) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		Error("Error reading PNG \"%s\"", filename)
		return nil, 0, 0
	}

	pngImage, err := png.Decode(f)
	if err != nil {
		Error("Error decoding PNG \"%s\"", filename)
		return nil, 0, 0
	}

	bounds := pngImage.Bounds()

	width = bounds.Dx()
	height = bounds.Dy()
	image = make([]Spectrum, width*height, width*height)

	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			rb, gb, bb, _ := pngImage.At(x, y).RGBA()
			r := float64(rb) / float64(0xffff)
			g := float64(gb) / float64(0xffff)
			b := float64(bb) / float64(0xffff)
			image[y*width+x] = *NewSpectrumRGB(r, g, b)
		}
	}
	return image, width, height
}

func writeImageJpeg(filename string, pixels []float32, xres, yres int) {
	outImage := image.NewNRGBA(image.Rect(0, 0, xres, yres))

	to_byte := func(v float32) uint8 {
		// apply gamma and convert to 0..255
		return uint8(Clamp(255.0*math.Pow(float64(v), 1.0/2.2), 0.0, 255.0))
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
		Error("Error writing JPEG \"%s\"", filename)
	} else {
		options := jpeg.Options{85}
		jpeg.Encode(f, outImage, &options)
	}
}

func readImageJpeg(filename string) (image []Spectrum, width, height int) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		Error("Error reading JPEG \"%s\"", filename)
		return nil, 0, 0
	}

	jpegImage, err := jpeg.Decode(f)
	if err != nil {
		Error("Error decoding JPEG \"%s\"", filename)
		return nil, 0, 0
	}

	bounds := jpegImage.Bounds()

	width = bounds.Dx()
	height = bounds.Dy()
	image = make([]Spectrum, width*height, width*height)

	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			rb, gb, bb, _ := jpegImage.At(x, y).RGBA()
			r := float64(rb) / float64(0xffff)
			g := float64(gb) / float64(0xffff)
			b := float64(bb) / float64(0xffff)
			image[y*width+x] = *NewSpectrumRGB(r, g, b)
		}
	}
	return image, width, height
}

func writeImageTga(filename string, pixels []float32, xres, yres int) {
	outImage := image.NewNRGBA(image.Rect(0, 0, xres, yres))

	to_byte := func(v float32) uint8 {
		// apply gamma and convert to 0..255
		return uint8(Clamp(255.0*math.Pow(float64(v), 1.0/2.2), 0.0, 255.0))
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
		Error("Error writing TGA \"%s\"", filename)
	} else {
		tga.Encode(f, outImage)
	}
}

func readImageTga(filename string) (image []Spectrum, width, height int) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		Error("Error reading TGA \"%s\"", filename)
		return nil, 0, 0
	}

	tgaImage, err := tga.Decode(f)
	if err != nil {
		Error("Error decoding TGA \"%s\"", filename)
		return nil, 0, 0
	}

	bounds := tgaImage.Bounds()

	width = bounds.Dx()
	height = bounds.Dy()
	image = make([]Spectrum, width*height, width*height)

	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			rb, gb, bb, _ := tgaImage.At(x, y).RGBA()
			r := float64(rb) / float64(0xffff)
			g := float64(gb) / float64(0xffff)
			b := float64(bb) / float64(0xffff)
			image[y*width+x] = *NewSpectrumRGB(r, g, b)
		}
	}
	return image, width, height
}

func readImageExr(filename string) (image []Spectrum, width, height int) {
	var planarRGBA []float32
	var channels int
	width, height, channels, planarRGBA = openexr.ReadImageEXR(filename)
	if planarRGBA != nil {
		Debug("Read EXR image (%dx%dx%d)", width, height, channels)
		image = make([]Spectrum, width*height, width*height)
		if channels == 3 || channels == 4 {
			for c := 0; c < 3; c++ { // ignore 4th channel
				chanStart := c * width * height
				for i, _ := range image {
					image[i].c[c] = float64(planarRGBA[i+chanStart])
				}
			}
		} else { // channels == 1 || channels == 2 (ignore second channel)
			for i, _ := range image {
				image[i].c[0] = float64(planarRGBA[i])
				image[i].c[1] = float64(planarRGBA[i])
				image[i].c[2] = float64(planarRGBA[i])
			}
		}
	}
	return image, width, height
}

func writeImageExr(filename string, pixels []float32, xres, yres int) {
	Debug("Write EXR image (%dx%dx%d)", xres, yres, 3)
	openexr.WriteImageEXR(filename, xres, yres, 3, pixels)
}
