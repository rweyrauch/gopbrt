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

const (
	filterTableSize = 16
)

type imagePixel struct {
	Lxyz      [3]float64
	weightSum float64
	splatXYZ  [3]float64
}

type ImageFilm struct {
	FilmData
	filter                                             Filter
	cropWindow                                         []float64
	filename                                           string
	xPixelStart, yPixelStart, xPixelCount, yPixelCount int

	pixels      []imagePixel
	filterTable [filterTableSize * filterTableSize]float64
}

func NewImageFilm(xres, yres int, filter Filter, crop []float64, filename string) *ImageFilm {
	film := new(ImageFilm)
	film.xResolution = xres
	film.yResolution = yres
	film.filter = filter
	film.cropWindow = crop
	film.filename = filename

	film.xPixelStart = Ceil2Int(float64(film.xResolution) * film.cropWindow[0])
	film.xPixelCount = Maxi(1, Ceil2Int(float64(film.xResolution)*film.cropWindow[1])-film.xPixelStart)
	film.yPixelStart = Ceil2Int(float64(film.yResolution) * film.cropWindow[2])
	film.yPixelCount = Maxi(1, Ceil2Int(float64(film.yResolution)*film.cropWindow[3])-film.yPixelStart)

	// Allocate film image storage
	film.pixels = make([]imagePixel, film.xPixelCount*film.yPixelCount, film.xPixelCount*film.yPixelCount)

	// Precompute filter weight table
	i := 0
	for y := 0; y < filterTableSize; y++ {
		fy := (float64(y) + 0.5) * filter.YWidth() / float64(filterTableSize)
		for x := 0; x < filterTableSize; x++ {
			fx := (float64(x) + 0.5) * filter.XWidth() / float64(filterTableSize)
			film.filterTable[i] = filter.Evaluate(fx, fy)
			i++
		}
	}

	return film
}

func (f *ImageFilm) pixelAt(x, y int) *imagePixel {
	return &f.pixels[y*f.xPixelCount+x]
}

func (f *ImageFilm) AddSample(sample *Sample, L *Spectrum) {
	// Compute sample's raster extent
	dimageX := sample.imageX - 0.5
	dimageY := sample.imageY - 0.5
	x0 := Ceil2Int(dimageX - f.filter.XWidth())
	x1 := Floor2Int(dimageX + f.filter.XWidth())
	y0 := Ceil2Int(dimageY - f.filter.YWidth())
	y1 := Floor2Int(dimageY + f.filter.YWidth())
	x0 = Maxi(x0, f.xPixelStart)
	x1 = Mini(x1, f.xPixelStart+f.xPixelCount-1)
	y0 = Maxi(y0, f.yPixelStart)
	y1 = Mini(y1, f.yPixelStart+f.yPixelCount-1)
	if (x1-x0) < 0 || (y1-y0) < 0 {
		//PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(const_cast<CameraSample *>(&sample))
		return
	}

	// Loop over filter support and add sample to pixel arrays
	xyz := L.ToXYZ()

	// Precompute $x$ and $y$ filter table offsets
	ifx := make([]int, x1-x0+1, x1-x0+1)
	for x := x0; x <= x1; x++ {
		fx := math.Abs((float64(x) - dimageX) * f.filter.InvXWidth() * filterTableSize)
		ifx[x-x0] = Mini(Floor2Int(fx), filterTableSize-1)
	}
	ify := make([]int, y1-y0+1, y1-y0+1)
	for y := y0; y <= y1; y++ {
		fy := math.Abs((float64(y) - dimageY) * f.filter.InvYWidth() * filterTableSize)
		ify[y-y0] = Mini(Floor2Int(fy), filterTableSize-1)
	}

	for y := y0; y <= y1; y++ {
		for x := x0; x <= x1; x++ {
			// Evaluate filter value at $(x,y)$ pixel
			offset := ify[y-y0]*filterTableSize + ifx[x-x0]
			filterWt := f.filterTable[offset]

			// Update pixel values with filtered sample contribution
			pixel := f.pixelAt(x-f.xPixelStart, y-f.yPixelStart)
			pixel.Lxyz[0] += filterWt * xyz[0]
			pixel.Lxyz[1] += filterWt * xyz[1]
			pixel.Lxyz[2] += filterWt * xyz[2]
			pixel.weightSum += filterWt
		}
	}
}

func (f *ImageFilm) Splat(sample *Sample, L *Spectrum) {
	if L.HasNaNs() {
		Warning("ImageFilm ignoring splatted spectrum with NaN values.")
		return
	}
	xyz := L.ToXYZ()

	x, y := Floor2Int(sample.imageX), Floor2Int(sample.imageY)
	if x < f.xPixelStart || x-f.xPixelStart >= f.xPixelCount ||
		y < f.yPixelStart || y-f.yPixelStart >= f.yPixelCount {
		return
	}

	pixel := f.pixelAt(x-f.xPixelStart, y-f.yPixelStart)
	pixel.splatXYZ[0] += xyz[0]
	pixel.splatXYZ[1] += xyz[1]
	pixel.splatXYZ[2] += xyz[2]
}

func (f *ImageFilm) GetSampleExtent() (xstart, xend, ystart, yend int) {
	xstart = Floor2Int(float64(f.xPixelStart) + 0.5 - f.filter.XWidth())
	xend = Ceil2Int(float64(f.xPixelStart) - 0.5 + float64(f.xPixelCount) + f.filter.XWidth())

	ystart = Floor2Int(float64(f.yPixelStart) + 0.5 - f.filter.YWidth())
	yend = Ceil2Int(float64(f.yPixelStart) - 0.5 + float64(f.yPixelCount) + f.filter.YWidth())

	return xstart, xend, ystart, yend
}

func (f *ImageFilm) GetPixelExtent() (xstart, xend, ystart, yend int) {
	xstart = f.xPixelStart
	xend = f.xPixelStart + f.xPixelCount
	ystart = f.yPixelStart
	yend = f.yPixelStart + f.yPixelCount

	return xstart, xend, ystart, yend
}

func (f *ImageFilm) UpdateDisplay(x0, y0, x1, y1 int, splatScale float64) {

}

func (f *ImageFilm) WriteImage(splatScale float64) {
	// Convert image to RGB and compute final pixel values
	nPix := f.xPixelCount * f.yPixelCount
	rgb := make([]float32, 3*nPix, 3*nPix)

	offset := 0
	for y := 0; y < f.yPixelCount; y++ {
		for x := 0; x < f.xPixelCount; x++ {
			// Convert pixel XYZ color to RGB
			pixel := f.pixelAt(x, y)

			tmprgb := XYZToRGB(pixel.Lxyz)
			rgb[3*offset] = float32(tmprgb[0])
			rgb[3*offset+1] = float32(tmprgb[1])
			rgb[3*offset+2] = float32(tmprgb[2])

			// Normalize pixel with weight sum
			weightSum := pixel.weightSum
			if weightSum != 0.0 {
				invWt := float32(1.0 / weightSum)
				rgb[3*offset] = float32(math.Max(0.0, float64(rgb[3*offset]*invWt)))
				rgb[3*offset+1] = float32(math.Max(0.0, float64(rgb[3*offset+1]*invWt)))
				rgb[3*offset+2] = float32(math.Max(0.0, float64(rgb[3*offset+2]*invWt)))
			}

			// Add splat value at pixel
			splatRGB := XYZToRGB(pixel.splatXYZ)
			rgb[3*offset] += float32(splatScale * splatRGB[0])
			rgb[3*offset+1] += float32(splatScale * splatRGB[1])
			rgb[3*offset+2] += float32(splatScale * splatRGB[2])

			offset++
		}
	}

	// Write RGB image
	WriteImage(f.filename, rgb, nil, f.xPixelCount, f.yPixelCount, f.xResolution, f.yResolution, f.xPixelStart, f.yPixelStart)
}

func (f *ImageFilm) XResolution() int {
	return f.xResolution
}
func (f *ImageFilm) YResolution() int {
	return f.yResolution
}

func CreateImageFilmFromParams(params *ParamSet, filter Filter) *ImageFilm {
	filename := params.FindStringParam("filename", "")
	xres := params.FindIntParam("xresolution", 640)
	yres := params.FindIntParam("yresolution", 480)
	crop := params.FindFloatArrayParam("cropwindow")
	if crop == nil {
		crop = []float64{0, 1, 0, 1}
	}
	if len(options.ImageFile) != 0 {
		if len(filename) != 0 {
			Warning("Output filename supplied on command line, \"%s\", ignored due to filename provided in scene description file, \"%s\".",
				options.ImageFile, filename)
		} else {
			filename = options.ImageFile
		}
	}
	if len(filename) == 0 {
		filename = "pbrt.exr"
	}
	if options.FastRender {
		xres = Maxi(1, xres/2)
		yres = Maxi(1, yres/2)
	} else if options.QuickRender {
		xres = Maxi(1, xres/4)
		yres = Maxi(1, yres/4)
	}

	return NewImageFilm(xres, yres, filter, crop, filename)
}
