package core

import (
	"math"
)

const (
	filterTableSize = 16
)

type imagePixel struct {
	Lxyz      [3]float32
	weightSum float32
	splatXYZ  [3]float32
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

func CreateImageFilm(xres, yres int, filter Filter, crop []float64, filename string, openWindow bool) *ImageFilm {
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
		}
	}

	// Possibly open window for image display
	if openWindow || options.OpenWindow {
		Warning("Support for opening image display window not available in this build.")
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
	ifx := make([]int, 0, x1-x0+1)
	for x := x0; x <= x1; x++ {
		fx := math.Abs((float64(x) - dimageX) * f.filter.InvXWidth() * filterTableSize)
		ifx[x-x0] = Mini(Floor2Int(fx), filterTableSize-1)
	}
	ify := make([]int, y1-y0+1)
	for y := y0; y <= y1; y++ {
		fy := math.Abs((float64(y) - dimageY) * f.filter.InvYWidth() * filterTableSize)
		ify[y-y0] = Mini(Floor2Int(fy), filterTableSize-1)
	}
	syncNeeded := (f.filter.XWidth() > 0.5 || f.filter.YWidth() > 0.5)
	for y := y0; y <= y1; y++ {
		for x := x0; x <= x1; x++ {
			// Evaluate filter value at $(x,y)$ pixel
			offset := ify[y-y0]*filterTableSize + ifx[x-x0]
			filterWt := float32(f.filterTable[offset])

			// Update pixel values with filtered sample contribution
			pixel := f.pixelAt(x-f.xPixelStart, y-f.yPixelStart)
			if !syncNeeded {
				pixel.Lxyz[0] += filterWt * xyz[0]
				pixel.Lxyz[1] += filterWt * xyz[1]
				pixel.Lxyz[2] += filterWt * xyz[2]
				pixel.weightSum += filterWt
			} else {
				// Safely update _Lxyz_ and _weightSum_ even with concurrency
				AtomicAddf(&pixel.Lxyz[0], filterWt*xyz[0])
				AtomicAddf(&pixel.Lxyz[1], filterWt*xyz[1])
				AtomicAddf(&pixel.Lxyz[2], filterWt*xyz[2])
				AtomicAddf(&pixel.weightSum, filterWt)
			}
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
	AtomicAddf(&pixel.splatXYZ[0], xyz[0])
	AtomicAddf(&pixel.splatXYZ[1], xyz[1])
	AtomicAddf(&pixel.splatXYZ[2], xyz[2])
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

func (f *ImageFilm) UpdateDisplay(x0, y0, x1, y1 int, splatScale float32) {

}

func (f *ImageFilm) WriteImage(splatScale float32) {
	// Convert image to RGB and compute final pixel values
	nPix := f.xPixelCount * f.yPixelCount
	rgb := make([]float32, 3*nPix, 3*nPix)

	offset := 0
	for y := 0; y < f.yPixelCount; y++ {
		for x := 0; x < f.xPixelCount; x++ {
			// Convert pixel XYZ color to RGB
			pixel := f.pixelAt(x, y)

			tmprgb := XYZToRGB(pixel.Lxyz)
			rgb[3*offset] = tmprgb[0]
			rgb[3*offset+1] = tmprgb[1]
			rgb[3*offset+2] = tmprgb[2]

			// Normalize pixel with weight sum
			weightSum := pixel.weightSum
			if weightSum != 0.0 {
				invWt := 1.0 / weightSum
				rgb[3*offset] = float32(math.Max(0.0, float64(rgb[3*offset]*invWt)))
				rgb[3*offset+1] = float32(math.Max(0.0, float64(rgb[3*offset+1]*invWt)))
				rgb[3*offset+2] = float32(math.Max(0.0, float64(rgb[3*offset+2]*invWt)))
			}

			// Add splat value at pixel
			splatRGB := XYZToRGB(pixel.splatXYZ)
			rgb[3*offset] += splatScale * splatRGB[0]
			rgb[3*offset+1] += splatScale * splatRGB[1]
			rgb[3*offset+2] += splatScale * splatRGB[2]
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
	openwin := params.FindBoolParam("display", false)
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
	if options.QuickRender {
		xres = Maxi(1, xres/4)
		yres = Maxi(1, yres/4)
	}

	return CreateImageFilm(xres, yres, filter, crop, filename, openwin)
}
