package core

type Film interface {
	AddSample(sample *Sample, L *Spectrum)
	Splat(sample *Sample, L *Spectrum)
	GetSampleExtent() (xstart, xend, ystart, yend int)
	GetPixelExtent() (xstart, xend, ystart, yend int)
	UpdateDisplay(x0, y0, x1, y1 int, splatScale float32)
	WriteImage(splatScale float32)
	XResolution() int
	YResolution() int
}

type FilmData struct {
    xResolution, yResolution int
}
