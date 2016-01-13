package pbrt

type Film interface {
    AddSample(sample *CameraSample, L *Spectrum)
    Splat(sample *CameraSample, L *Spectrum) 
    GetSampleExtent() (xstart, xend, ystart, yend int)
    GetPixelExtent() (xstart, xend, ystart, yend int)
    UpdateDisplay(x0, y0, x1, y1 int, splatScale float64)
    WriteImage(splatScale float64)
    XResolution() int
    YResolution() int
}