package pbrt

type Filter interface {
	Evaluate(x, y float64) float64
    XWidth() float64
    YWidth() float64
    InvXWidth() float64
    InvYWidth() float64
}

type FilterData struct {
    xWidth, yWidth float64
    invXWidth, invYWidth float64
}