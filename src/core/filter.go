package pbrt

type Filter interface {
	Evaulate(x, y float64) float64
}