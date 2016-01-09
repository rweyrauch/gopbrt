package pbrt

type Vector struct {
    x, y, z float64
}

type Point struct {
    x, y, z float64
}

type Normal struct {
    x, y, z float64
}

type Ray struct {
    o Point
    d Vector
    mint, maxt float64
    time float64
    depth int
}

type RayDifferential struct {
    Ray
    hasDifferentials bool
    rxOrigin, ryOrigin Point
    rxDirection, ryDirection Vector
}

type BBox struct {
    pMin, pMax Point
}
