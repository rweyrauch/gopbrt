package pbrt

type Scene struct {
    aggregate Primitive
    lights []Light
    volumeRegion VolumeRegion
    bound BBox 
}