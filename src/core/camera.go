package pbrt

type Camera interface {
    GenerateRay(sample *CameraSample) (ray *Ray, weight float32)
    GenerateRayDifferential(sample *CameraSample) (ray *RayDifferential, weight float32)
}

type CameraCore struct {
    CameraToWorld *AnimatedTransform
    ShutterOpen, ShutterClose float32
    Film *Film
}

type ProjectiveCamera struct {
    CameraCore
    CameraToScreen, RasterToCamera *Transform
    ScreenToRaster, RasterToScreen *Transform
    LensRadius, FocalDistance float32
}

func CreateProjectiveCamera(cam2world *AnimatedTransform, proj *Transform, screenWindow [4]float32,
    sopen, sclose, lensr, focald float32, film *Film) *ProjectiveCamera {

    camera := &ProjectiveCamera{ cam2world, sopen, sclose, film, proj, nil, nil, nil, lensr, focald }

    return camera
}
