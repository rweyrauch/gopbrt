package pbrt

type Camera interface {
	GenerateRay(sample *CameraSample) (ray *Ray, weight float64)
	GenerateRayDifferential(sample *CameraSample) (ray *RayDifferential, weight float64)
}

type CameraCore struct {
	CameraToWorld             *AnimatedTransform
	ShutterOpen, ShutterClose float64
	Film                      Film
}


func GenerateRayDifferential(camera Camera, sample *CameraSample) (rd *RayDifferential, weight float64) {
    ray, wt := camera.GenerateRay(sample)
	rd = CreateRayDifferentialFromRay(ray)
    // Find ray after shifting one pixel in the $x$ direction
    sshift := sample
    sshift.imageX++
    	
    rx, wtx := camera.GenerateRay(sshift)
    rd.rxOrigin = rx.origin
    rd.rxDirection = rx.dir

    // Find ray after shifting one pixel in the $y$ direction
    sshift.imageX--
    sshift.imageY++
    ry, wty := camera.GenerateRay(sshift)
    rd.ryOrigin = ry.origin
    rd.ryDirection = ry.dir
    if wtx == 0.0 || wty == 0.0 { return rd, 0.0 }
    rd.hasDifferentials = true
    
    return rd, wt
}

type ProjectiveCamera struct {
	CameraCore
	CameraToScreen, RasterToCamera *Transform
	ScreenToRaster, RasterToScreen *Transform
	LensRadius, FocalDistance      float64
}

func CreateProjectiveCamera(cam2world *AnimatedTransform, proj *Transform, screenWindow [4]float64, sopen, sclose, lensr, focald float64, film Film) *ProjectiveCamera {

	camera := &ProjectiveCamera{CameraCore{cam2world, sopen, sclose, film}, proj, nil, nil, nil, lensr, focald}

	// Compute projective camera screen transformations
	camera.ScreenToRaster = ScaleTransform(float64(film.XResolution()), float64(film.YResolution()), 1.0).MultTransform(
		ScaleTransform(1.0/(screenWindow[1]-screenWindow[0]), 1.0/(screenWindow[2]-screenWindow[3]), 1.0)).MultTransform(
		TranslateTransform(&Vector{-screenWindow[0], -screenWindow[3], 0.0}))
	camera.RasterToScreen = InverseTransform(camera.ScreenToRaster)
	camera.RasterToCamera = InverseTransform(camera.CameraToScreen).MultTransform(camera.RasterToScreen)

	return camera
}