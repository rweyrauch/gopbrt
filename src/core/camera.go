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

type Camera interface {
	GenerateRay(sample *Sample) (ray *Ray, weight float64)
	GenerateRayDifferential(sample *Sample) (ray *RayDifferential, weight float64)
    CameraToWorld() *AnimatedTransform
    ShutterOpen() float64
    ShutterClose() float64
    Film() Film
}

type CameraCore struct {
	cameraToWorld             *AnimatedTransform
	shutterOpen, shutterClose float64
	film                      Film
}

func GenerateRayDifferential(camera Camera, sample *Sample) (rd *RayDifferential, weight float64) {
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
	if wtx == 0.0 || wty == 0.0 {
		return rd, 0.0
	}
	rd.hasDifferentials = true

	return rd, wt
}

type ProjectiveCamera struct {
	CameraCore
	cameraToScreen, rasterToCamera *Transform
	screenToRaster, rasterToScreen *Transform
	lensRadius, focalDistance      float64
}

type (
	PerspectiveCamera struct {
		ProjectiveCamera
		dxCamera, dyCamera Vector
	}

	EnvironmentCamera struct {
		CameraCore
	}

	OrthoCamera struct {
		ProjectiveCamera
		dxCamera, dyCamera Vector
	}
)

func NewPerspectiveCamera(cam2world *AnimatedTransform, screenWindow [4]float64, shutteropen,
        shutterclose, lensradius, focaldistance, fov float64, film Film) *PerspectiveCamera {
	camera := new(PerspectiveCamera)
	
	// camera core members
	camera.cameraToWorld = cam2world
	camera.shutterOpen = shutteropen
	camera.shutterClose = shutterclose
	camera.film = film
	
	// Compute projective camera screen transformations
	camera.cameraToScreen = PerspectiveTransform(fov, 1.0e-2, 1000.0)
	camera.screenToRaster = ScaleTransform(float64(film.XResolution()), float64(film.YResolution()), 1.0).MultTransform(
		ScaleTransform(1.0/(screenWindow[1]-screenWindow[0]), 1.0/(screenWindow[2]-screenWindow[3]), 1.0)).MultTransform(
		TranslateTransform(&Vector{-screenWindow[0], -screenWindow[3], 0.0}))
	camera.rasterToScreen = InverseTransform(camera.screenToRaster)
	camera.rasterToCamera = InverseTransform(camera.cameraToScreen).MultTransform(camera.rasterToScreen)
	camera.lensRadius = lensradius
	camera.focalDistance = focaldistance
	
    // Compute differential changes in origin for perspective camera rays
    camera.dxCamera = *PointTransform(camera.rasterToCamera, CreatePoint(1,0,0)).Sub(PointTransform(camera.rasterToCamera, CreatePoint(0,0,0)))
    camera.dyCamera = *PointTransform(camera.rasterToCamera, CreatePoint(0,1,0)).Sub(PointTransform(camera.rasterToCamera, CreatePoint(0,0,0)))
    
	return camera
}

func (c *PerspectiveCamera) GenerateRay(sample *Sample) (ray *Ray, weight float64) {
    // Generate raster and camera samples
    Pras := CreatePoint(sample.imageX, sample.imageY, 0)
    Pcamera := PointTransform(c.rasterToCamera, Pras)
    ray = CreateRay(CreatePoint(0,0,0), NormalizeVector(CreateVectorFromPoint(Pcamera)), 0.0, INFINITY, 0.0, 0)
    // Modify ray for depth of field
    if c.lensRadius > 0.0 {
        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        // Compute point on plane of focus
        ft := c.focalDistance / ray.dir.Z
        Pfocus := ray.PointAt(ft)

        // Update ray for effect of lens
        ray.origin = *CreatePoint(lensU, lensV, 0.0)
        ray.dir = *NormalizeVector(Pfocus.Sub(&ray.origin))
    }
    ray.time = sample.time
    ray = RayAnimatedTransform(c.cameraToWorld, ray)
    
    return ray, 1.0
}

func (c *PerspectiveCamera) GenerateRayDifferential(sample *Sample) (ray *RayDifferential, weight float64) {
   // Generate raster and camera samples
    Pras := CreatePoint(sample.imageX, sample.imageY, 0)
    Pcamera := PointTransform(c.rasterToCamera, Pras)
    dir := NormalizeVector(CreateVector(Pcamera.X, Pcamera.Y, Pcamera.Z))
    ray = CreateRayDifferential(CreatePoint(0,0,0), dir, 0.0, INFINITY, 0.0, 0)
    // Modify ray for depth of field
    if c.lensRadius > 0.0 {
        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        // Compute point on plane of focus
        ft := c.focalDistance / ray.dir.Z
        Pfocus := ray.PointAt(ft)

        // Update ray for effect of lens
        ray.origin = *CreatePoint(lensU, lensV, 0.0)
        ray.dir = *NormalizeVector(Pfocus.Sub(&ray.origin))
    }

    // Compute offset rays for _PerspectiveCamera_ ray differentials
    if c.lensRadius > 0.0 {
        // Compute _PerspectiveCamera_ ray differentials with defocus blur

        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        dx := NormalizeVector(CreateVectorFromPoint(Pcamera).Add(&c.dxCamera))
        ft := c.focalDistance / dx.Z
        pFocus := CreatePoint(0,0,0).Add(dx.Scale(ft))
        ray.rxOrigin = *CreatePoint(lensU, lensV, 0.0)
        ray.rxDirection = *NormalizeVector(pFocus.Sub(&ray.rxOrigin))

        dy := NormalizeVector(CreateVectorFromPoint(Pcamera).Add(&c.dyCamera))
        ft = c.focalDistance / dy.Z
        pFocus = CreatePoint(0,0,0).Add(dy.Scale(ft))
        ray.ryOrigin = *CreatePoint(lensU, lensV, 0.0)
        ray.ryDirection = *NormalizeVector(pFocus.Sub(&ray.ryOrigin))
    } else {
        ray.rxOrigin = ray.origin
        ray.ryOrigin = ray.origin
        ray.rxDirection = *NormalizeVector(CreateVectorFromPoint(Pcamera).Add(&c.dxCamera))
        ray.ryDirection = *NormalizeVector(CreateVectorFromPoint(Pcamera).Add(&c.dyCamera))
    }

    ray.time = sample.time
    ray = RayDifferentialAnimatedTransform(c.cameraToWorld, ray)
    ray.hasDifferentials = true
    
    return ray, 1.0
}

func (c *PerspectiveCamera) CameraToWorld() *AnimatedTransform { return c.cameraToWorld }
func (c *PerspectiveCamera) ShutterOpen() float64 { return c.shutterOpen }
func (c *PerspectiveCamera) ShutterClose() float64 { return c.shutterClose }
func (c *PerspectiveCamera) Film() Film { return c.film }

func NewEnvironmentCamera(cam2world *AnimatedTransform, screenWindow [4]float64, shutteropen,
        shutterclose float64, film Film) *EnvironmentCamera {
	camera := new(EnvironmentCamera)
	
	// camera core members
	camera.cameraToWorld = cam2world
	camera.shutterOpen = shutteropen
	camera.shutterClose = shutterclose
	camera.film = film
    
	return camera
}

func (c *EnvironmentCamera) GenerateRay(sample *Sample) (ray *Ray, weight float64) {
    // Compute environment camera ray direction
    theta := math.Pi * sample.imageY / float64(c.film.YResolution())
    phi := 2 * math.Pi * sample.imageX / float64(c.film.XResolution())
    dir := CreateVector(math.Sin(theta) * math.Cos(phi), math.Cos(theta), math.Sin(theta) * math.Sin(phi))
    ray = CreateRay(CreatePoint(0,0,0), dir, 0.0, INFINITY, sample.time, 0)
    ray = RayAnimatedTransform(c.cameraToWorld, ray)
    
    return ray, 1.0
}

func (c *EnvironmentCamera) GenerateRayDifferential(sample *Sample) (ray *RayDifferential, weight float64) {
	return nil, 0.0
}

func (c *EnvironmentCamera) CameraToWorld() *AnimatedTransform { return c.cameraToWorld }
func (c *EnvironmentCamera) ShutterOpen() float64 { return c.shutterOpen }
func (c *EnvironmentCamera) ShutterClose() float64 { return c.shutterClose }
func (c *EnvironmentCamera) Film() Film { return c.film }

func NewOrthoCamera(cam2world *AnimatedTransform, screenWindow [4]float64, shutteropen,
        shutterclose, lensradius, focaldistance float64, film Film) *OrthoCamera {
	camera := new(OrthoCamera)
	
	// camera core members
	camera.cameraToWorld = cam2world
	camera.shutterOpen = shutteropen
	camera.shutterClose = shutterclose
	camera.film = film
	
	// Compute projective camera screen transformations
	camera.cameraToScreen = OrthographicTransform(0.0, 1.0)
	camera.screenToRaster = ScaleTransform(float64(film.XResolution()), float64(film.YResolution()), 1.0).MultTransform(
		ScaleTransform(1.0/(screenWindow[1]-screenWindow[0]), 1.0/(screenWindow[2]-screenWindow[3]), 1.0)).MultTransform(
		TranslateTransform(&Vector{-screenWindow[0], -screenWindow[3], 0.0}))
	camera.rasterToScreen = InverseTransform(camera.screenToRaster)
	camera.rasterToCamera = InverseTransform(camera.cameraToScreen).MultTransform(camera.rasterToScreen)
	camera.lensRadius = lensradius
	camera.focalDistance = focaldistance
	
    // Compute differential changes in origin for perspective camera rays
    camera.dxCamera = *VectorTransform(camera.rasterToCamera, CreateVector(1,0,0))
    camera.dyCamera = *VectorTransform(camera.rasterToCamera, CreateVector(0,1,0))
    
	return camera
}

func (c *OrthoCamera) GenerateRay(sample *Sample) (ray *Ray, weight float64) { 
    // Generate raster and camera samples
    Pras := CreatePoint(sample.imageX, sample.imageY, 0)
    Pcamera := PointTransform(c.rasterToCamera, Pras)
    ray = CreateRay(Pcamera, CreateVector(0,0,1), 0.0, INFINITY, 0.0, 0)
    // Modify ray for depth of field
    if c.lensRadius > 0.0 {
        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        // Compute point on plane of focus
        ft := c.focalDistance / ray.dir.Z
        Pfocus := ray.PointAt(ft)

        // Update ray for effect of lens
        ray.origin = *CreatePoint(lensU, lensV, 0.0)
        ray.dir = *NormalizeVector(Pfocus.Sub(&ray.origin))
    }
    ray.time = sample.time
    ray = RayAnimatedTransform(c.cameraToWorld, ray)
    
    return ray, 1.0
}
	
func (c *OrthoCamera) GenerateRayDifferential(sample *Sample) (ray *RayDifferential, weight float64) {
    // Compute main orthographic viewing ray

    // Generate raster and camera samples
    Pras := CreatePoint(sample.imageX, sample.imageY, 0)
    Pcamera := PointTransform(c.rasterToCamera, Pras)
    ray = CreateRayDifferential(Pcamera, CreateVector(0,0,1), 0.0, INFINITY, 0.0, 0)

    // Modify ray for depth of field
    if c.lensRadius > 0.0 {
        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        // Compute point on plane of focus
        ft := c.focalDistance / ray.dir.Z
        Pfocus := ray.PointAt(ft)

        // Update ray for effect of lens
        ray.origin = *CreatePoint(lensU, lensV, 0.0)
        ray.dir = *NormalizeVector(Pfocus.Sub(&ray.origin))
    }
    ray.time = sample.time
    // Compute ray differentials for _OrthoCamera_
    if c.lensRadius > 0 {
        // Compute _OrthoCamera_ ray differentials with defocus blur

        // Sample point on lens
        lensU, lensV := ConcentricSampleDisk(sample.lensU, sample.lensV)
        lensU *= c.lensRadius
        lensV *= c.lensRadius

        ft := c.focalDistance / ray.dir.Z

        pFocus := Pcamera.Add(c.dxCamera.Add(CreateVector(0, 0, 1).Scale(ft)))
        ray.rxOrigin = *CreatePoint(lensU, lensV, 0.0)
        ray.rxDirection = *NormalizeVector(pFocus.Sub(&ray.rxOrigin))

        pFocus = Pcamera.Add(c.dyCamera.Add(CreateVector(0, 0, 1).Scale(ft)))
        ray.ryOrigin = *CreatePoint(lensU, lensV, 0.0)
        ray.ryDirection = *NormalizeVector(pFocus.Sub(&ray.ryOrigin))
    } else {
        ray.rxOrigin = *ray.origin.Add(&c.dxCamera)
        ray.ryOrigin = *ray.origin.Add(&c.dyCamera)
        ray.rxDirection = ray.dir
        ray.ryDirection = ray.dir
    }
    ray.hasDifferentials = true
	ray = RayDifferentialAnimatedTransform(c.cameraToWorld, ray)
    
    return ray, 1.0
}

func (c *OrthoCamera) CameraToWorld() *AnimatedTransform { return c.cameraToWorld }
func (c *OrthoCamera) ShutterOpen() float64 { return c.shutterOpen }
func (c *OrthoCamera) ShutterClose() float64 { return c.shutterClose }
func (c *OrthoCamera) Film() Film { return c.film }

func CreatePerspectiveCamera(params *ParamSet, cam2world *AnimatedTransform, film Film) *PerspectiveCamera {
    // Extract common camera parameters from _ParamSet_
    shutteropen := params.FindFloatParam("shutteropen", 0.0)
    shutterclose := params.FindFloatParam("shutterclose", 1.0)
    if shutterclose < shutteropen {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.", shutterclose, shutteropen)
        shutterclose, shutteropen = shutteropen, shutterclose
    }
    lensradius := params.FindFloatParam("lensradius", 0.0)
    focaldistance := params.FindFloatParam("focaldistance", 1.0e30)
    frame := params.FindFloatParam("frameaspectratio", float64(film.XResolution())/float64(film.YResolution()))
    var screen [4]float64
    if frame > 1.0 {
        screen[0] = -frame
        screen[1] =  frame
        screen[2] = -1.0
        screen[3] =  1.0
    } else {
        screen[0] = -1.0
        screen[1] =  1.0
        screen[2] = -1.0 / frame
        screen[3] =  1.0 / frame
    }
    sw := params.FindFloatArrayParam("screenwindow")
    if sw != nil && len(sw) == 4 {
		for i, v := range sw {
			screen[i] = v
        }
	}
    fov := params.FindFloatParam("fov", 90.0)
    halffov := params.FindFloatParam("halffov", -1.0)
    if halffov > 0.0 {
        // hack for structure synth, which exports half of the full fov
        fov = 2.0 * halffov
	}
    return NewPerspectiveCamera(cam2world, screen, shutteropen,
        shutterclose, lensradius, focaldistance, fov, film)
}

func CreateEnvironmentCamera(params *ParamSet, cam2world *AnimatedTransform, film Film) *EnvironmentCamera {
    // Extract common camera parameters from _ParamSet_
    shutteropen := params.FindFloatParam("shutteropen", 0.0)
    shutterclose := params.FindFloatParam("shutterclose", 1.0)
    if shutterclose < shutteropen {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.", shutterclose, shutteropen)
        shutterclose, shutteropen = shutteropen, shutterclose
    }
    frame := params.FindFloatParam("frameaspectratio", float64(film.XResolution())/float64(film.YResolution()))
    var screen [4]float64
    if frame > 1.0 {
        screen[0] = -frame
        screen[1] =  frame
        screen[2] = -1.0
        screen[3] =  1.0
    } else {
        screen[0] = -1.0
        screen[1] =  1.0
        screen[2] = -1.0 / frame
        screen[3] =  1.0 / frame
    }
    sw := params.FindFloatArrayParam("screenwindow")
    if sw != nil && len(sw) == 4 {
		for i, v := range sw {
			screen[i] = v
        }
	}
        
    return NewEnvironmentCamera(cam2world, screen, shutteropen, shutterclose, film)
}

func CreateOrthographicCamera(params *ParamSet, cam2world *AnimatedTransform, film Film) *OrthoCamera {
    // Extract common camera parameters from _ParamSet_
   shutteropen := params.FindFloatParam("shutteropen", 0.0)
    shutterclose := params.FindFloatParam("shutterclose", 1.0)
    if shutterclose < shutteropen {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.", shutterclose, shutteropen)
        shutterclose, shutteropen = shutteropen, shutterclose
    }
    lensradius := params.FindFloatParam("lensradius", 0.0)
    focaldistance := params.FindFloatParam("focaldistance", 1.0e30)
    frame := params.FindFloatParam("frameaspectratio", float64(film.XResolution())/float64(film.YResolution()))
    var screen [4]float64
    if frame > 1.0 {
        screen[0] = -frame
        screen[1] =  frame
        screen[2] = -1.0
        screen[3] =  1.0
    } else {
        screen[0] = -1.0
        screen[1] =  1.0
        screen[2] = -1.0 / frame
        screen[3] =  1.0 / frame
    }
    sw := params.FindFloatArrayParam("screenwindow")
    if sw != nil && len(sw) == 4 {
		for i, v := range sw {
			screen[i] = v
        }
	}
        
    return NewOrthoCamera(cam2world, screen, shutteropen, shutterclose,
        lensradius, focaldistance, film)
}
