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

type Scene struct {
	aggregate    Primitive
	lights       []Light
	volumeRegion VolumeRegion
	bound        BBox
}

func NewScene(accel Primitive, lights []Light, vr VolumeRegion) *Scene {
	scene := new(Scene)
	scene.aggregate = accel
	scene.lights = lights
	scene.volumeRegion = vr
	scene.bound = *accel.WorldBound()
	if scene.volumeRegion != nil {
		scene.bound = *UnionBBoxes(&scene.bound, scene.volumeRegion.WorldBound())
	}
	return scene
}
func (s *Scene) Intersect(ray RayBase) (hit bool, isect *Intersection) {
	return s.aggregate.Intersect(ray)
}

func (s *Scene) IntersectP(ray RayBase) bool {
	return s.aggregate.IntersectP(ray)
}

func (s *Scene) WorldBound() *BBox { return &s.bound }
