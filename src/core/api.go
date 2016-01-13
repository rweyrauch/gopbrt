package api

import (
	"fmt"
	"strings"
    "github.com/rweyrauch/gopbrt/src/shapes"
    "github.com/rweyrauch/gopbrt/src/core"
)

var (
	options *pbrt.Options
)

const (
	MAX_TRANSFORMS       = 2
	START_TRANSFORM_BITS = 1 << 0
	END_TRANSFORM_BITS   = 1 << 1
	ALL_TRANSFORMS_BITS  = ((1 << MAX_TRANSFORMS) - 1)
)

type TransformSet struct {
	t [MAX_TRANSFORMS]pbrt.Transform
}

func (ts *TransformSet) get(i int) pbrt.Transform {
	return ts.t[i]
}
func inverseTransformSet(ts TransformSet) TransformSet {
	var t2 TransformSet
	for i := 0; i < MAX_TRANSFORMS; i++ {
		t2.t[i] = *pbrt.InverseTransform(&ts.t[i])
	}
	return t2
}
func (ts *TransformSet) isAnimated() bool {
	for i := 0; i < MAX_TRANSFORMS-1; i++ {
		if pbrt.NotEqualTransform(&ts.t[i], &ts.t[i+1]) {
			return true
		}
	}
	return false
}

type TransformPair struct {
	first, second *pbrt.Transform
}
type TransformCache struct {
	cache map[pbrt.Transform]TransformPair
	arena pbrt.MemoryArena
}

func (tc *TransformCache) lookup(t pbrt.Transform) (tCache, tInvCache *pbrt.Transform) {
	tpair := tc.cache[t]
	// TODO: handle case of adding a missing transform to cache
	return tpair.first, tpair.second
}
func (tc *TransformCache) clear() {
	// TODO: clear the cache map
}

type RenderOptions struct {
	// RenderOptions Public Data
	transformStartTime, transformEndTime      float64
	FilterName                                string
	FilterParams                              pbrt.ParamSet
	FilmName                                  string
	FilmParams                                pbrt.ParamSet
	SamplerName                               string
	SamplerParams                             pbrt.ParamSet
	AcceleratorName                           string
	AcceleratorParams                         pbrt.ParamSet
	RendererName                              string
	SurfIntegratorName, VolIntegratorName     string
	RendererParams                            pbrt.ParamSet
	SurfIntegratorParams, VolIntegratorParams pbrt.ParamSet
	CameraName                                string
	CameraParams                              pbrt.ParamSet
	CameraToWorld                             TransformSet
	lights                                    []pbrt.Light
	primitives                                []pbrt.Primitive
	volumeRegions                             []pbrt.VolumeRegion
	instances                                 map[string]pbrt.Primitive
	currentInstance                           []pbrt.Primitive
}

func CreateRenderOptions() *RenderOptions {
	opts := new(RenderOptions)
	opts.transformStartTime = 0.0
	opts.transformEndTime = 1.0
	opts.FilterName = "box"
	opts.FilmName = "image"
	opts.SamplerName = "lowdiscrepancy"
	opts.AcceleratorName = "bvh"
	opts.RendererName = "sampler"
	opts.SurfIntegratorName = "directlighting"
	opts.VolIntegratorName = "emission"
	opts.CameraName = "perspective"
	opts.currentInstance = nil
	return opts
}

type GraphicsState struct {
	// Graphics State
	floatTextures        map[string]pbrt.TextureFloat
	spectrumTextures     map[string]pbrt.TextureSpectrum
	materialParams       pbrt.ParamSet
	material             string
	namedMaterials       map[string]pbrt.Material
	currentNamedMaterial string
	areaLightParams      pbrt.ParamSet
	areaLight            string
	reverseOrientation   bool
}

func CreateGraphicsState() *GraphicsState {
	return &GraphicsState{material: "matte", reverseOrientation: false}
}

// API Static Data
const (
	STATE_UNINITIALIZED = 0
	STATE_OPTIONS_BLOCK = 1
	STATE_WORLD_BLOCK   = 2
)

var (
	currentApiState           int = STATE_UNINITIALIZED
	curTransform              TransformSet
	activeTransformBits       int = ALL_TRANSFORMS_BITS
	namedCoordinateSystems    map[string]TransformSet
	renderOptions             *RenderOptions
	graphicsState             GraphicsState
	pushedGraphicsStates      []GraphicsState
	pushedTransforms          []TransformSet
	pushedActiveTransformBits []uint32
	transformCache            TransformCache
)

// Object Creation Function Definitions
func MakeShape(name string, object2world, world2object *pbrt.Transform,
        reverseOrientation bool, paramSet *pbrt.ParamSet) pbrt.Shape {
    var s pbrt.Shape = nil

    if strings.Compare(name, "sphere") == 0 {
        s = shapes.CreateSphereShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "cylinder") == 0 {
        s = shapes.CreateCylinderShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "disk") == 0 {
        s = shapes.CreateDiskShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "cone") == 0 {
        s = shapes.CreateConeShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "paraboloid") == 0 {
        s = shapes.CreateParaboloidShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "hyperboloid") == 0 {
        s = shapes.CreateHyperboloidShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "trianglemesh") == 0 {
        s = shapes.CreateTriangleMeshShape(object2world, world2object, reverseOrientation, paramSet, graphicsState.floatTextures)
    } else if strings.Compare(name, "heightfield") == 0 {
        s = shapes.CreateHeightfieldShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "loopsubdiv") == 0 {
        s = shapes.CreateLoopSubdivShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "nurbs") == 0 {
        s = shapes.CreateNURBSShape(object2world, world2object, reverseOrientation, paramSet)
    } else {
        fmt.Printf("Shape \"%s\" unknown.", name)
    }
    return s
}

// API Function Declarations
func PbrtInit(opt *pbrt.Options) {
    options = opt
    // API Initialization
    if currentApiState != STATE_UNINITIALIZED {
        fmt.Printf("PbrtInit() has already been called.\n")
    }
    currentApiState = STATE_OPTIONS_BLOCK
    renderOptions = CreateRenderOptions()
    graphicsState = CreateGraphicsState()
    //SampledSpectrum::Init()	
}

func PbrtCleanup()                                          {
    //ProbesCleanup()
    // API Cleanup
    if currentApiState == STATE_UNINITIALIZED {
        fmt.Printf("pbrtCleanup() called without pbrtInit().\n")
    } else if currentApiState == STATE_WORLD_BLOCK {
        fmt.Printf("pbrtCleanup() called while inside world block.\n")
    }    
    currentApiState = STATE_UNINITIALIZED
    renderOptions = nil
}

func PbrtIdentity()                                         {}
func PbrtTranslate(dx, dy, dz float64)                      {}
func PbrtRotate(angle, ax, ay, az float64)                  {}
func PbrtScale(sx, sy, sz float64)                          {}
func PbrtLookAt(ex, ey, ez, lx, ly, lz, ux, uy, uz float64) {}
func PbrtConcatTransform(transform pbrt.Matrix4x4)               {}
func PbrtTransform(transform pbrt.Matrix4x4)                     {}
func PbrtCoordinateSystem(name string)                      {}
func PbrtCoordSysTransform(name string)                     {}
func PbrtActiveTransformAll()                               {}
func PbrtActiveTransformEndTime()                           {}
func PbrtActiveTransformStartTime()                         {}

func PbrtTransformTimes(start, end float64) {
	renderOptions.transformStartTime = start
    renderOptions.transformEndTime = end
}

func PbrtPixelFilter(name string, params *pbrt.ParamSet) {
    renderOptions.FilterName = name
    renderOptions.FilterParams = params	
}

func PbrtFilm(filmtype string, params *pbrt.ParamSet) {
    renderOptions.FilmParams = params
    renderOptions.FilmName = filmtype
}

func PbrtSampler(name string, params *pbrt.ParamSet)             {}
func PbrtAccelerator(name string, params *pbrt.ParamSet)         {}
func PbrtSurfaceIntegrator(name string, params *pbrt.ParamSet)   {}
func PbrtVolumeIntegrator(name string, params *pbrt.ParamSet)    {}
func PbrtRenderer(name string, params *pbrt.ParamSet)            {}
func PbrtCamera(camtype string, cameraParams *pbrt.ParamSet)     {}
func PbrtWorldBegin() {
	fmt.Printf("WorldBegin...")
}
func PbrtAttributeBegin()                                                       {}
func PbrtAttributeEnd()                                                         {}
func PbrtTransformBegin()                                                       {}
func PbrtTransformEnd()                                                         {}
func PbrtTexture(name string, textype string, texname string, params *pbrt.ParamSet) {}
func PbrtMaterial(name string, params *pbrt.ParamSet)                                {}
func PbrtMakeNamedMaterial(name string, params *pbrt.ParamSet)                       {}
func PbrtNamedMaterial(name string)                                             {}
func PbrtLightSource(name string, params *pbrt.ParamSet)                             {}
func PbrtAreaLightSource(name string, params *pbrt.ParamSet)                         {}
func PbrtShape(name string, params *pbrt.ParamSet)                                   {}
func PbrtReverseOrientation()                                                   {}
func PbrtVolume(name string, params *pbrt.ParamSet)                                  {}
func PbrtObjectBegin(name string)                                               {}
func PbrtObjectEnd()                                                            {}
func PbrtObjectInstance(name string)                                            {}
func PbrtWorldEnd() {
	fmt.Printf("WorldEnd\n")
}
