package pbrt

import (
	"fmt"
	"strings"
)

var (
	options *Options
)

const (
	MAX_TRANSFORMS       = 2
	START_TRANSFORM_BITS = 1 << 0
	END_TRANSFORM_BITS   = 1 << 1
	ALL_TRANSFORMS_BITS  = ((1 << MAX_TRANSFORMS) - 1)
)

type TransformSet struct {
	t [MAX_TRANSFORMS]Transform
}

func (ts *TransformSet) get(i int) Transform {
	return ts.t[i]
}
func inverseTransformSet(ts TransformSet) TransformSet {
	var t2 TransformSet
	for i := 0; i < MAX_TRANSFORMS; i++ {
		t2.t[i] = *InverseTransform(&ts.t[i])
	}
	return t2
}
func (ts *TransformSet) isAnimated() bool {
	for i := 0; i < MAX_TRANSFORMS-1; i++ {
		if NotEqualTransform(&ts.t[i], &ts.t[i+1]) {
			return true
		}
	}
	return false
}

type TransformPair struct {
	first, second *Transform
}
type TransformCache struct {
	cache map[Transform]TransformPair
	arena MemoryArena
}

func (tc *TransformCache) lookup(t Transform) (tCache, tInvCache *Transform) {
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
	FilterParams                              *ParamSet
	FilmName                                  string
	FilmParams                                *ParamSet
	SamplerName                               string
	SamplerParams                             *ParamSet
	AcceleratorName                           string
	AcceleratorParams                         *ParamSet
	RendererName                              string
	SurfIntegratorName, VolIntegratorName     string
	RendererParams                            *ParamSet
	SurfIntegratorParams, VolIntegratorParams *ParamSet
	CameraName                                string
	CameraParams                              *ParamSet
	CameraToWorld                             *TransformSet
	lights                                    []Light
	primitives                                []Primitive
	volumeRegions                             []VolumeRegion
	instances                                 map[string]Primitive
	currentInstance                           []Primitive
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
	floatTextures        map[string]TextureFloat
	spectrumTextures     map[string]TextureSpectrum
	materialParams       *ParamSet
	material             string
	namedMaterials       map[string]Material
	currentNamedMaterial string
	areaLightParams      *ParamSet
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
	graphicsState             *GraphicsState
	pushedGraphicsStates      []*GraphicsState
	pushedTransforms          []TransformSet
	pushedActiveTransformBits []uint32
	transformCache            TransformCache
)

// Object Creation Function Definitions
func MakeShape(name string, object2world, world2object *Transform,
        reverseOrientation bool, paramSet *ParamSet) Shape {
    var s Shape = nil

    if strings.Compare(name, "sphere") == 0 {
        s = CreateSphereShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "cylinder") == 0 {
        s = CreateCylinderShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "disk") == 0 {
        s = CreateDiskShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "cone") == 0 {
        s = CreateConeShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "paraboloid") == 0 {
        s = CreateParaboloidShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "hyperboloid") == 0 {
        s = CreateHyperboloidShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "trianglemesh") == 0 {
        s = CreateTriangleMeshShape(object2world, world2object, reverseOrientation, paramSet, graphicsState.floatTextures)
    } else if strings.Compare(name, "heightfield") == 0 {
        s = CreateHeightfieldShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "loopsubdiv") == 0 {
        s = CreateLoopSubdivShape(object2world, world2object, reverseOrientation, paramSet)
    } else if strings.Compare(name, "nurbs") == 0 {
        s = CreateNURBSShape(object2world, world2object, reverseOrientation, paramSet)
    } else {
        fmt.Printf("Shape \"%s\" unknown.", name)
    }
    return s
}

// API Function Declarations
func PbrtInit(opt *Options) {
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
func PbrtConcatTransform(transform Matrix4x4)               {}
func PbrtTransform(transform Matrix4x4)                     {}
func PbrtCoordinateSystem(name string)                      {}
func PbrtCoordSysTransform(name string)                     {}
func PbrtActiveTransformAll()                               {}
func PbrtActiveTransformEndTime()                           {}
func PbrtActiveTransformStartTime()                         {}

func PbrtTransformTimes(start, end float64) {
	renderOptions.transformStartTime = start
    renderOptions.transformEndTime = end
}

func PbrtPixelFilter(name string, params *ParamSet) {
    renderOptions.FilterName = name
    renderOptions.FilterParams = params	
}

func PbrtFilm(filmtype string, params *ParamSet) {
    renderOptions.FilmParams = params
    renderOptions.FilmName = filmtype
}

func PbrtSampler(name string, params *ParamSet)             {}
func PbrtAccelerator(name string, params *ParamSet)         {}
func PbrtSurfaceIntegrator(name string, params *ParamSet)   {}
func PbrtVolumeIntegrator(name string, params *ParamSet)    {}
func PbrtRenderer(name string, params *ParamSet)            {}
func PbrtCamera(camtype string, cameraParams *ParamSet)     {}
func PbrtWorldBegin() {
	fmt.Printf("WorldBegin...")
}
func PbrtAttributeBegin()                                                       {}
func PbrtAttributeEnd()                                                         {}
func PbrtTransformBegin()                                                       {}
func PbrtTransformEnd()                                                         {}
func PbrtTexture(name string, textype string, texname string, params *ParamSet) {}
func PbrtMaterial(name string, params *ParamSet)                                {}
func PbrtMakeNamedMaterial(name string, params *ParamSet)                       {}
func PbrtNamedMaterial(name string)                                             {}
func PbrtLightSource(name string, params *ParamSet)                             {}
func PbrtAreaLightSource(name string, params *ParamSet)                         {}
func PbrtShape(name string, params *ParamSet)                                   {}
func PbrtReverseOrientation()                                                   {}
func PbrtVolume(name string, params *ParamSet)                                  {}
func PbrtObjectBegin(name string)                                               {}
func PbrtObjectEnd()                                                            {}
func PbrtObjectInstance(name string)                                            {}
func PbrtWorldEnd() {
	fmt.Printf("WorldEnd\n")
}
