package pbrt

var (
	options *Options
)

// TODO: replace these stubs
type TextureFloat struct {
}
type TextureSpectrum struct {
}
type Material struct {
}

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
	FilterParams                              ParamSet
	FilmName                                  string
	FilmParams                                ParamSet
	SamplerName                               string
	SamplerParams                             ParamSet
	AcceleratorName                           string
	AcceleratorParams                         ParamSet
	RendererName                              string
	SurfIntegratorName, VolIntegratorName     string
	RendererParams                            ParamSet
	SurfIntegratorParams, VolIntegratorParams ParamSet
	CameraName                                string
	CameraParams                              ParamSet
	CameraToWorld                             TransformSet
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
	materialParams       ParamSet
	material             string
	namedMaterials       map[string]Material
	currentNamedMaterial string
	areaLightParams      ParamSet
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

// API Function Declarations
func PbrtInit(opt *Options)                                                     {}
func PbrtCleanup()                                                              {}
func PbrtIdentity()                                                             {}
func PbrtTranslate(dx, dy, dz float64)                                          {}
func PbrtRotate(angle, ax, ay, az float64)                                      {}
func PbrtScale(sx, sy, sz float64)                                              {}
func PbrtLookAt(ex, ey, ez, lx, ly, lz, ux, uy, uz float64)                     {}
func PbrtConcatTransform(transform [16]float64)                                 {}
func PbrtTransform(transform [16]float64)                                       {}
func PbrtCoordinateSystem(name string)                                          {}
func PbrtCoordSysTransform(name string)                                         {}
func PbrtActiveTransformAll()                                                   {}
func PbrtActiveTransformEndTime()                                               {}
func PbrtActiveTransformStartTime()                                             {}
func PbrtTransformTimes(start, end float64)                                     {}
func PbrtPixelFilter(name string, params *ParamSet)                             {}
func PbrtFilm(filmtype string, params *ParamSet)                                {}
func PbrtSampler(name string, params *ParamSet)                                 {}
func PbrtAccelerator(name string, params *ParamSet)                             {}
func PbrtSurfaceIntegrator(name string, params *ParamSet)                       {}
func PbrtVolumeIntegrator(name string, params *ParamSet)                        {}
func PbrtRenderer(name string, params *ParamSet)                                {}
func PbrtCamera(camtype string, cameraParams *ParamSet)                         {}
func PbrtWorldBegin()                                                           {}
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
func PbrtWorldEnd()                                                             {}
