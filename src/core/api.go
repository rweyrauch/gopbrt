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
	t [MAX_TRANSFORMS]*Transform
}

func CreateTransformSet() *TransformSet {
	ts := new(TransformSet)
	ts.t[0] = NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4())
	ts.t[1] = NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4())
	return ts
}
func (ts *TransformSet) get(i int) *Transform {
	return ts.t[i]
}
func inverseTransformSet(ts *TransformSet) *TransformSet {
	t2 := new(TransformSet)
	for i := 0; i < MAX_TRANSFORMS; i++ {
		t2.t[i] = InverseTransform(ts.t[i])
	}
	return t2
}
func (ts *TransformSet) isAnimated() bool {
	for i := 0; i < MAX_TRANSFORMS-1; i++ {
		if NotEqualTransform(ts.t[i], ts.t[i+1]) {
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
	arena *MemoryArena
}

func CreateTransformCache() *TransformCache {
	tc := new(TransformCache)
	tc.cache = make(map[Transform]TransformPair, 4)
	tc.arena = nil
	return tc
}
func (tc *TransformCache) Lookup(t *Transform) (tCache, tInvCache *Transform) {
	tpair := tc.cache[*t]
	if tpair.first == nil && tpair.second == nil {
		tinv := InverseTransform(t)
		tpair.first = t
		tpair.second = tinv
		tc.cache[*t] = tpair
	}
	return tpair.first, tpair.second
}
func (tc *TransformCache) Clear() {
	tc.cache = make(map[Transform]TransformPair)
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
	instances                                 map[string][]Primitive
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
	opts.CameraToWorld = CreateTransformSet()
	opts.instances = make(map[string][]Primitive, 4)
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
	gs := &GraphicsState{material: "matte", reverseOrientation: false}
	gs.floatTextures = make(map[string]TextureFloat, 4)
	gs.spectrumTextures = make(map[string]TextureSpectrum, 4)
	gs.namedMaterials = make(map[string]Material, 4)
	gs.materialParams = &ParamSet{nil, nil}
	gs.areaLightParams = &ParamSet{nil, nil}
	return gs
}

func (gs *GraphicsState) CreateMaterial(params *ParamSet) Material {
	mp := CreateTextureParams(params, gs.materialParams, gs.floatTextures, gs.spectrumTextures)

	var mtl Material
	if len(gs.currentNamedMaterial) != 0 {
		mtl = gs.namedMaterials[gs.currentNamedMaterial]
	}
	if mtl == nil {
		mtl = MakeMaterial(gs.material, curTransform.t[0], mp)
	}
	if mtl == nil {
		mtl = MakeMaterial("matte", curTransform.t[0], mp)
	}
	if mtl == nil {
		Severe("Unable to create \"matte\" material?!")
	}
	return mtl
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
	activeTransformBits       uint32 = ALL_TRANSFORMS_BITS
	namedCoordinateSystems    map[string]*TransformSet
	renderOptions             RenderOptions
	graphicsState             GraphicsState
	pushedGraphicsStates      []GraphicsState
	pushedTransforms          []TransformSet
	pushedActiveTransformBits []uint32
	transformCache            TransformCache
)

func init() {
	curTransform = *CreateTransformSet()

	namedCoordinateSystems = make(map[string]*TransformSet, 4)
	pushedGraphicsStates = make([]GraphicsState, 0, 8)
	pushedTransforms = make([]TransformSet, 0, 8)
	pushedActiveTransformBits = make([]uint32, 0, 8)

	transformCache = *CreateTransformCache()
}

// API "Macros"
func VERIFY_INITIALIZED(funcname string) {
	Debug("Trace: %s", funcname)
	if currentApiState == STATE_UNINITIALIZED {
		Error("pbrtInit() must be before calling \"%s()\".  Ignoring.", funcname)
	}
}

func VERIFY_OPTIONS(funcname string) {
	VERIFY_INITIALIZED(funcname)
	if currentApiState == STATE_WORLD_BLOCK {
		Error("Options cannot be set inside world block; \"%s\" not allowed.  Ignoring.", funcname)
	}
}

func VERIFY_WORLD(funcname string) {
	VERIFY_INITIALIZED(funcname)
	if currentApiState == STATE_OPTIONS_BLOCK {
		Error("Scene description must be inside world block; \"%s\" not allowed. Ignoring.", funcname)
	}
}

func FOR_ACTIVE_TRANSFORMS(action func(ndx uint)) {
	var i uint
	for i = 0; i < MAX_TRANSFORMS; i++ {
		if activeTransformBits&(1<<i) != 0 {
			action(i)
		}
	}
}

func WARN_IF_ANIMATED_TRANSFORM(funcname string) {
	if curTransform.isAnimated() {
		Warning("Animated transformations set; ignoring for \"%s\" and using the start transform only", funcname)
	}
}

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
		Warning("Shape \"%s\" unknown.", name)
	}
	return s
}

func MakeMaterial(name string, mtl2world *Transform, mp *TextureParams) Material {
	var material Material = nil
	if strings.Compare(name, "matte") == 0 {
		material = CreateMatteMaterial(mtl2world, mp)
	} else if strings.Compare(name, "plastic") == 0 {
		material = CreatePlasticMaterial(mtl2world, mp)
	} else if strings.Compare(name, "translucent") == 0 {
		material = CreateTranslucentMaterial(mtl2world, mp)
	} else if strings.Compare(name, "glass") == 0 {
		material = CreateGlassMaterial(mtl2world, mp)
	} else if strings.Compare(name, "mirror") == 0 {
		material = CreateMirrorMaterial(mtl2world, mp)
	} else if strings.Compare(name, "mix") == 0 {
		m1 := mp.FindString("namedmaterial1", "")
		m2 := mp.FindString("namedmaterial2", "")
		mat1 := graphicsState.namedMaterials[m1]
		mat2 := graphicsState.namedMaterials[m2]
		if mat1 == nil {
			Error("Named material \"%s\" undefined.  Using \"matte\"", m1)
			mat1 = MakeMaterial("matte", curTransform.t[0], mp)
		}
		if mat2 == nil {
			Error("Named material \"%s\" undefined.  Using \"matte\"", m2)
			mat2 = MakeMaterial("matte", curTransform.t[0], mp)
		}

		material = CreateMixMaterial(mtl2world, mp, mat1, mat2)
	} else if strings.Compare(name, "metal") == 0 {
		material = CreateMetalMaterial(mtl2world, mp)
	} else if strings.Compare(name, "substrate") == 0 {
		material = CreateSubstrateMaterial(mtl2world, mp)
	} else if strings.Compare(name, "uber") == 0 {
		material = CreateUberMaterial(mtl2world, mp)
	} else if strings.Compare(name, "subsurface") == 0 {
		material = CreateSubsurfaceMaterial(mtl2world, mp)
	} else if strings.Compare(name, "kdsubsurface") == 0 {
		material = CreateKdSubsurfaceMaterial(mtl2world, mp)
	} else if strings.Compare(name, "measured") == 0 {
		material = CreateMeasuredMaterial(mtl2world, mp)
	} else if strings.Compare(name, "shinymetal") == 0 {
		material = CreateShinyMetalMaterial(mtl2world, mp)
	} else {
		Warning("Material \"%s\" unknown.", name)
	}
	if material == nil {
		Error("Unable to create material \"%s\"", name)
	}
	return material
}

func MakeFloatTexture(name string, tex2world *Transform, tp *TextureParams) TextureFloat {
	var tex TextureFloat = nil
	if strings.Compare(name, "constant") == 0 {
		tex = CreateConstantFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "scale") == 0 {
		tex = CreateScaleFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "mix") == 0 {
		tex = CreateMixFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "bilerp") == 0 {
		tex = CreateBilerpFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "imagemap") == 0 {
		tex = CreateImageFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "uv") == 0 {
		tex = CreateUVFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "checkerboard") == 0 {
		tex = CreateCheckerboardFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "dots") == 0 {
		tex = CreateDotsFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "fbm") == 0 {
		tex = CreateFBmFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "wrinkled") == 0 {
		tex = CreateWrinkledFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "marble") == 0 {
		tex = CreateMarbleFloatTexture(tex2world, tp)
	} else if strings.Compare(name, "windy") == 0 {
		tex = CreateWindyFloatTexture(tex2world, tp)
	} else {
		Warning("Float texture \"%s\" unknown.", name)
	}
	return tex
}

func MakeSpectrumTexture(name string, tex2world *Transform, tp *TextureParams) TextureSpectrum {
	var tex TextureSpectrum = nil
	if strings.Compare(name, "constant") == 0 {
		tex = CreateConstantSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "scale") == 0 {
		tex = CreateScaleSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "mix") == 0 {
		tex = CreateMixSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "bilerp") == 0 {
		tex = CreateBilerpSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "imagemap") == 0 {
		tex = CreateImageSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "uv") == 0 {
		tex = CreateUVSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "checkerboard") == 0 {
		tex = CreateCheckerboardSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "dots") == 0 {
		tex = CreateDotsSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "fbm") == 0 {
		tex = CreateFBmSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "wrinkled") == 0 {
		tex = CreateWrinkledSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "marble") == 0 {
		tex = CreateMarbleSpectrumTexture(tex2world, tp)
	} else if strings.Compare(name, "windy") == 0 {
		tex = CreateWindySpectrumTexture(tex2world, tp)
	} else {
		Warning("Spectrum texture \"%s\" unknown.", name)
	}
	return tex
}

func MakeLight(name string, light2world *Transform, paramSet *ParamSet) Light {
	var light Light = nil
	if strings.Compare(name, "point") == 0 {
		light = CreatePointLight(light2world, paramSet)
	} else if strings.Compare(name, "spot") == 0 {
		light = CreateSpotLight(light2world, paramSet)
	} else if strings.Compare(name, "goniometric") == 0 {
		light = CreateGoniometricLight(light2world, paramSet)
	} else if strings.Compare(name, "projection") == 0 {
		light = CreateProjectionLight(light2world, paramSet)
	} else if strings.Compare(name, "distant") == 0 {
		light = CreateDistantLight(light2world, paramSet)
	} else if strings.Compare(name, "infinite") == 0 || strings.Compare(name, "exinfinite") == 0 {
		light = CreateInfiniteLight(light2world, paramSet)
	} else {
		Warning("Light \"%s\" unknown.", name)
	}
	return light
}

func MakeAreaLight(name string, light2world *Transform, paramSet *ParamSet, shape Shape) AreaLight {
	var area AreaLight = nil
	if strings.Compare(name, "area") == 0 || strings.Compare(name, "diffuse") == 0 {
		area = CreateDiffuseAreaLight(light2world, paramSet, shape)
	} else {
		Warning("Area light \"%s\" unknown.", name)
	}
	return area
}

func MakeVolumeRegion(name string, volume2world *Transform, paramSet *ParamSet) VolumeRegion {
	var vr VolumeRegion = nil
	if strings.Compare(name, "homogeneous") == 0 {
		vr = CreateHomogeneousVolumeDensityRegion(volume2world, paramSet)
	} else if strings.Compare(name, "volumegrid") == 0 {
		vr = CreateGridVolumeRegion(volume2world, paramSet)
	} else if strings.Compare(name, "exponential") == 0 {
		vr = CreateExponentialVolumeRegion(volume2world, paramSet)
	} else {
		Warning("Volume region \"%s\" unknown.", name)
	}
	return vr
}

func MakeSurfaceIntegrator(name string, paramSet *ParamSet) SurfaceIntegrator {
	var si SurfaceIntegrator = nil
	if strings.Compare(name, "whitted") == 0 {
		si = CreateWhittedSurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "directlighting") == 0 {
		si = CreateDirectLightingIntegrator(paramSet)
	} else if strings.Compare(name, "path") == 0 {
		si = CreatePathSurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "photonmap") == 0 || strings.Compare(name, "exphotonmap") == 0 {
		si = CreatePhotonMapSurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "irradiancecache") == 0 {
		si = CreateIrradianceCacheIntegrator(paramSet)
	} else if strings.Compare(name, "igi") == 0 {
		si = CreateIGISurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "dipolesubsurface") == 0 {
		si = CreateDipoleSubsurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "ambientocclusion") == 0 {
		si = CreateAmbientOcclusionIntegrator(paramSet)
	} else if strings.Compare(name, "useprobes") == 0 {
		si = CreateRadianceProbesSurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "diffuseprt") == 0 {
		si = CreateDiffusePRTIntegratorSurfaceIntegrator(paramSet)
	} else if strings.Compare(name, "glossyprt") == 0 {
		si = CreateGlossyPRTIntegratorSurfaceIntegrator(paramSet)
	} else {
		Warning("Surface integrator \"%s\" unknown.", name)
	}
	return si
}

func MakeVolumeIntegrator(name string, paramSet *ParamSet) VolumeIntegrator {
	var vi VolumeIntegrator = nil
	if strings.Compare(name, "single") == 0 {
		vi = CreateSingleScatteringIntegrator(paramSet)
	} else if strings.Compare(name, "emission") == 0 {
		vi = CreateEmissionVolumeIntegrator(paramSet)
	} else {
		Warning("Volume integrator \"%s\" unknown.", name)
	}
	return vi
}

func MakeAccelerator(name string, prims []Primitive, paramSet *ParamSet) Primitive {
	var accel Primitive = nil
	if strings.Compare(name, "bvh") == 0 {
		accel = CreateBVHAccelerator(prims, paramSet)
	} else if strings.Compare(name, "grid") == 0 {
		accel = CreateGridAccelerator(prims, paramSet)
	} else if strings.Compare(name, "kdtree") == 0 {
		accel = CreateKdTreeAccelerator(prims, paramSet)
	} else if strings.Compare(name, "none") == 0 {
		accel = CreateNoneAccelerator(prims, paramSet)
	} else {
		Warning("Accelerator \"%s\" unknown.", name)
	}
	return accel
}

func MakeCamera(name string, paramSet *ParamSet, cam2worldSet *TransformSet, transformStart, transformEnd float64, film Film) Camera {
	var camera Camera = nil
	var cam2worldStart, cam2worldEnd *Transform
	cam2worldStart, _ = transformCache.Lookup(cam2worldSet.get(0))
	cam2worldEnd, _ = transformCache.Lookup(cam2worldSet.get(1))
	animatedCam2World := NewAnimatedTransform(cam2worldStart, transformStart, cam2worldEnd, transformEnd)
	if strings.Compare(name, "perspective") == 0 {
		camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film)
	} else if strings.Compare(name, "orthographic") == 0 {
		camera = CreateOrthographicCamera(paramSet, animatedCam2World, film)
	} else if strings.Compare(name, "environment") == 0 {
		camera = CreateEnvironmentCamera(paramSet, animatedCam2World, film)
	} else {
		Warning("Camera \"%s\" unknown.", name)
	}
	return camera
}

func MakeSampler(name string, paramSet *ParamSet, film Film, camera Camera) Sampler {
	var sampler Sampler = nil
	if strings.Compare(name, "adaptive") == 0 {
		sampler = CreateAdaptiveSampler(paramSet, film, camera)
	} else if strings.Compare(name, "bestcandidate") == 0 {
		sampler = CreateBestCandidateSampler(paramSet, film, camera)
	} else if strings.Compare(name, "halton") == 0 {
		sampler = CreateHaltonSampler(paramSet, film, camera)
	} else if strings.Compare(name, "lowdiscrepancy") == 0 {
		sampler = CreateLowDiscrepancySampler(paramSet, film, camera)
	} else if strings.Compare(name, "random") == 0 {
		sampler = CreateRandomSampler(paramSet, film, camera)
	} else if strings.Compare(name, "stratified") == 0 {
		sampler = CreateStratifiedSampler(paramSet, film, camera)
	} else {
		Warning("Sampler \"%s\" unknown.", name)
	}
	return sampler
}

func MakeFilter(name string, paramSet *ParamSet) Filter {
	var filter Filter = nil
	if strings.Compare(name, "box") == 0 {
		filter = CreateBoxFilter(paramSet)
	} else if strings.Compare(name, "gaussian") == 0 {
		filter = CreateGaussianFilter(paramSet)
	} else if strings.Compare(name, "mitchell") == 0 {
		filter = CreateMitchellFilter(paramSet)
	} else if strings.Compare(name, "sinc") == 0 {
		filter = CreateLanczosSincFilter(paramSet)
	} else if strings.Compare(name, "triangle") == 0 {
		filter = CreateTriangleFilter(paramSet)
	} else {
		Warning("Filter \"%s\" unknown.", name)
	}
	return filter
}

func MakeFilm(name string, paramSet *ParamSet, filter Filter) Film {
	var film Film = nil
	if strings.Compare(name, "image") == 0 {
		film = CreateImageFilmFromParams(paramSet, filter)
	} else {
		Warning("Film \"%s\" unknown.", name)
	}
	return film
}

// API Function Declarations
func PbrtInit(opt *Options) {
	options = opt
	// API Initialization
	if currentApiState != STATE_UNINITIALIZED {
		Error("pbrtInit() has already been called.")
	}
	currentApiState = STATE_OPTIONS_BLOCK
	renderOptions = *CreateRenderOptions()
	graphicsState = *CreateGraphicsState()
}

func PbrtCleanup() {
	// API Cleanup
	if currentApiState == STATE_UNINITIALIZED {
		Error("pbrtCleanup() called without pbrtInit().")
	} else if currentApiState == STATE_WORLD_BLOCK {
		Error("pbrtCleanup() called while inside world block.")
	}
	currentApiState = STATE_UNINITIALIZED
}

func PbrtIdentity() {
	VERIFY_INITIALIZED("Identity")
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i] = new(Transform) })
}

func PbrtTranslate(dx, dy, dz float64) {
	VERIFY_INITIALIZED("Translate")
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		curTransform.t[i] = curTransform.t[i].MultTransform(TranslateTransform(&Vector{dx, dy, dz}))
	})
}

func PbrtRotate(angle, ax, ay, az float64) {
	VERIFY_INITIALIZED("Rotate")
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		curTransform.t[i] = curTransform.t[i].MultTransform(RotateTransform(angle, &Vector{ax, ay, az}))
	})
}

func PbrtScale(sx, sy, sz float64) {
	VERIFY_INITIALIZED("Scale")
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i] = curTransform.t[i].MultTransform(ScaleTransform(sx, sy, sz)) })
}

func PbrtLookAt(ex, ey, ez, lx, ly, lz, ux, uy, uz float64) {
	VERIFY_INITIALIZED("LookAt")
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		look, _ := LookAtTransform(&Point{ex, ey, ez}, &Point{lx, ly, lz}, &Vector{ux, uy, uz})
		curTransform.t[i] = curTransform.t[i].MultTransform(look)
	})
}

func PbrtConcatTransform(transform Matrix4x4) {
	VERIFY_INITIALIZED("ConcatTransform")
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		xform, _ := NewTransform(&transform)
		curTransform.t[i] = curTransform.t[i].MultTransform(xform)
	})
}

func PbrtTransform(transform Matrix4x4) {
	VERIFY_INITIALIZED("Transform")
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i], _ = NewTransform(&transform) })
}

func PbrtCoordinateSystem(name string) {
	VERIFY_INITIALIZED("CoordinateSystem")
	namedCoordinateSystems[name] = &curTransform
}

func PbrtCoordSysTransform(name string) {
	VERIFY_INITIALIZED("CoordSysTransform")
	if namedCoordinateSystems[name] != nil {
		curTransform = *namedCoordinateSystems[name]
	} else {
		Warning("Could't find named coordinate system \"%s\"", name)
	}
}

func PbrtActiveTransformAll() {
	activeTransformBits = ALL_TRANSFORMS_BITS
}

func PbrtActiveTransformEndTime() {
	activeTransformBits = END_TRANSFORM_BITS
}

func PbrtActiveTransformStartTime() {
	activeTransformBits = START_TRANSFORM_BITS
}

func PbrtTransformTimes(start, end float64) {
	VERIFY_OPTIONS("TransformTimes")
	renderOptions.transformStartTime = start
	renderOptions.transformEndTime = end
}

func PbrtPixelFilter(name string, params *ParamSet) {
	VERIFY_OPTIONS("PixelFilter")
	renderOptions.FilterName = name
	renderOptions.FilterParams = params
}

func PbrtFilm(filmtype string, params *ParamSet) {
	VERIFY_OPTIONS("Film")
	renderOptions.FilmParams = params
	renderOptions.FilmName = filmtype
}

func PbrtSampler(name string, params *ParamSet) {
	VERIFY_OPTIONS("Sampler")
	renderOptions.SamplerName = name
	renderOptions.SamplerParams = params
}

func PbrtAccelerator(name string, params *ParamSet) {
	VERIFY_OPTIONS("Accelerator")
	renderOptions.AcceleratorName = name
	renderOptions.AcceleratorParams = params
}

func PbrtSurfaceIntegrator(name string, params *ParamSet) {
	VERIFY_OPTIONS("SurfaceIntegrator")
	renderOptions.SurfIntegratorName = name
	renderOptions.SurfIntegratorParams = params
}

func PbrtVolumeIntegrator(name string, params *ParamSet) {
	VERIFY_OPTIONS("VolumeIntegrator")
	renderOptions.VolIntegratorName = name
	renderOptions.VolIntegratorParams = params
}

func PbrtRenderer(name string, params *ParamSet) {
	VERIFY_OPTIONS("Renderer")
	renderOptions.RendererName = name
	renderOptions.RendererParams = params
}

func PbrtCamera(camtype string, cameraParams *ParamSet) {
	VERIFY_OPTIONS("Camera")
	renderOptions.CameraName = camtype
	renderOptions.CameraParams = cameraParams
	renderOptions.CameraToWorld = inverseTransformSet(&curTransform)
	namedCoordinateSystems["camera"] = renderOptions.CameraToWorld
}

func PbrtWorldBegin() {
	VERIFY_OPTIONS("WorldBegin")
	currentApiState = STATE_WORLD_BLOCK
	for i := 0; i < MAX_TRANSFORMS; i++ {
		curTransform.t[i] = NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4())
	}
	activeTransformBits = ALL_TRANSFORMS_BITS
	namedCoordinateSystems["world"] = &curTransform
}

func PbrtAttributeBegin() {
	VERIFY_WORLD("AttributeBegin")
	pushedGraphicsStates = append(pushedGraphicsStates, graphicsState)
	pushedTransforms = append(pushedTransforms, curTransform)
	pushedActiveTransformBits = append(pushedActiveTransformBits, activeTransformBits)
}

func PbrtAttributeEnd() {
	VERIFY_WORLD("AttributeEnd")
	if len(pushedGraphicsStates) == 0 {
		Error("Unmatched pbrtAttributeEnd() encountered.  Ignoring it.")
		return
	}

	graphicsState = pushedGraphicsStates[len(pushedGraphicsStates)-1]
	pushedGraphicsStates = pushedGraphicsStates[:len(pushedGraphicsStates)-1]

	curTransform = pushedTransforms[len(pushedTransforms)-1]
	pushedTransforms = pushedTransforms[:len(pushedTransforms)-1]

	activeTransformBits = pushedActiveTransformBits[len(pushedActiveTransformBits)-1]
	pushedActiveTransformBits = pushedActiveTransformBits[:len(pushedActiveTransformBits)-1]
}

func PbrtTransformBegin() {
	VERIFY_WORLD("TransformBegin")
	pushedTransforms = append(pushedTransforms, curTransform)
	pushedActiveTransformBits = append(pushedActiveTransformBits, activeTransformBits)
}

func PbrtTransformEnd() {
	VERIFY_WORLD("TransformEnd")
	if len(pushedTransforms) == 0 {
		Error("Unmatched pbrtTransformEnd() encountered.  Ignoring it.")
		return
	}

	curTransform = pushedTransforms[len(pushedTransforms)-1]
	pushedTransforms = pushedTransforms[:len(pushedTransforms)-1]
	activeTransformBits = pushedActiveTransformBits[len(pushedActiveTransformBits)-1]
	pushedActiveTransformBits = pushedActiveTransformBits[:len(pushedActiveTransformBits)-1]
}

func PbrtTexture(name string, textype string, texname string, params *ParamSet) {
	VERIFY_WORLD("Texture")
	tp := CreateTextureParams(params, params, graphicsState.floatTextures, graphicsState.spectrumTextures)
	if strings.Compare(textype, "float") == 0 {
		// Create _float_ texture and store in _floatTextures_
		if graphicsState.floatTextures[name] != nil {
			Info("Texture \"%s\" being redefined", name)
		}
		WARN_IF_ANIMATED_TRANSFORM("Texture")
		ft := MakeFloatTexture(texname, curTransform.t[0], tp)
		if ft != nil {
			graphicsState.floatTextures[name] = ft
		}
	} else if strings.Compare(textype, "color") == 0 || strings.Compare(textype, "spectrum") == 0 {
		// Create _color_ texture and store in _spectrumTextures_
		if graphicsState.spectrumTextures[name] != nil {
			Info("Texture \"%s\" being redefined", name)
		}
		WARN_IF_ANIMATED_TRANSFORM("Texture")
		st := MakeSpectrumTexture(texname, curTransform.t[0], tp)
		if st != nil {
			graphicsState.spectrumTextures[name] = st
		}
	} else {
		Error("Texture type \"%s\" unknown.", textype)
	}
}

func PbrtMaterial(name string, params *ParamSet) {
	VERIFY_WORLD("Material")
	graphicsState.material = name
	graphicsState.materialParams = params
	graphicsState.currentNamedMaterial = ""
}

func PbrtMakeNamedMaterial(name string, params *ParamSet) {
	VERIFY_WORLD("MakeNamedMaterial")
	// error checking, warning if replace, what to use for transform?
	mp := CreateTextureParams(params, graphicsState.materialParams, graphicsState.floatTextures, graphicsState.spectrumTextures)
	matName := mp.FindString("type", "")
	WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial")
	if len(matName) == 0 {
		Error("No parameter string \"type\" found in MakeNamedMaterial")
	} else {
		mtl := MakeMaterial(matName, curTransform.t[0], mp)
		if mtl != nil {
			graphicsState.namedMaterials[name] = mtl
		}
	}
}

func PbrtNamedMaterial(name string) {
	VERIFY_WORLD("NamedMaterial")
	graphicsState.currentNamedMaterial = name
}

func PbrtLightSource(name string, params *ParamSet) {
	VERIFY_WORLD("LightSource")
	WARN_IF_ANIMATED_TRANSFORM("LightSource")
	lt := MakeLight(name, curTransform.t[0], params)
	if lt == nil {
		Error("pbrtLightSource: light type \"%s\" unknown.", name)
	} else {
		renderOptions.lights = append(renderOptions.lights, lt)
	}
}

func PbrtAreaLightSource(name string, params *ParamSet) {
	VERIFY_WORLD("AreaLightSource")
	graphicsState.areaLight = name
	graphicsState.areaLightParams = params
}

func PbrtShape(name string, params *ParamSet) {
	VERIFY_WORLD("Shape")
	var area AreaLight
	var prim Primitive
	if !curTransform.isAnimated() {
		// Create primitive for static shape
		obj2world, world2obj := transformCache.Lookup(curTransform.t[0])
		shape := MakeShape(name, obj2world, world2obj, graphicsState.reverseOrientation, params)
		if shape == nil {
			return
		}

		mtl := graphicsState.CreateMaterial(params)

		// Possibly create area light for shape
		if len(graphicsState.areaLight) != 0 {
			area = MakeAreaLight(graphicsState.areaLight, curTransform.t[0], graphicsState.areaLightParams, shape)
		}
		prim = NewGeometricPrimitive(shape, mtl, area)
	} else {
		// Create primitive for animated shape

		// Create initial _Shape_ for animated shape
		if len(graphicsState.areaLight) != 0 {
			Warning("Ignoring currently set area light when creating animated shape.")
		}
		identity, _ := transformCache.Lookup(NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4()))
		shape := MakeShape(name, identity, identity, graphicsState.reverseOrientation, params)
		if shape == nil {
			return
		}

		mtl := graphicsState.CreateMaterial(params)

		// Get _animatedWorldToObject_ transform for shape
		Assert(MAX_TRANSFORMS == 2)
		var world2obj [2]*Transform
		_, world2obj[0] = transformCache.Lookup(curTransform.t[0])
		_, world2obj[1] = transformCache.Lookup(curTransform.t[1])
		animatedWorldToObject := NewAnimatedTransform(world2obj[0], renderOptions.transformStartTime, world2obj[1], renderOptions.transformEndTime)

		//Info("AnimW2O: %v", animatedWorldToObject)

		var baseprim Primitive
		baseprim = NewGeometricPrimitive(shape, mtl, nil)
		if !baseprim.CanIntersect() {
			// Refine animated shape and create BVH if more than one shape created
			refinedPrimitives := make([]Primitive, 0)
			refinedPrimitives = baseprim.FullyRefine(refinedPrimitives)
			if len(refinedPrimitives) == 0 {
				return
			}
			if len(refinedPrimitives) > 1 {
				baseprim = NewBVHAccel(refinedPrimitives, 1, "sah")
			} else {
				baseprim = refinedPrimitives[0]
			}
		}
		prim = NewTransformedPrimitive(baseprim, animatedWorldToObject)
	}
	// Add primitive to scene or current instance
	if renderOptions.currentInstance != nil {
		if area != nil {
			Warning("Area lights not supported with object instancing.")
		}
		renderOptions.currentInstance = append(renderOptions.currentInstance, prim)
	} else {
		renderOptions.primitives = append(renderOptions.primitives, prim)
		if area != nil {
			renderOptions.lights = append(renderOptions.lights, area)
		}
	}
}

func PbrtReverseOrientation() {
	VERIFY_WORLD("ReverseOrientation")
	graphicsState.reverseOrientation = !graphicsState.reverseOrientation
}

func PbrtVolume(name string, params *ParamSet) {
	VERIFY_WORLD("Volume")
	WARN_IF_ANIMATED_TRANSFORM("Volume")
	vr := MakeVolumeRegion(name, curTransform.t[0], params)
	if vr != nil {
		renderOptions.volumeRegions = append(renderOptions.volumeRegions, vr)
	}
}

func PbrtObjectBegin(name string) {
	VERIFY_WORLD("ObjectBegin")
	PbrtAttributeBegin()
	if renderOptions.currentInstance != nil {
		Error("ObjectBegin called inside of instance definition")
	}
	renderOptions.instances[name] = make([]Primitive, 0, 1)
	renderOptions.currentInstance = renderOptions.instances[name]
}

func PbrtObjectEnd() {
	VERIFY_WORLD("ObjectEnd")
	if renderOptions.currentInstance == nil {
		Error("ObjectEnd called outside of instance definition")
	}
	renderOptions.currentInstance = nil
	PbrtAttributeEnd()
}

func PbrtObjectInstance(name string) {
	VERIFY_WORLD("ObjectInstance")
	// Object instance error checking
	if renderOptions.currentInstance != nil {
		Error("ObjectInstance can't be called inside instance definition")
		return
	}
	if renderOptions.instances[name] == nil {
		Error("Unable to find instance named \"%s\"", name)
		return
	}
	in := renderOptions.instances[name]
	if len(in) == 0 {
		return
	}
	if len(in) > 1 || !in[0].CanIntersect() {
		// Refine instance _Primitive_s and create aggregate
		accel := MakeAccelerator(renderOptions.AcceleratorName, in, renderOptions.AcceleratorParams)
		if accel == nil {
			accel = MakeAccelerator("bvh", in, &ParamSet{nil, nil})
		}
		if accel == nil {
			Severe("Unable to create \"bvh\" accelerator")
		}
		in = make([]Primitive, 1, 1)
		in[0] = accel
	}
	Assert(MAX_TRANSFORMS == 2)
	var world2instance [2]*Transform
	world2instance[0], _ = transformCache.Lookup(curTransform.t[0])
	world2instance[1], _ = transformCache.Lookup(curTransform.t[1])
	animatedWorldToInstance := NewAnimatedTransform(world2instance[0], renderOptions.transformStartTime,
		world2instance[1], renderOptions.transformEndTime)
	prim := NewTransformedPrimitive(in[0], animatedWorldToInstance)
	renderOptions.primitives = append(renderOptions.primitives, prim)
}

func PbrtWorldEnd() {
	VERIFY_WORLD("WorldEnd")
	// Ensure there are no pushed graphics states
	if len(pushedGraphicsStates) != 0 {
		Warning("Missing end to pbrtAttributeBegin()")
		pushedGraphicsStates = make([]GraphicsState, 0, 8)
	}
	if len(pushedTransforms) != 0 {
		Warning("Missing end to pbrtTransformBegin()")
		pushedTransforms = make([]TransformSet, 0, 8)
		pushedActiveTransformBits = make([]uint32, 0, 8)
	}

	// Create scene and render
	renderer := renderOptions.MakeRenderer()
	scene := renderOptions.MakeScene()
	if scene != nil && renderer != nil {
		renderer.Render(scene)
	}
	renderer = nil
	scene = nil

	// Clean up after rendering
	graphicsState = *CreateGraphicsState()
	transformCache.Clear()
	currentApiState = STATE_OPTIONS_BLOCK
	for i := 0; i < MAX_TRANSFORMS; i++ {
		curTransform.t[i] = NewTransformExplicit(NewIdentityMatrix4x4(), NewIdentityMatrix4x4())
	}
	activeTransformBits = ALL_TRANSFORMS_BITS
	namedCoordinateSystems = make(map[string]*TransformSet)

	ImageTextureFloatClearCache()
	ImageTextureSpectrumClearCache()
}

func (ro *RenderOptions) MakeScene() *Scene {
	// Initialize _volumeRegion_ from volume region(s)
	var volumeRegion VolumeRegion
	if len(ro.volumeRegions) == 0 {
		volumeRegion = nil
	} else if len(ro.volumeRegions) == 1 {
		volumeRegion = ro.volumeRegions[0]
	} else {
		volumeRegion = NewAggregateVolume(ro.volumeRegions)
	}
	accelerator := MakeAccelerator(ro.AcceleratorName, ro.primitives, ro.AcceleratorParams)
	if accelerator == nil {
		accelerator = MakeAccelerator("bvh", ro.primitives, &ParamSet{nil, nil})
	}
	if accelerator == nil {
		Severe("Unable to create \"bvh\" accelerator.")
	}
	scene := NewScene(accelerator, ro.lights, volumeRegion)
	// Erase primitives, lights, and volume regions from _RenderOptions_
	ro.primitives = make([]Primitive, 0, 8)
	ro.lights = make([]Light, 0, 2)
	ro.volumeRegions = make([]VolumeRegion, 0, 2)
	return scene
}

func (ro *RenderOptions) MakeRenderer() Renderer {
	var renderer Renderer
	camera := ro.MakeCamera()
	if strings.Compare(ro.RendererName, "metropolis") == 0 {
		renderer = CreateMetropolisRenderer(ro.RendererParams, camera)
		// Warn if no light sources are defined
		if len(ro.lights) == 0 {
			Warning("No light sources defined in scene; possibly rendering a black image.")
		}
	} else if strings.Compare(ro.RendererName, "createprobes") == 0 { // Create remaining _Renderer_ types
		// Create surface and volume integrators
		surfaceIntegrator := MakeSurfaceIntegrator(ro.SurfIntegratorName, ro.SurfIntegratorParams)
		if surfaceIntegrator == nil {
			Severe("Unable to create surface integrator.")
		}
		volumeIntegrator := MakeVolumeIntegrator(ro.VolIntegratorName, ro.VolIntegratorParams)
		if volumeIntegrator == nil {
			Severe("Unable to create volume integrator.")
		}
		renderer = CreateRadianceProbesRenderer(camera, surfaceIntegrator, volumeIntegrator, ro.RendererParams)
		// Warn if no light sources are defined
		if len(ro.lights) == 0 {
			Warning("No light sources defined in scene; possibly rendering a black image.")
		}
	} else if strings.Compare(ro.RendererName, "aggregatetest") == 0 {
		renderer = CreateAggregateTestRenderer(ro.RendererParams, ro.primitives)
	} else if strings.Compare(ro.RendererName, "testrenderer") == 0 {
		sampler := MakeSampler(ro.SamplerName, ro.SamplerParams, camera.Film(), camera)
		if sampler == nil {
			Severe("Unable to create sampler.")
		}
		// Create surface and volume integrators
		surfaceIntegrator := MakeSurfaceIntegrator(ro.SurfIntegratorName, ro.SurfIntegratorParams)
		if surfaceIntegrator == nil {
			Severe("Unable to create surface integrator.")
		}
		volumeIntegrator := MakeVolumeIntegrator(ro.VolIntegratorName, ro.VolIntegratorParams)
		if volumeIntegrator == nil {
			Severe("Unable to create volume integrator.")
		}
		renderer = &TestRenderer{sampler, surfaceIntegrator, volumeIntegrator, camera}
	} else if strings.Compare(ro.RendererName, "surfacepoints") == 0 {
		pCamera := PointAnimatedTransform(camera.CameraToWorld(), camera.ShutterOpen(), CreatePoint(0, 0, 0))
		renderer = CreateSurfacePointsRenderer(ro.RendererParams, pCamera, camera.ShutterOpen())
	} else {
		if strings.Compare(ro.RendererName, "sampler") != 0 {
			Warning("Renderer type \"%s\" unknown.  Using \"sampler\".", ro.RendererName)
		}
		visIds := ro.RendererParams.FindBoolParam("visualizeobjectids", false)
		sampler := MakeSampler(ro.SamplerName, ro.SamplerParams, camera.Film(), camera)
		if sampler == nil {
			Severe("Unable to create sampler.")
		}
		// Create surface and volume integrators
		surfaceIntegrator := MakeSurfaceIntegrator(ro.SurfIntegratorName, ro.SurfIntegratorParams)
		if surfaceIntegrator == nil {
			Severe("Unable to create surface integrator.")
		}
		volumeIntegrator := MakeVolumeIntegrator(ro.VolIntegratorName, ro.VolIntegratorParams)
		if volumeIntegrator == nil {
			Severe("Unable to create volume integrator.")
		}
		renderer = NewSamplerRenderer(sampler, camera, surfaceIntegrator, volumeIntegrator, visIds)
		// Warn if no light sources are defined
		if len(ro.lights) == 0 {
			Warning("No light sources defined in scene;  possibly rendering a black image.")
		}
	}
	return renderer
}

func (ro *RenderOptions) MakeCamera() Camera {
	filter := MakeFilter(ro.FilterName, ro.FilterParams)
	film := MakeFilm(ro.FilmName, ro.FilmParams, filter)
	if film == nil {
		Severe("Unable to create film.")
	}
	camera := MakeCamera(ro.CameraName, ro.CameraParams, ro.CameraToWorld, ro.transformStartTime, ro.transformEndTime, film)
	if camera == nil {
		Severe("Unable to create camera.")
	}
	return camera
}
