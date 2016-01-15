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
	t [MAX_TRANSFORMS]*Transform
}

func (ts *TransformSet) get(i int) *Transform {
	return ts.t[i]
}
func inverseTransformSet(ts *TransformSet) *TransformSet {
	var t2 *TransformSet
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
	cache map[*Transform]TransformPair
	arena MemoryArena
}

func (tc *TransformCache) lookup(t *Transform) (tCache, tInvCache *Transform) {
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
func (gs *GraphicsState) CreateMaterial(params *ParamSet) Material {
	return nil
}

// API Static Data
const (
	STATE_UNINITIALIZED = 0
	STATE_OPTIONS_BLOCK = 1
	STATE_WORLD_BLOCK   = 2
)

var (
	currentApiState           int = STATE_UNINITIALIZED
	curTransform              *TransformSet
	activeTransformBits       uint32 = ALL_TRANSFORMS_BITS
	namedCoordinateSystems    map[string]*TransformSet
	renderOptions             *RenderOptions
	graphicsState             *GraphicsState
	pushedGraphicsStates      []*GraphicsState
	pushedTransforms          []*TransformSet
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
    	// TODO: implement this
        //string m1 = mp.FindString("namedmaterial1", "")
        //string m2 = mp.FindString("namedmaterial2", "")
        var m1, m2 string
        mat1 := graphicsState.namedMaterials[m1]
        mat2 := graphicsState.namedMaterials[m2]
        if mat1 == nil {
            fmt.Printf("Named material \"%s\" undefined.  Using \"matte\"\n", m1)
            mat1 = MakeMaterial("matte", curTransform.t[0], mp)
        }
        if mat2 == nil {
            fmt.Printf("Named material \"%s\" undefined.  Using \"matte\"\n", m2)
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
        fmt.Printf("Material \"%s\" unknown.\n", name)
    }
    //mp.ReportUnused()
    if material == nil { fmt.Printf("Unable to create material \"%s\"\n", name) }
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
    } else if strings.Compare(name, "marble") == 0{
        tex = CreateMarbleFloatTexture(tex2world, tp)
    } else if strings.Compare(name, "windy") == 0 {
        tex = CreateWindyFloatTexture(tex2world, tp)
    } else {
        fmt.Printf("Float texture \"%s\" unknown.\n", name)
    }
    //tp.ReportUnused()
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
        fmt.Printf("Spectrum texture \"%s\" unknown.\n", name)
	}        
    //tp.ReportUnused()
    return tex;
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
        fmt.Printf("Light \"%s\" unknown.\n", name)
	}        
    //paramSet.ReportUnused()
    return light
}

func MakeAreaLight(name string, object2world, world2object *Transform, shape Shape) AreaLight {
    var area AreaLight = nil
    if strings.Compare(name, "area") == 0 || strings.Compare(name, "diffuse") == 0 {
        area = CreateDiffuseAreaLight(light2world, paramSet, shape);
    } else {
        fmt.Printf("Area light \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused()
    return area
}

func MakeVolumeRegion(name string, volume2world *Transform, paramSet *ParamSet) *VolumeRegion {
    var vr *VolumeRegion = nil
    if strings.Compare(name, "homogeneous") == 0 {
        vr = CreateHomogeneousVolumeDensityRegion(volume2world, paramSet)
    } else if strings.Compare(name, "volumegrid") == 0 {
        vr = CreateGridVolumeRegion(volume2world, paramSet)
    } else if strings.Compare(name, "exponential") == 0 {
        vr = CreateExponentialVolumeRegion(volume2world, paramSet)
    } else {
        fmt.Printf("Volume region \"%s\" unknown.\n", name)
	}        
    //paramSet.ReportUnused()
    return vr
}


func MakeSurfaceIntegrator(name string, paramSet *ParamSet) *SurfaceIntegrator {
    var si *SurfaceIntegrator = nil
    if strings.Compare(name, "whitted") == 0 {
        si = CreateWhittedSurfaceIntegrator(paramSet)
    } else if strings.Compare(name, "directlighting") == 0 {
        si = CreateDirectLightingIntegrator(paramSet)
    } else if strings.Compare(name, "path") == 0 {
        si = CreatePathSurfaceIntegrator(paramSet);
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
        fmt.Printf("Surface integrator \"%s\" unknown.\n", name)
	}
    //paramSet.ReportUnused()
    return si
}


func MakeVolumeIntegrator(name string, paramSet *ParamSet) *VolumeIntegrator {
    var vi *VolumeIntegrator = nil
    if strings.Compare(name, "single") == 0 {
        vi = CreateSingleScatteringIntegrator(paramSet)
    } else if strings.Compare(name, "emission") == 0 {
        vi = CreateEmissionVolumeIntegrator(paramSet)
    } else {
        fmt.Printf("Volume integrator \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused();
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
    } else {
        fmt.Printf("Accelerator \"%s\" unknown.\n", name.c)
    }
    //paramSet.ReportUnused()
    return accel
}


func MakeCamera(name string, paramSet *ParamSet, cam2worldSet *TransformSet, transformStart, transformEnd float64, film Film) Camera {
    var camera Camera = nil
    var cam2worldStart, cam2worldEnd *Transform
    cam2worldStart, _ = transformCache.lookup(cam2worldSet[0])
    cam2worldEnd, _ = transformCache.lookup(cam2worldSet[1])
    animatedCam2World := CreateAnimatedTransform(cam2worldStart, transformStart, cam2worldEnd, transformEnd)
    if strings.Compare(name, "perspective") == 0 {
        camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film)
    } else if strings.Compare(name, "orthographic") == 0 {
        camera = CreateOrthographicCamera(paramSet, animatedCam2World, film)
    } else if strings.Compare(name, "environment") == 0 {
        camera = CreateEnvironmentCamera(paramSet, animatedCam2World, film)
    } else {
        fmt.Printf("Camera \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused()
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
        fmt.Printf("Sampler \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused()
    return sampler
}

func MakeFilter(name string, paramSet *ParamSet) Filter {
    var filter Filter = nil
    if strings.Compare(name, "box") == 0 {
        filter = CreateBoxFilter(paramSet);
    } else if strings.Compare(name, "gaussian") == 0 {
        filter = CreateGaussianFilter(paramSet)
    } else if strings.Compare(name, "mitchell") == 0 {
        filter = CreateMitchellFilter(paramSet)
    } else if strings.Compare(name, "sinc") == 0 {
        filter = CreateSincFilter(paramSet)
    } else if strings.Compare(name, "triangle") == 0 {
        filter = CreateTriangleFilter(paramSet)
    } else {
        fmt.Printf("Filter \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused()
    return filter
}

func MakeFilm(name string, paramSet *ParamSet, filter Filter) Film {
    var film Film = nil
    if strings.Compare(name, "image") == 0 {
        film = CreateImageFilm(paramSet, filter)
    } else {
        fmt.Printf("Film \"%s\" unknown.\n", name)
    }
    //paramSet.ReportUnused()
    return film
}

func FOR_ACTIVE_TRANSFORMS(action func(ndx uint)) {
	var i uint
	for i = 0; i < MAX_TRANSFORMS; i++ {
		if activeTransformBits&(1<<i) != 0 {
			action(i)
		}
	}
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

func PbrtCleanup() {
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

func PbrtIdentity() {
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i] = new(Transform) })
}

func PbrtTranslate(dx, dy, dz float64) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		curTransform.t[i] = curTransform.t[i].MultTransform(TranslateTransform(&Vector{dx, dy, dz}))
	})
}

func PbrtRotate(angle, ax, ay, az float64) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		curTransform.t[i] = curTransform.t[i].MultTransform(RotateTransform(angle, &Vector{ax, ay, az}))
	})
}

func PbrtScale(sx, sy, sz float64) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i] = curTransform.t[i].MultTransform(ScaleTransform(sx, sy, sz)) })
}

func PbrtLookAt(ex, ey, ez, lx, ly, lz, ux, uy, uz float64) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		look, _ := LookAtTransform(&Point{ex, ey, ez}, &Point{lx, ly, lz}, &Vector{ux, uy, uz})
		curTransform.t[i] = curTransform.t[i].MultTransform(look)
	})
}

func PbrtConcatTransform(transform Matrix4x4) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) {
		xform, _ := CreateTransform(&transform)
		curTransform.t[i] = curTransform.t[i].MultTransform(xform)
	})
}

func PbrtTransform(transform Matrix4x4) {
	FOR_ACTIVE_TRANSFORMS(func(i uint) { curTransform.t[i], _ = CreateTransform(&transform) })
}

func PbrtCoordinateSystem(name string) {
	namedCoordinateSystems[name] = curTransform
}

func PbrtCoordSysTransform(name string) {
	if namedCoordinateSystems[name] != nil {
		curTransform = namedCoordinateSystems[name]
	} else {
		fmt.Printf("Could't find named coordinate system \"%s\"\n", name)
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

func PbrtSampler(name string, params *ParamSet) {
	renderOptions.SamplerName = name
	renderOptions.SamplerParams = params
}

func PbrtAccelerator(name string, params *ParamSet) {
	renderOptions.AcceleratorName = name
	renderOptions.AcceleratorParams = params
}

func PbrtSurfaceIntegrator(name string, params *ParamSet) {
	renderOptions.SurfIntegratorName = name
	renderOptions.SurfIntegratorParams = params
}

func PbrtVolumeIntegrator(name string, params *ParamSet) {
	renderOptions.VolIntegratorName = name
	renderOptions.VolIntegratorParams = params
}

func PbrtRenderer(name string, params *ParamSet) {
	renderOptions.RendererName = name
	renderOptions.RendererParams = params
}

func PbrtCamera(camtype string, cameraParams *ParamSet) {
	renderOptions.CameraName = camtype
	renderOptions.CameraParams = cameraParams
	renderOptions.CameraToWorld = inverseTransformSet(curTransform)
	namedCoordinateSystems["camera"] = renderOptions.CameraToWorld
}

func PbrtWorldBegin() {
	currentApiState = STATE_WORLD_BLOCK
	for i := 0; i < MAX_TRANSFORMS; i++ {
		curTransform.t[i] = new(Transform)
	}
	activeTransformBits = ALL_TRANSFORMS_BITS
	namedCoordinateSystems["world"] = curTransform
}

func PbrtAttributeBegin() {
	pushedGraphicsStates = append(pushedGraphicsStates, graphicsState)
	pushedTransforms = append(pushedTransforms, curTransform)
	pushedActiveTransformBits = append(pushedActiveTransformBits, activeTransformBits)
}

func PbrtAttributeEnd() {
	if len(pushedGraphicsStates) == 0 {
		fmt.Printf("Unmatched pbrtAttributeEnd() encountered.  Ignoring it.\n")
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
	pushedTransforms = append(pushedTransforms, curTransform)
	pushedActiveTransformBits = append(pushedActiveTransformBits, activeTransformBits)
}

func PbrtTransformEnd() {
	if len(pushedTransforms) == 0 {
		fmt.Printf("Unmatched pbrtTransformEnd() encountered.  Ignoring it.\n")
		return
	}

	curTransform = pushedTransforms[len(pushedTransforms)-1]
	pushedTransforms = pushedTransforms[:len(pushedTransforms)-1]
	activeTransformBits = pushedActiveTransformBits[len(pushedActiveTransformBits)-1]
	pushedActiveTransformBits = pushedActiveTransformBits[:len(pushedActiveTransformBits)-1]
}

func PbrtTexture(name string, textype string, texname string, params *ParamSet) {}
func PbrtMaterial(name string, params *ParamSet) {
	graphicsState.material = name
	graphicsState.materialParams = params
	graphicsState.currentNamedMaterial = ""
}

func PbrtMakeNamedMaterial(name string, params *ParamSet) {}
func PbrtNamedMaterial(name string) {
	graphicsState.currentNamedMaterial = name
}

func PbrtLightSource(name string, params *ParamSet) {}
func PbrtAreaLightSource(name string, params *ParamSet) {
	graphicsState.areaLight = name
	graphicsState.areaLightParams = params
}

func PbrtShape(name string, params *ParamSet) {
	var area *AreaLight
     if !curTransform.isAnimated() {
        // Create primitive for static shape
        obj2world, world2obj := transformCache.lookup(curTransform.t[0])
        shape := MakeShape(name, obj2world, world2obj, graphicsState.reverseOrientation, params)
        if shape == nil { return }
        
        mtl := graphicsState.CreateMaterial(params)
        //params.ReportUnused()

        // Possibly create area light for shape
        if len(graphicsState.areaLight) != 0 {
            area = MakeAreaLight(graphicsState.areaLight, curTransform.t[0], graphicsState.areaLightParams, shape)
        }
        prim := CreateGeometricPrimitive(shape, mtl, area)
    } else {
        // Create primitive for animated shape

        // Create initial _Shape_ for animated shape
        if len(graphicsState.areaLight) != 0 {
            fmt.Printf("Ignoring currently set area light when creating animated shape.\n")
        }
        identity, _ := transformCache.lookup(Transform())
        shape := MakeShape(name, identity, identity,graphicsState.reverseOrientation, params)
        if shape == nil { return }
        
        mtl := graphicsState.CreateMaterial(params)
        //params.ReportUnused()

        // Get _animatedWorldToObject_ transform for shape
        //Assert(MAX_TRANSFORMS == 2);
        var world2obj [2]*Transform
        _, world2obj[0] = transformCache.lookup(curTransform.t[0])
        _, world2obj[1] = transformCache.lookup(curTransform.t[1])
        animatedWorldToObject := CreateAnimatedTransform(world2obj[0], renderOptions.transformStartTime, world2obj[1], renderOptions.transformEndTime)
        var baseprim Primitive
        baseprim = CreateGeometricPrimitive(shape, mtl, nil)
        if !baseprim.CanIntersect() {
            // Refine animated shape and create BVH if more than one shape created
            refinedPrimitives := baseprim.FullyRefine(nil)
            if len(refinedPrimitives) == 0 { return }
            if len(refinedPrimitives) > 1 {
                baseprim = CreateBVHAccel(refinedPrimitives)
            } else {
                baseprim = refinedPrimitives[0]
            }
        }
        prim := CreateTransformedPrimitive(baseprim, animatedWorldToObject)
    }
    // Add primitive to scene or current instance
    if renderOptions.currentInstance {
        if area != nil {
            fmt.Printf("Area lights not supported with object instancing.\n")
        }
        renderOptions.currentInstance = append(renderOptions.currentInstance, prim)
    } else {
        renderOptions.primitives = append(renderOptions.primitives, prim)
        if (area != nil) {
            renderOptions.lights = append(renderOptions.lights, area)
        }
    }
}

func PbrtReverseOrientation() {
	graphicsState.reverseOrientation = !graphicsState.reverseOrientation
}

func PbrtVolume(name string, params *ParamSet) {}
func PbrtObjectBegin(name string)              {}
func PbrtObjectEnd()                           {}
func PbrtObjectInstance(name string)           {}
func PbrtWorldEnd() {
	fmt.Printf("WorldEnd\n")
}
