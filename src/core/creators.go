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

/*
	Shapes
*/

func CreateConeShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cone {
	radius := params.FindFloatParam("radius", 1.0)
	height := params.FindFloatParam("height", radius)
	phimax := params.FindFloatParam("phimax", 360)
	return NewCone(o2w, w2o, reverseOrientation, radius, height, phimax)
}

func CreateCylinderShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Cylinder {
	radius := params.FindFloatParam("radius", 1.0)
	zmin := params.FindFloatParam("zmin", -1.0)
	zmax := params.FindFloatParam("zmax", 1.0)
	phimax := params.FindFloatParam("phimax", 360)
	return NewCylinder(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}

func CreateDiskShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Disk {
	height := params.FindFloatParam("height", 0.0)
	radius := params.FindFloatParam("radius", 1.0)
	innerradius := params.FindFloatParam("innerradius", 0.0)
	phimax := params.FindFloatParam("phimax", 360)
	return NewDisk(o2w, w2o, reverseOrientation, height, radius, innerradius, phimax)
}

func CreateHeightfieldShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Heightfield {
	nu := params.FindIntParam("nu", -1)
	nv := params.FindIntParam("nv", -1)
	Pz := params.FindFloatArrayParam("Pz")
	Assert(nu != -1 && nv != -1 && Pz != nil)
	Assert(len(Pz) == nu*nv)
	return NewHeightfield(o2w, w2o, reverseOrientation, nu, nv, Pz)
}

func CreateHyperboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Hyperboloid {
	p1 := params.FindPointParam("p1", Point{0, 0, 0})
	p2 := params.FindPointParam("p2", Point{1, 1, 1})
	phimax := params.FindFloatParam("phimax", 360)
	return NewHyperboloid(o2w, w2o, reverseOrientation, p1, p2, phimax)
}

func CreateLoopSubdivShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *LoopSubdiv {
	nlevels := params.FindIntParam("nlevels", 3)
	vi := params.FindIntArrayParam("indices")
	P := params.FindPointArrayParam("P")
	if vi == nil || P == nil {
		return nil
	}

	// don't actually use this for now...
	//scheme := params.FindStringParam("scheme", "loop")

	return NewLoopSubdiv(o2w, w2o, reverseOrientation, len(vi)/3, len(P),
		vi, P, nlevels)
}

func CreateNURBSShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *NURBS {
	nu := params.FindIntParam("nu", -1)
	uorder := params.FindIntParam("uorder", -1)
	uknots := params.FindFloatArrayParam("uknots")
	Assert(nu != -1 && uorder != -1 && uknots != nil)
	Assert(len(uknots) == nu+uorder)
	u0 := params.FindFloatParam("u0", uknots[uorder-1])
	u1 := params.FindFloatParam("u1", uknots[nu])

	nv := params.FindIntParam("nv", -1)
	vorder := params.FindIntParam("vorder", -1)
	vknots := params.FindFloatArrayParam("vknots")
	Assert(nv != -1 && vorder != -1 && vknots != nil)
	Assert(len(vknots) == nv+vorder)
	v0 := params.FindFloatParam("v0", vknots[vorder-1])
	v1 := params.FindFloatParam("v1", vknots[nv])

	isHomogeneous := false
	npts := 0
	var P []float64
	Pnts := params.FindPointArrayParam("P")
	if Pnts == nil {
		P = params.FindFloatArrayParam("Pw")
		if P == nil {
			Error("Must provide control points via \"P\" or \"Pw\" parameter to NURBS shape.")
			return nil
		}
		if len(P)%4 != 0 {
			Error("Number of \"Pw\" control points provided to NURBS shape must be multiple of four")
			return nil
		}

		npts = len(P) / 4
		isHomogeneous = true
	} else {
		for _, p := range Pnts {
			P = append(P, p.X, p.Y, p.Z)
		}
		npts = len(Pnts)
	}
	if npts != nu*nv {
		Error("NURBS shape was expecting %dx%d=%d control points, was given %d", nu, nv, nu*nv, npts)
		return nil
	}

	return NewNURBS(o2w, w2o, reverseOrientation, nu, uorder, uknots, u0, u1,
		nv, vorder, vknots, v0, v1, P,
		isHomogeneous)
}

func CreateParaboloidShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Paraboloid {
	radius := params.FindFloatParam("radius", 1)
	zmin := params.FindFloatParam("zmin", 0)
	zmax := params.FindFloatParam("zmax", 1)
	phimax := params.FindFloatParam("phimax", 360)
	return NewParaboloid(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}

func CreateSphereShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Sphere {
	radius := params.FindFloatParam("radius", 1.0)
	zmin := params.FindFloatParam("zmin", -radius)
	zmax := params.FindFloatParam("zmax", radius)
	phimax := params.FindFloatParam("phimax", 360)
	return CreateSphere(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax)
}

func CreateTriangleMeshShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet, textures map[string]TextureFloat) *TriangleMesh {
	vi := params.FindIntArrayParam("indices")
	P := params.FindPointArrayParam("P")
	uvs := params.FindFloatArrayParam("uv")
	if uvs == nil {
		uvs = params.FindFloatArrayParam("st")
	}
	discardDegnerateUVs := params.FindBoolParam("discarddegenerateUVs", false)

	// XXX should complain if uvs aren't an array of 2...
	if uvs != nil {
		if len(uvs) < 2*len(P) {
			Error("Not enough of \"uv\"s for triangle mesh.  Expencted %d, found %d.  Discarding.", 2*len(P), len(uvs))
			uvs = nil
		} else if len(uvs) > 2*len(P) {
			Warning("More \"uv\"s provided than will be used for triangle  mesh.  (%d expcted, %d found)", 2*len(P), len(uvs))
		}
	}
	if vi == nil || P == nil {
		return nil
	}

	S := params.FindVectorArrayParam("S")
	if S != nil && len(S) != len(P) {
		Error("Number of \"S\"s for triangle mesh must match \"P\"s")
		S = nil
	}
	N := params.FindNormalArrayParam("N")
	if N != nil && len(N) != len(P) {
		Error("Number of \"N\"s for triangle mesh must match \"P\"s")
		N = nil
	}
	if discardDegnerateUVs && uvs != nil && N != nil {
		// if there are normals, check for bad uv's that
		// give degenerate mappings; discard them if so
		for i := 0; i < len(vi); i = i + 3 {
			area := 0.5 * CrossVector(P[vi[i+0]].Sub(&P[vi[i+1]]), P[vi[i+2]].Sub(&P[vi[i+1]])).Length()
			if area < 1.0e-7 {
				continue
			} // ignore degenerate tris.
			if (uvs[2*vi[i+0]] == uvs[2*vi[i+1]] &&
				uvs[2*vi[i+0]+1] == uvs[2*vi[i+1]+1]) ||
				(uvs[2*vi[i+1]] == uvs[2*vi[i+2]] &&
					uvs[2*vi[i+1]+1] == uvs[2*vi[i+2]+1]) ||
				(uvs[2*vi[i+2]] == uvs[2*vi[i+0]] &&
					uvs[2*vi[i+2]+1] == uvs[2*vi[i+0]+1]) {
				Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.")
				uvs = nil
				break
			}
		}
	}

	for i := 0; i < len(vi); i++ {
		if vi[i] >= len(P) {
			Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
				vi[i], len(P))
			return nil
		}
	}

	var alphaTex TextureFloat
	alphaTexName := params.FindTextureParam("alpha")
	if len(alphaTexName) != 0 {
		alphaTex = textures[alphaTexName]
		if alphaTex == nil {
			Error("Couldn't find float texture \"%s\" for \"alpha\" parameter", alphaTexName)
		}
	} else if params.FindFloatParam("alpha", 1.0) == 0.0 {
		alphaTex = &ConstantTextureFloat{0.0}
	}

	return NewTriangleMesh(o2w, w2o, reverseOrientation, vi, P,
		N, S, uvs, alphaTex)
}

/*
	Materials
*/

func CreateGlassMaterial(xform *Transform, mp *TextureParams) *GlassMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	Kt := mp.GetSpectrumTexture("Kt", *NewSpectrum1(1.0))
	index := mp.GetFloatTexture("index", 1.5)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &GlassMaterial{Kr, Kt, index, bumpMap}
}

func CreateKdSubsurfaceMaterial(xform *Transform, mp *TextureParams) *KdSubsurfaceMaterial {
	kd := mp.GetSpectrumTexture("Kd", *NewSpectrumRGB(0.5, 0.5, 0.5))
	mfp := mp.GetFloatTexture("meanfreepath", 1.0)
	ior := mp.GetFloatTexture("index", 1.3)
	kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &KdSubsurfaceMaterial{kd, kr, mfp, ior, bumpMap}
}

func CreateMatteMaterial(xform *Transform, mp *TextureParams) *MatteMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.5))
	sigma := mp.GetFloatTexture("sigma", 0.0)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &MatteMaterial{Kd, sigma, bumpMap}
}

func CreateMirrorMaterial(xform *Transform, mp *TextureParams) *MirrorMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(0.9))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &MirrorMaterial{Kr, bumpMap}
}

func CreateMixMaterial(xform *Transform, mp *TextureParams, m1, m2 Material) *MixMaterial {
	scale := mp.GetSpectrumTexture("amount", *NewSpectrum1(0.5))
	return &MixMaterial{m1, m2, scale}
}

func CreatePlasticMaterial(xform *Transform, mp *TextureParams) *PlasticMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &PlasticMaterial{Kd, Ks, roughness, bumpMap}
}

func CreateShinyMetalMaterial(xform *Transform, mp *TextureParams) *ShinyMetalMaterial {
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(1.0))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &ShinyMetalMaterial{Kr, Ks, roughness, bumpMap}
}

func CreateSubstrateMaterial(xform *Transform, mp *TextureParams) *SubstrateMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.5))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.5))
	uroughness := mp.GetFloatTexture("uroughness", 0.1)
	vroughness := mp.GetFloatTexture("vroughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &SubstrateMaterial{Kd, Ks, uroughness, vroughness, bumpMap}
}

func CreateSubsurfaceMaterial(xform *Transform, mp *TextureParams) *SubsurfaceMaterial {
	sa := NewSpectrumRGB(0.0011, 0.0024, 0.014)
	sps := NewSpectrumRGB(2.55, 3.21, 3.77)
	name := mp.FindString("name", "")
	if len(name) != 0 {
		found, sap, spsp := GetVolumeScatteringProperties(name)
		if !found {
			Warning("Named material \"%s\" not found.  Using defaults.", name)
		} else {
			sa = sap
			sps = spsp
		}
	}
	scale := mp.FindFloat("scale", 1.0)

	sigma_a := mp.GetSpectrumTexture("sigma_a", *sa)
	sigma_prime_s := mp.GetSpectrumTexture("sigma_prime_s", *sps)
	ior := mp.GetFloatTexture("index", 1.3)
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &SubsurfaceMaterial{scale, Kr, sigma_a, sigma_prime_s, ior, bumpMap}
}

func CreateTranslucentMaterial(xform *Transform, mp *TextureParams) *TranslucentMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	reflect := mp.GetSpectrumTexture("reflect", *NewSpectrum1(0.5))
	transmit := mp.GetSpectrumTexture("transmit", *NewSpectrum1(0.5))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &TranslucentMaterial{Kd, Ks, reflect, transmit, roughness, bumpMap}
}

func CreateUberMaterial(xform *Transform, mp *TextureParams) *UberMaterial {
	Kd := mp.GetSpectrumTexture("Kd", *NewSpectrum1(0.25))
	Ks := mp.GetSpectrumTexture("Ks", *NewSpectrum1(0.25))
	Kr := mp.GetSpectrumTexture("Kr", *NewSpectrum1(0.0))
	Kt := mp.GetSpectrumTexture("Kt", *NewSpectrum1(0.0))
	roughness := mp.GetFloatTexture("roughness", 0.1)
	eta := mp.GetFloatTexture("index", 1.5)
	opacity := mp.GetSpectrumTexture("opacity", *NewSpectrum1(1.0))
	bumpMap := mp.GetFloatTextureOrNil("bumpmap")
	return &UberMaterial{Kd, Ks, Kr, Kt, opacity, roughness, eta, bumpMap}
}

/*
	Textures
*/

func CreateBilerpFloatTexture(tex2world *Transform, tp *TextureParams) *BilerpTextureFloat {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1, 1, 0, 0)
	}
	return NewBilerpTextureFloat(mapping,
		tp.FindFloat("v00", 0.0), tp.FindFloat("v01", 1.0),
		tp.FindFloat("v10", 0.0), tp.FindFloat("v11", 1.0))
}

func CreateBilerpSpectrumTexture(tex2world *Transform, tp *TextureParams) *BilerpTextureSpectrum {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1, 1, 0, 0)
	}
	return NewBilerpTextureSpectrum(mapping,
		tp.FindSpectrum("v00", *NewSpectrum1(0)), tp.FindSpectrum("v01", *NewSpectrum1(1)),
		tp.FindSpectrum("v10", *NewSpectrum1(0)), tp.FindSpectrum("v11", *NewSpectrum1(1)))
}

func CreateCheckerboardFloatTexture(tex2world *Transform, tp *TextureParams) TextureFloat {
	dim := tp.FindInt("dimension", 2)
	if dim != 2 && dim != 3 {
		Error("%d dimensional checkerboard texture not supported", dim)
		return nil
	}
	tex1 := tp.GetFloatTexture("tex1", 1.0)
	tex2 := tp.GetFloatTexture("tex2", 0.0)
	if dim == 2 {
		// Initialize 2D texture mapping _map_ from _tp_
		var mapping TextureMapping2D
		maptype := tp.FindString("mapping", "uv")
		if strings.Compare(maptype, "uv") == 0 {
			su := tp.FindFloat("uscale", 1.0)
			sv := tp.FindFloat("vscale", 1.0)
			du := tp.FindFloat("udelta", 0.0)
			dv := tp.FindFloat("vdelta", 0.0)
			mapping = NewUVMapping2D(su, sv, du, dv)
		} else if strings.Compare(maptype, "spherical") == 0 {
			mapping = NewSphericalMapping2D(InverseTransform(tex2world))
		} else if strings.Compare(maptype, "cylindrical") == 0 {
			mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
		} else if strings.Compare(maptype, "planar") == 0 {
			mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
				tp.FindVector("v2", Vector{0, 1, 0}),
				tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
		} else {
			Error("2D texture mapping \"%s\" unknown", maptype)
			mapping = NewUVMapping2D(1.0, 1.0, 0.0, 0.0)
		}
		aamode := tp.FindString("aamode", "closedform")
		return NewCheckerboard2DTextureFloat(mapping, tex1, tex2, aamode)
	} else {
		// Initialize 3D texture mapping _map_ from _tp_
		mapping := NewIdentityMapping3D(tex2world)
		return NewCheckerboard3DTextureFloat(mapping, tex1, tex2)
	}
}

func CreateCheckerboardSpectrumTexture(tex2world *Transform, tp *TextureParams) TextureSpectrum {
	dim := tp.FindInt("dimension", 2)
	if dim != 2 && dim != 3 {
		Error("%d dimensional checkerboard texture not supported", dim)
		return nil
	}
	tex1 := tp.GetSpectrumTexture("tex1", *NewSpectrum1(1.0))
	tex2 := tp.GetSpectrumTexture("tex2", *NewSpectrum1(0.0))
	if dim == 2 {
		// Initialize 2D texture mapping _map_ from _tp_
		var mapping TextureMapping2D
		maptype := tp.FindString("mapping", "uv")
		if strings.Compare(maptype, "uv") == 0 {
			su := tp.FindFloat("uscale", 1.0)
			sv := tp.FindFloat("vscale", 1.0)
			du := tp.FindFloat("udelta", 0.0)
			dv := tp.FindFloat("vdelta", 0.0)
			mapping = NewUVMapping2D(su, sv, du, dv)
		} else if strings.Compare(maptype, "spherical") == 0 {
			mapping = NewSphericalMapping2D(InverseTransform(tex2world))
		} else if strings.Compare(maptype, "cylindrical") == 0 {
			mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
		} else if strings.Compare(maptype, "planar") == 0 {
			mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
				tp.FindVector("v2", Vector{0, 1, 0}),
				tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
		} else {
			Error("2D texture mapping \"%s\" unknown", maptype)
			mapping = NewUVMapping2D(1.0, 1.0, 0.0, 0.0)
		}
		aamode := tp.FindString("aamode", "closedform")
		return NewCheckerboard2DTextureSpectrum(mapping, tex1, tex2, aamode)
	} else {
		// Initialize 3D texture mapping _map_ from _tp_
		mapping := NewIdentityMapping3D(tex2world)
		return NewCheckerboard3DTextureSpectrum(mapping, tex1, tex2)
	}
}

func CreateConstantFloatTexture(tex2world *Transform, tp *TextureParams) *ConstantTextureFloat {
	return &ConstantTextureFloat{tp.FindFloat("value", 1.0)}
}
func CreateConstantSpectrumTexture(tex2world *Transform, tp *TextureParams) *ConstantTextureSpectrum {
	return &ConstantTextureSpectrum{tp.FindSpectrum("value", *NewSpectrum1(1.0))}
}

func CreateDotsFloatTexture(tex2world *Transform, tp *TextureParams) *DotsTextureFloat {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1, 1, 0, 0)
	}
	return NewDotsTextureFloat(mapping, tp.GetFloatTexture("inside", 1.0), tp.GetFloatTexture("outside", 0.0))
}

func CreateDotsSpectrumTexture(tex2world *Transform, tp *TextureParams) *DotsTextureSpectrum {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1, 1, 0, 0)
	}
	return NewDotsTextureSpectrum(mapping, tp.GetSpectrumTexture("inside", *NewSpectrum1(1.0)), tp.GetSpectrumTexture("outside", *NewSpectrum1(0.0)))
}

func CreateFBmFloatTexture(tex2world *Transform, tp *TextureParams) *FBmTextureFloat {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return NewFBmTextureFloat(tp.FindInt("octaves", 8), tp.FindFloat("roughness", 0.5), mapping)
}
func CreateFBmSpectrumTexture(tex2world *Transform, tp *TextureParams) *FBmTextureSpectrum {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return NewFBmTextureSpectrum(tp.FindInt("octaves", 8), tp.FindFloat("roughness", 0.5), mapping)
}

func CreateImageFloatTexture(tex2world *Transform, tp *TextureParams) *ImageTextureFloat {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1.0, 1.0, 0.0, 0.0)
	}

	// Initialize _ImageTexture_ parameters
	maxAniso := tp.FindFloat("maxanisotropy", 8.0)
	trilerp := tp.FindBool("trilinear", false)
	wrap := tp.FindString("wrap", "repeat")
	wrapMode := TEXTURE_REPEAT
	if strings.Compare(wrap, "black") == 0 {
		wrapMode = TEXTURE_BLACK
	} else if strings.Compare(wrap, "clamp") == 0 {
		wrapMode = TEXTURE_CLAMP
	}
	scale := tp.FindFloat("scale", 1.0)
	gamma := tp.FindFloat("gamma", 1.0)
	return NewImageTextureFloat(mapping, tp.FindFilename("filename", ""),
		trilerp, maxAniso, wrapMode, scale, gamma)
}

func CreateImageSpectrumTexture(tex2world *Transform, tp *TextureParams) *ImageTextureSpectrum {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1.0, 1.0, 0.0, 0.0)
	}

	// Initialize _ImageTexture_ parameters
	maxAniso := tp.FindFloat("maxanisotropy", 8.0)
	trilerp := tp.FindBool("trilinear", false)
	wrap := tp.FindString("wrap", "repeat")
	wrapMode := TEXTURE_REPEAT
	if strings.Compare(wrap, "black") == 0 {
		wrapMode = TEXTURE_BLACK
	} else if strings.Compare(wrap, "clamp") == 0 {
		wrapMode = TEXTURE_CLAMP
	}
	scale := tp.FindFloat("scale", 1.0)
	gamma := tp.FindFloat("gamma", 1.0)
	return NewImageTextureSpectrum(mapping, tp.FindFilename("filename", ""),
		trilerp, maxAniso, wrapMode, scale, gamma)
}

func CreateMarbleFloatTexture(tex2world *Transform, tp *TextureParams) *MarbleTextureFloat {
	return nil
}
func CreateMarbleSpectrumTexture(tex2world *Transform, tp *TextureParams) *MarbleTextureSpectrum {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return NewMarbleTextureSpectrum(tp.FindInt("octaves", 8), tp.FindFloat("roughness", 0.5), tp.FindFloat("scale", 1.0),
		tp.FindFloat("variation", 0.2), mapping)
}

func CreateMixFloatTexture(tex2world *Transform, tp *TextureParams) *MixTextureFloat {
	return &MixTextureFloat{tp.GetFloatTexture("tex1", 0.0), tp.GetFloatTexture("tex2", 1.0), tp.GetFloatTexture("amount", 0.5)}
}
func CreateMixSpectrumTexture(tex2world *Transform, tp *TextureParams) *MixTextureSpectrum {
	return &MixTextureSpectrum{tp.GetSpectrumTexture("tex1", *NewSpectrum1(0.0)), tp.GetSpectrumTexture("tex2", *NewSpectrum1(1.0)), tp.GetFloatTexture("amount", 0.50)}
}

func CreateScaleFloatTexture(tex2world *Transform, tp *TextureParams) *ScaleTextureFloat {
	tex1 := tp.GetFloatTexture("tex1", 1.0)
	tex2 := tp.GetFloatTexture("tex2", 1.0)
	Assert(tex1 != nil)
	Assert(tex2 != nil)
	return &ScaleTextureFloat{tex1, tex2}
}
func CreateScaleSpectrumTexture(tex2world *Transform, tp *TextureParams) *ScaleTextureSpectrum {
	return &ScaleTextureSpectrum{tp.GetSpectrumTexture("tex1", *NewSpectrum1(1.0)), tp.GetSpectrumTexture("tex2", *NewSpectrum1(1.0))}
}
func CreateUVFloatTexture(tex2world *Transform, tp *TextureParams) *UVTextureFloat {
	return nil
}
func CreateUVSpectrumTexture(tex2world *Transform, tp *TextureParams) *UVTextureSpectrum {
	// Initialize 2D texture mapping _map_ from _tp_
	var mapping TextureMapping2D
	maptype := tp.FindString("mapping", "uv")
	if strings.Compare(maptype, "uv") == 0 {
		su := tp.FindFloat("uscale", 1.0)
		sv := tp.FindFloat("vscale", 1.0)
		du := tp.FindFloat("udelta", 0.0)
		dv := tp.FindFloat("vdelta", 0.0)
		mapping = NewUVMapping2D(su, sv, du, dv)
	} else if strings.Compare(maptype, "spherical") == 0 {
		mapping = NewSphericalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "cylindrical") == 0 {
		mapping = NewCylindricalMapping2D(InverseTransform(tex2world))
	} else if strings.Compare(maptype, "planar") == 0 {
		mapping = NewPlanarMapping2D(tp.FindVector("v1", Vector{1, 0, 0}),
			tp.FindVector("v2", Vector{0, 1, 0}),
			tp.FindFloat("udelta", 0.0), tp.FindFloat("vdelta", 0.0))
	} else {
		Error("2D texture mapping \"%s\" unknown", maptype)
		mapping = NewUVMapping2D(1.0, 1.0, 0.0, 0.0)
	}
	return &UVTextureSpectrum{mapping}
}

func CreateWindyFloatTexture(tex2world *Transform, tp *TextureParams) *WindyTextureFloat {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return &WindyTextureFloat{mapping}
}
func CreateWindySpectrumTexture(tex2world *Transform, tp *TextureParams) *WindyTextureSpectrum {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return &WindyTextureSpectrum{mapping}
}

func CreateWrinkledFloatTexture(tex2world *Transform, tp *TextureParams) *WrinkledTextureFloat {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return &WrinkledTextureFloat{tp.FindInt("octaves", 8), tp.FindFloat("roughness", 0.5), mapping}
}
func CreateWrinkledSpectrumTexture(tex2world *Transform, tp *TextureParams) *WrinkledTextureSpectrum {
	// Initialize 3D texture mapping _map_ from _tp_
	mapping := NewIdentityMapping3D(tex2world)
	return &WrinkledTextureSpectrum{tp.FindInt("octaves", 8), tp.FindFloat("roughness", 0.5), mapping}
}

/*
	Integrators
*/

func CreateAmbientOcclusionIntegrator(params *ParamSet) *AmbientOcclusionIntegrator {
	nSamples := params.FindIntParam("nsamples", 2048)
	maxDist := params.FindFloatParam("maxdist", INFINITY)
	if options.FastRender {
		nSamples = Maxi(1, nSamples/2)
	} else if options.QuickRender {
		nSamples = Maxi(1, nSamples/4)
	}
	return &AmbientOcclusionIntegrator{nSamples, maxDist}
}

func CreateDiffusePRTIntegratorSurfaceIntegrator(params *ParamSet) *DiffusePRTIntegrator {
	lmax := params.FindIntParam("lmax", 4)
	ns := params.FindIntParam("nsamples", 4096)
	integrator := &DiffusePRTIntegrator{lmax, int(RoundUpPow2(uint32(ns))), nil}
	integrator.c_in = make([]Spectrum, SHTerms(lmax), SHTerms(lmax))
	return integrator
}

func CreateGlossyPRTIntegratorSurfaceIntegrator(params *ParamSet) *GlossyPRTIntegrator {
	lmax := params.FindIntParam("lmax", 4)
	ns := params.FindIntParam("nsamples", 4096)
	Kd := params.FindSpectrumParam("Kd", *NewSpectrum1(0.5))
	Ks := params.FindSpectrumParam("Ks", *NewSpectrum1(0.25))
	roughness := params.FindFloatParam("roughness", 0.1)
	return NewGlossyPRTIntegrator(Kd, Ks, roughness, lmax, ns)
}

func CreateIrradianceCacheIntegrator(params *ParamSet) *IrradianceCacheIntegrator {
	minWeight := params.FindFloatParam("minweight", 0.5)
	minSpacing := params.FindFloatParam("minpixelspacing", 2.5)
	maxSpacing := params.FindFloatParam("maxpixelspacing", 15.0)
	maxAngle := params.FindFloatParam("maxangledifference", 10.0)
	maxSpecularDepth := params.FindIntParam("maxspeculardepth", 5)
	maxIndirectDepth := params.FindIntParam("maxindirectdepth", 3)
	nSamples := params.FindIntParam("nsamples", 4096)
	if options.FastRender {
		nSamples = Maxi(1, nSamples/4)
	} else if options.QuickRender {
		nSamples = Maxi(1, nSamples/16)
	}
	return NewIrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing, maxAngle,
		maxSpecularDepth, maxIndirectDepth, nSamples)
}

func CreatePathSurfaceIntegrator(params *ParamSet) *PathIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	return NewPathIntegrator(maxDepth)
}

func CreatePhotonMapSurfaceIntegrator(params *ParamSet) *PhotonIntegrator {
	nCaustic := params.FindIntParam("causticphotons", 20000)
	nIndirect := params.FindIntParam("indirectphotons", 100000)
	nUsed := params.FindIntParam("nused", 50)
	if options.FastRender {
		nCaustic = nCaustic / 5
		nIndirect = nIndirect / 5
		nUsed = Maxi(1, nUsed / 5)		
	} else if options.QuickRender {
		nCaustic = nCaustic / 10
		nIndirect = nIndirect / 10
		nUsed = Maxi(1, nUsed/10)
	}
	maxSpecularDepth := params.FindIntParam("maxspeculardepth", 5)
	maxPhotonDepth := params.FindIntParam("maxphotondepth", 5)
	finalGather := params.FindBoolParam("finalgather", true)
	gatherSamples := params.FindIntParam("finalgathersamples", 32)
	if options.FastRender {
		gatherSamples = Maxi(1, gatherSamples/2)		
	} else if options.QuickRender {
		gatherSamples = Maxi(1, gatherSamples/4)
	}
	maxDist := params.FindFloatParam("maxdist", 0.1)
	gatherAngle := params.FindFloatParam("gatherangle", 10.0)
	return NewPhotonIntegrator(nCaustic, nIndirect, nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples, gatherAngle)
}

func CreateRadianceProbesSurfaceIntegrator(params *ParamSet) *UseRadianceProbes {
	filename := params.FindFilenameParam("filename", "probes.out")
	return NewUseRadianceProbes(filename)
}

func CreateDirectLightingIntegrator(params *ParamSet) *DirectLightingIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	var strategy LightStrategy
	st := params.FindStringParam("strategy", "all")
	if strings.Compare(st, "one") == 0 {
		strategy = SAMPLE_ONE_UNIFORM
	} else if strings.Compare(st, "all") == 0 {
		strategy = SAMPLE_ALL_UNIFORM
	} else {
		Warning("Strategy \"%s\" for direct lighting unknown.  Using \"all\".", st)
		strategy = SAMPLE_ALL_UNIFORM
	}
	return &DirectLightingIntegrator{strategy, maxDepth, nil, nil, 0}
}

func CreateWhittedSurfaceIntegrator(params *ParamSet) *WhittedIntegrator {
	maxDepth := params.FindIntParam("maxdepth", 5)
	return &WhittedIntegrator{maxDepth}
}
