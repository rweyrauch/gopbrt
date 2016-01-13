%{
package pbrt

import (
	"fmt"	
)

/*
If this file is not pbrtparser.y, it was generated from pbrtparser.y and
should not be edited directly.
*/

func init() {
	yyDebug = 0
	yyErrorVerbose = true
}

type pbrtParameter struct {
	name string
	values Object
}

func splitParamList(paramlist []pbrtParameter) (*ParamSet, bool) {
	// an empty list is valid
	if paramlist == nil {
		return nil, true
	}
	
	pset := new(ParamSet)
	pset.tokens = make([]string, 0, len(paramlist))
	pset.params = make([]Object, 0, len(paramlist))
	
	for _, p := range paramlist {
		pset.tokens = append(pset.tokens, p.name)
		pset.params = append(pset.params, p.values)		
	}
	ok := false
	if len(pset.tokens) == len(pset.params) {
		ok = true
	} else {
		fmt.Errorf("Token and param lists do not have the same length. %d vs. %d\n.", len(pset.tokens), len(pset.params))
	}
	return pset, ok
}

%}

%union {
	value float64
	id string
	params []pbrtParameter
	numbers []float64
	tokens []string
	param pbrtParameter
	objects []Object
}

%token <id> STRING IDENTIFIER
%token <value> NUMBER 
%token LBRACK RBRACK

%token ACCELERATOR ACTIVETRANSFORM ALL AREALIGHTSOURCE ATTRIBUTEBEGIN
%token ATTRIBUTEEND CAMERA CONCATTRANSFORM COORDINATESYSTEM COORDSYSTRANSFORM
%token ENDTIME FILM IDENTITY INCLUDE LIGHTSOURCE LOOKAT MAKENAMEDMATERIAL
%token MATERIAL NAMEDMATERIAL OBJECTBEGIN OBJECTEND OBJECTINSTANCE PIXELFILTER
%token RENDERER REVERSEORIENTATION ROTATE SAMPLER SCALE SHAPE STARTTIME
%token SURFACEINTEGRATOR TEXTURE TRANSFORMBEGIN TRANSFORMEND TRANSFORMTIMES
%token TRANSFORM TRANSLATE VOLUME VOLUMEINTEGRATOR WORLDBEGIN WORLDEND

%type <numbers> number_list number_array
%type <tokens> string_list string_array
%type <id> single_element_string_array
%type <value> single_element_number_array

%type <params> param_list
%type <param> param_list_entry
%type <objects> array

%%

start
	: pbrt_stmt_list
	;

array
	: string_array
	{
		// convert array of stings to array of interfaces
		objects := make([]Object, len($1), len($1))
		for i,v := range $1 {
			objects[i] = v
		}
		$$ = objects
	}
	| number_array
	{
		// convert array of floats to array of interfaces
		objects := make([]Object, len($1), len($1))
		for i,v := range $1 {
			objects[i] = v
		}
		$$ = objects
	}
	;

string_array
	: LBRACK string_list RBRACK
	{
		$$ = $2
	}
	| single_element_string_array
	{
		$$ = append($$, $1)
	}
	;

single_element_string_array
	: STRING
	{
		$$ = $1
	} 
	;
	
string_list
	: string_list STRING
	{
		$$ = append($1, $2)
	}
	| STRING
	{
		$$ = make([]string, 1, 16)
		$$[0] = $1
	}
	;

number_array
	: LBRACK number_list RBRACK
	{
		$$ = $2
	}
	| single_element_number_array
	{
		$$ = append($$, $1)
	}
	;

single_element_number_array
	: NUMBER
	{
		$$ = $1
	}
	;
	
number_list
	: number_list NUMBER
	{
		$$ = append($1, $2)
	}
	| NUMBER
	{
		$$ = make([]float64, 1, 16)
		$$[0] = float64($1)
	}
	;
	
param_list
	: param_list_entry param_list
	{
		$$ = append($2, $1)
	}
	|
	{
		// empty list
		$$ = nil
	}
	;

param_list_entry
	: STRING array
	{
		$$ = pbrtParameter{$1, $2}
	}
	;

pbrt_stmt_list
	: pbrt_stmt_list pbrt_stmt
	| pbrt_stmt
	;

pbrt_stmt
	: ACCELERATOR STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_REFLECTANCE)
			PbrtAccelerator($2, params)
		}
	}
	| ACTIVETRANSFORM ALL
	{
		PbrtActiveTransformAll()
	}
	| ACTIVETRANSFORM ENDTIME
	{
		PbrtActiveTransformEndTime()
	}
	| ACTIVETRANSFORM STARTTIME
	{
		PbrtActiveTransformStartTime()
	}
	| AREALIGHTSOURCE STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_ILLUMINANT)
			PbrtAreaLightSource($2, params)
		}
	}	
	| ATTRIBUTEBEGIN
	{
		PbrtAttributeBegin()
	}
	| ATTRIBUTEEND
	{
		PbrtAttributeEnd()
	}
	| CAMERA STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_REFLECTANCE)
			PbrtCamera($2, params)
		}	
	}
	| CONCATTRANSFORM number_array
	{
		values := $2
		if len(values) == 16 {		
			var matrix Matrix4x4
			i := 0
			for r := range matrix.m {
				for c := range matrix.m[r] {
					matrix.m[r][c] = values[i]
					i++
				}
			}		
			PbrtConcatTransform(matrix)
		} else {
			// TODO: error - require 16 values
			fmt.Printf("Array argument to ConcatTransform requires 16 values.\n")
		}
	}
	| COORDINATESYSTEM STRING
	{
		PbrtCoordinateSystem($2)
	}
	| COORDSYSTRANSFORM STRING
	{
		PbrtCoordSysTransform($2)
	}
	| FILM STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_REFLECTANCE)
			PbrtFilm($2, params)
		}		
	}
	| IDENTITY
	{
		PbrtIdentity()
	}
	| INCLUDE STRING
	{
		//include_push($2)
	}
	| LIGHTSOURCE STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_ILLUMINANT)
			PbrtLightSource($2, params)
		} else {
			fmt.Printf("Failed to parse parameter list for light source.\n")
		}
	}
	| LOOKAT NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER
	{
    	PbrtLookAt($2, $3, $4, $5, $6, $7, $8, $9, $10)
	}
	| MAKENAMEDMATERIAL STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtMakeNamedMaterial($2, params)
    	}
	}
	| MATERIAL STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
			//InitParamSet(params, SPECTRUM_REFLECTANCE)
	    	PbrtMaterial($2, params)
	   	}
	}
	| NAMEDMATERIAL STRING
	{
    	PbrtNamedMaterial($2)
	}
	| OBJECTBEGIN STRING
	{
		PbrtObjectBegin($2)
	}
	| OBJECTEND
	{
		PbrtObjectEnd()
	}
	| OBJECTINSTANCE STRING
	{
		PbrtObjectInstance($2)
	}
	| PIXELFILTER STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
    		//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtPixelFilter($2, params)
    	}
	}
	| RENDERER STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
    		//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtRenderer($2, params)
    	}
	}
	| REVERSEORIENTATION
	{
		PbrtReverseOrientation()
	}
	| ROTATE NUMBER NUMBER NUMBER NUMBER
	{
		PbrtRotate($2, $3, $4, $5)
	}	
	| SAMPLER STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
    		//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtSampler($2, params)
    	}
	}
	| SCALE NUMBER NUMBER NUMBER
	{
		PbrtScale($2, $3, $4)
	}	
	| SHAPE STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
    		//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtShape($2, params)
    	}
	}
	| SURFACEINTEGRATOR STRING param_list
	{
 		params, ok := splitParamList($3)
		if ok {
  	  		//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtSurfaceIntegrator($2, params)
    	}
	}
	| TEXTURE STRING STRING STRING param_list
	{
		params, ok := splitParamList($5)
		if ok {
	   	 	//InitParamSet(params, SPECTRUM_REFLECTANCE)
			PbrtTexture($2, $3, $4, params)
		}
	}
	| TRANSFORMBEGIN
	{
		PbrtTransformBegin()
	}
	| TRANSFORMEND
	{
		PbrtTransformEnd()
	}
	| TRANSFORMTIMES NUMBER NUMBER
	{
		PbrtTransformTimes($2, $3)
	}
	| TRANSFORM number_array
	{
		values := $2
		if len(values) == 16 {		
			var matrix Matrix4x4
			i := 0
			for r := range matrix.m {
				for c := range matrix.m[r] {
					matrix.m[r][c] = values[i]
					i++
				}
			}		
			PbrtTransform(matrix)
		} else {
			// TODO: error - require 16 values
			fmt.Printf("Array argument to Transform requires 16 values.\n")
		}
	}
	| TRANSLATE NUMBER NUMBER NUMBER
	{
		PbrtTranslate($2, $3, $4)
	}
	| VOLUME STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
	   	 	//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtVolume($2, params)
		}
	}
	| VOLUMEINTEGRATOR STRING param_list
	{
		params, ok := splitParamList($3)
		if ok {
	   	 	//InitParamSet(params, SPECTRUM_REFLECTANCE)
    		PbrtVolumeIntegrator($2, params)
		}
	}
	| WORLDBEGIN
	{
		PbrtWorldBegin()
	}
	| WORLDEND
	{
		PbrtWorldEnd()
	}
	;

%%
