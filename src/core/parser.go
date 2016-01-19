package core

// Go commands to generate lexer/parser.
//
//go:generate nex -o pbrtlex.go pbrtlex.nn
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

import (
	"os"
	"strings"
	"path/filepath"
)

type Object interface{}

type ParamSet struct {
	tokens []string
	params []Object
}

func (ps *ParamSet) FindStringParam(name, defval string) string {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "string " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if s, ok := values[0].(string); ok {
					value = s
				}
			}
		}
	}
	return value
} 

func (ps *ParamSet) FindFilenameParam(name, defval string) string {
	filename := ps.FindStringParam(name, "")
	if len(filename) == 0 {
		return defval
	}
	return filepath.Clean(filename)
} 

func (ps *ParamSet) FindTextureParam(name string) string {
	if ps == nil {
		return ""
	}
	var value string
	fullparamname := "texture " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if s, ok := values[0].(string); ok {
					value = s
				}
			}
		}
	}
	return value
} 

func (ps *ParamSet) FindFloatParam(name string, defval float64) float64 {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "float " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if val, ok := values[0].(float64); ok {
					value = val					
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindFloatArrayParam(name string) []float64 {
	if ps == nil {
		return nil
	}
	var value []float64
	fullparamname := "float " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				value = make([]float64, len(values), len(values))
				for i, vs := range values {
					if v, ok := vs.(float64); ok {
						value[i] = v
					}
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindIntParam(name string, defval int) int {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "integer " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if val, ok := values[0].(float64); ok {
					value = int(val)					
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindIntArrayParam(name string) []int {
	if ps == nil {
		return nil
	}
	var value []int
	fullparamname := "integer " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				value = make([]int, len(values), len(values))
				for i, vs := range values {
					if v, ok := vs.(float64); ok {
						value[i] = int(v)
					}
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindBoolParam(name string, defval bool) bool {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "bool " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if b, ok := values[0].(bool); ok {
					value = b
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindPointParam(name string, defval Point) Point {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "point " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                if len(values) == 3 {
					if v, ok := values[i*3+0].(float64); ok {
						value.x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						value.y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						value.z = v
					}
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindPointArrayParam(name string) []Point {
	if ps == nil {
		return nil
	}
	var array []Point
	fullparamname := "point " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                numPoints := len(values)/3
				array = make([]Point, numPoints, numPoints)
				for i, _ := range array {
					if v, ok := values[i*3+0].(float64); ok {
						array[i].x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						array[i].y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						array[i].z = v
					}
				}
			}
		}
	}
	return array	
}

func (ps *ParamSet) FindVectorParam(name string, defval Vector) Vector {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "vector " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                if len(values) == 3 {
					if v, ok := values[i*3+0].(float64); ok {
						value.x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						value.y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						value.z = v
					}
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindVectorArrayParam(name string) []Vector {
	if ps == nil {
		return nil
	}
	var array []Vector
	fullparamname := "vector " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                numVectors := len(values)/3
				array = make([]Vector, numVectors, numVectors)
				for i, _ := range array {
					if v, ok := values[i*3+0].(float64); ok {
						array[i].x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						array[i].y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						array[i].z = v
					}
				}
			}
		}
	}
	return array	
}

func (ps *ParamSet) FindNormalParam(name string, defval Normal) Normal {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "normal " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                if len(values) == 3 {
					if v, ok := values[i*3+0].(float64); ok {
						value.x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						value.y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						value.z = v
					}
				}
			}
		}
	}
	return value	
}

func (ps *ParamSet) FindNormalArrayParam(name string) []Normal {
	if ps == nil {
		return nil
	}
	var array []Normal
	fullparamname := "normal " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                numNormals := len(values)/3
				array = make([]Normal, numNormals, numNormals)
				for i, _ := range array {
					if v, ok := values[i*3+0].(float64); ok {
						array[i].x = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						array[i].y = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						array[i].z = v
					}
				}
			}
		}
	}
	return array	
}

func (ps *ParamSet) FindSpectrumParam(name string, defval Spectrum) Spectrum {
	if ps == nil {
		return defval
	}
	value := defval
	fullparamname := "spectrum " + name
	for i, p := range ps.tokens {
		if strings.Compare(p, fullparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
                if len(values) == 3 {
					if v, ok := values[i*3+0].(float64); ok {
						value.c[0] = v
					}
					if v, ok := values[i*3+1].(float64); ok {
						value.c[1] = v
					}
					if v, ok := values[i*3+2].(float64); ok {
						value.c[2] = v
					}
				}
			}
		}
	}
	return value	
}

type TextureParams struct {
	floatTextures        map[string]TextureFloat
	spectrumTextures     map[string]TextureSpectrum
    geomParams, materialParams *ParamSet	
}

func CreateTextureParams(gp, mp *ParamSet, tf map[string]TextureFloat, ts map[string]TextureSpectrum) *TextureParams {
    return &TextureParams{tf, ts, gp, mp}
}

func (tp *TextureParams) FindFloat(name string, defval float64) float64 {
	return tp.geomParams.FindFloatParam(name, tp.materialParams.FindFloatParam(name, defval))	
}
func (tp *TextureParams) FindString(name string, defval string) string {
	return tp.geomParams.FindStringParam(name, tp.materialParams.FindStringParam(name, defval))	
}
func (tp *TextureParams) FindFilename(name string, defval string) string {
	return tp.geomParams.FindFilenameParam(name, tp.materialParams.FindFilenameParam(name, defval))	
}
func (tp *TextureParams) FindInt(name string, defval int) int {
	return tp.geomParams.FindIntParam(name, tp.materialParams.FindIntParam(name, defval))	
}
func (tp *TextureParams) FindBool(name string, defval bool) bool {
	return tp.geomParams.FindBoolParam(name, tp.materialParams.FindBoolParam(name, defval))	
}
func (tp *TextureParams) FindPoint(name string, defval Point) Point {
	return tp.geomParams.FindPointParam(name, tp.materialParams.FindPointParam(name, defval))	
}
func (tp *TextureParams) FindVector(name string, defval Vector) Vector {
	return tp.geomParams.FindVectorParam(name, tp.materialParams.FindVectorParam(name, defval))	
}
func (tp *TextureParams) FindNormal(name string, defval Normal) Normal {
	return tp.geomParams.FindNormalParam(name, tp.materialParams.FindNormalParam(name, defval))	
}
func (tp *TextureParams) FindSpectrum(name string, defval Spectrum) Spectrum {
	return tp.geomParams.FindSpectrumParam(name, tp.materialParams.FindSpectrumParam(name, defval))	
}

func extractFloatParam(parm []Object) (float64, bool) {
	v, ok := parm[0].(float64)
	return v, ok
}
func extractColorParam(parm []Object) (*Spectrum, bool) {
	if len(parm) != 3 {
		return nil, false
	}
	if r, ok := parm[0].(float64); ok {
		if g, ok := parm[1].(float64); ok {
			if b, ok := parm[2].(float64); ok {
				return CreateSpectrumRGB(r, g, b), true
			}
		}
	}
	return nil, false
}
func extractPointParam(parm []Object) (*Vector, bool) {
	if len(parm) != 3 {
		return nil, false
	}
	if x, ok := parm[0].(float64); ok {
		if y, ok := parm[1].(float64); ok {
			if z, ok := parm[2].(float64); ok {
				return &Vector{x, y, z}, true
			}
		}
	}
	return nil, false
}

func ParseFile(filename string) bool {
	fi, err := os.Open(filename)
	defer fi.Close()

	if err != nil {
		return false
	}

	yylex := NewLexer(fi)
	yyParse(yylex)
	return true
}