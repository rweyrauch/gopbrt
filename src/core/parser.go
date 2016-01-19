package core

// Go commands to generate lexer/parser.
//
//go:generate nex -o pbrtlex.go pbrtlex.nn
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

import (
	"os"
	"strings"
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

type TextureParams struct {
	floatTextures        map[string]TextureFloat
	spectrumTextures     map[string]TextureSpectrum
    geomParams, materialParams *ParamSet	
}

func CreateTextureParams(gp, mp *ParamSet, tf map[string]TextureFloat, ts map[string]TextureSpectrum) *TextureParams {
    return &TextureParams{tf, ts, gp, mp}
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
