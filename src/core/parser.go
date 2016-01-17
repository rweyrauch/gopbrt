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

func (ps *ParamSet) FindFloatArrayParam(name string, defval []float64) []float64 {
	if ps == nil {
		return defval
	}
	value := defval
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

func (ps *ParamSet) FindIntArrayParam(name string, defval []int) []int {
	if ps == nil {
		return defval
	}
	value := defval
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
