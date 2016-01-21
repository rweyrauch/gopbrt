package core

// Go commands to generate lexer/parser.
//
//go:generate nex -o pbrtlex.go pbrtlex.nn
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

import (
	"os"
	"path/filepath"
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
				for ii, vs := range values {
					if v, ok := vs.(float64); ok {
						value[ii] = v
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
				for ii, vs := range values {
					if v, ok := vs.(float64); ok {
						value[ii] = int(v)
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
					if v, ok := values[0].(float64); ok {
						value.x = v
					}
					if v, ok := values[1].(float64); ok {
						value.y = v
					}
					if v, ok := values[2].(float64); ok {
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
				numPoints := len(values) / 3
				array = make([]Point, numPoints, numPoints)
				for ii, _ := range array {
					if v, ok := values[ii*3+0].(float64); ok {
						array[ii].x = v
					}
					if v, ok := values[ii*3+1].(float64); ok {
						array[ii].y = v
					}
					if v, ok := values[ii*3+2].(float64); ok {
						array[ii].z = v
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
					if v, ok := values[0].(float64); ok {
						value.x = v
					}
					if v, ok := values[1].(float64); ok {
						value.y = v
					}
					if v, ok := values[2].(float64); ok {
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
				numVectors := len(values) / 3
				array = make([]Vector, numVectors, numVectors)
				for ii, _ := range array {
					if v, ok := values[ii*3+0].(float64); ok {
						array[ii].x = v
					}
					if v, ok := values[ii*3+1].(float64); ok {
						array[ii].y = v
					}
					if v, ok := values[ii*3+2].(float64); ok {
						array[ii].z = v
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
					if v, ok := values[0].(float64); ok {
						value.x = v
					}
					if v, ok := values[1].(float64); ok {
						value.y = v
					}
					if v, ok := values[2].(float64); ok {
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
				numNormals := len(values) / 3
				array = make([]Normal, numNormals, numNormals)
				for ii, _ := range array {
					if v, ok := values[ii*3+0].(float64); ok {
						array[ii].x = v
					}
					if v, ok := values[ii*3+1].(float64); ok {
						array[ii].y = v
					}
					if v, ok := values[ii*3+2].(float64); ok {
						array[ii].z = v
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
						value.c[0] = float32(v)
					}
					if v, ok := values[i*3+1].(float64); ok {
						value.c[1] = float32(v)
					}
					if v, ok := values[i*3+2].(float64); ok {
						value.c[2] = float32(v)
					}
				}
			}
		}
	}
	return value
}

type TextureParams struct {
	floatTextures              map[string]TextureFloat
	spectrumTextures           map[string]TextureSpectrum
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

func (tp *TextureParams) GetFloatTexture(name string, defval float64) TextureFloat {
	texname := tp.geomParams.FindTextureParam(name)
	if len(texname) == 0 {
		texname = tp.materialParams.FindTextureParam(name)
	}
	if len(texname) == 0 {
		if tp.floatTextures[texname] != nil {
			return tp.floatTextures[name]
		} else {
			Error("Couldn't find float texture named \"%s\" for parameter \"%s\"", texname, name)
		}
	}
	val := tp.geomParams.FindFloatParam(name, tp.materialParams.FindFloatParam(name, defval))
	return NewConstantTextureFloat(float32(val))
}

func (tp *TextureParams) GetSpectrumTexture(name string, defval Spectrum) TextureSpectrum {
	texname := tp.geomParams.FindTextureParam(name)
	if len(texname) == 0 {
		texname = tp.materialParams.FindTextureParam(name)
	}
	if len(texname) == 0 {
		if tp.spectrumTextures[texname] != nil {
			return tp.spectrumTextures[texname]
		} else {
			Error("Couldn't find spectrum texture named \"%s\" for parameter \"%s\"", texname, name)
		}
	}
	val := tp.geomParams.FindSpectrumParam(name, tp.materialParams.FindSpectrumParam(name, defval))
	return NewConstantTextureSpectrum(val)
}

func (tp *TextureParams) GetFloatTextureOrNil(name string) TextureFloat {
	texname := tp.geomParams.FindTextureParam(name)
	if len(texname) == 0 {
		texname = tp.materialParams.FindTextureParam(name)
	}
	if len(texname) == 0 {
		return nil
	}
	if tp.floatTextures[texname] != nil {
		return tp.floatTextures[texname]
	} else {
		Error("Couldn't find float texture named \"%s\" for parameter \"%s\"", texname, name)
		return nil
	}
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
