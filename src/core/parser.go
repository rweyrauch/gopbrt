package core

// Go commands to generate lexer/parser.
//
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

import (
	"os"
	"io/ioutil"
	"path/filepath"
	"strings"
	"strconv"
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
				if b, ok := values[0].(string); ok {
					if strings.Compare(b, "true") == 0 {
						value = true
					} else if strings.Compare(b, "false") == 0 {
						value = false
					}
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

	spectrumparamname := "spectrum " + name
	colorparamname := "color " + name
	rgbparamname := "rgb " + name
	bbodyparamname := "blockbody " + name
	xyzparamname := "xyz " + name

	for i, p := range ps.tokens {
		if strings.Compare(p, spectrumparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if len(values) > 0 {
					if _, ok := values[0].(string); ok {
						// TODO: load and read spectrum data file.
						Unimplemented()
					} else {
						// TODO: read pairs of sampld spectrum data
						Unimplemented()
					}
				}
			}
		} else if strings.Compare(p, colorparamname) == 0 ||
			strings.Compare(p, rgbparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if len(values) == 3 {
					if v, ok := values[0].(float64); ok {
						value.c[0] = v
					}
					if v, ok := values[1].(float64); ok {
						value.c[1] = v
					}
					if v, ok := values[2].(float64); ok {
						value.c[2] = v
					}
				} else {
					Error("RGB values given with parameter \"%s\" expected 3 value, got %d.", name, len(values))
				}
			}
		} else if strings.Compare(p, bbodyparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if len(values)%2 != 0 {
					Warning("Excess value given with blackbody parameter \"%s\".  Ignoring extra one.", name)
				}
				lambdas := make([]float64, len(values)/2, len(values)/2)
				vals := make([]float64, len(values)/2, len(values)/2)
				for vi := 0; vi < len(values); vi = vi + 2 {
					if v, ok := values[vi].(float64); ok {
						lambdas[vi/2] = v
					}
					if v, ok := values[vi+1].(float64); ok {
						vals[vi/2] = v
					}
				}
				value = *SpectrumFromSampled(lambdas, vals)
			}
		} else if strings.Compare(p, xyzparamname) == 0 {
			if values, ok := ps.params[i].([]Object); ok {
				if len(values) == 3 {
					xyz := [3]float64{0.0, 0.0, 0.0}
					if v, ok := values[0].(float64); ok {
						xyz[0] = v
					}
					if v, ok := values[1].(float64); ok {
						xyz[1] = v
					}
					if v, ok := values[2].(float64); ok {
						xyz[2] = v
					}
					value = *SpectrumFromXYZ(xyz)
				} else {
					Error("XYZ values given with parameter \"%s\" expected 3 value, got %d.", name, len(values))
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
	if len(texname) != 0 {
		if tp.floatTextures[texname] != nil {
			return tp.floatTextures[name]
		} else {
			Error("Couldn't find float texture named \"%s\" for parameter \"%s\"", texname, name)
		}
	}
	val := tp.geomParams.FindFloatParam(name, tp.materialParams.FindFloatParam(name, defval))
	return NewConstantTextureFloat(val)
}

func (tp *TextureParams) GetSpectrumTexture(name string, defval Spectrum) TextureSpectrum {
	texname := tp.geomParams.FindTextureParam(name)
	if len(texname) == 0 {
		texname = tp.materialParams.FindTextureParam(name)
	}
	if len(texname) != 0 {
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

type Lexer struct {
	scanner *Scanner
}

func (yylex Lexer) Error(e string) {
	panic(e)
}
func (yylex Lexer) Lex(lval *yySymType) int {
scanagain:	
	_, token, val := yylex.scanner.Scan()
	
	if token == EOF { return 0 }
	switch token {
	case NUMBER:
		value, _ := strconv.ParseFloat(val, 64) 
		lval.value = value 	
	case IDENTIFIER:
		lval.id = val
	case STRING:
		lval.id = val
	case COMMENT:
		goto scanagain				
	}
	
	return token
}

func ParseFile(filename string) bool {
	fi, err := os.Open(filename)
	defer fi.Close()

	if err != nil {
		return false
	}

	var yylex Lexer
	src, err := ioutil.ReadAll(fi)
	yylex.scanner = new(Scanner)
	yylex.scanner.Init(filename, src, nil)
	yyParse(yylex)
	return true
}
