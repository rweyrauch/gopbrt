package api

// Go commands to generate lexer/parser.
//
//go:generate nex -o pbrtlex.go pbrtlex.nn
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

import (
	"os"
    "github.com/rweyrauch/gopbrt/src/core"
)

func extractFloatParam(parm []pbrt.Object) (float64, bool) {
	v, ok := parm[0].(float64)
	return v, ok
}
func extractColorParam(parm []pbrt.Object) (*pbrt.Spectrum, bool) {
	if len(parm) != 3 {
		return nil, false
	}
	if r, ok := parm[0].(float64); ok {
		if g, ok := parm[1].(float64); ok {
			if b, ok := parm[2].(float64); ok {
				return pbrt.CreateSpectrumRGB(r, g, b), true
			}
		}
	}
	return nil, false
}
func extractPointParam(parm []pbrt.Object) (*pbrt.Vector, bool) {
	if len(parm) != 3 {
		return nil, false
	}
	if x, ok := parm[0].(float64); ok {
		if y, ok := parm[1].(float64); ok {
			if z, ok := parm[2].(float64); ok {
				return &pbrt.Vector{x, y, z}, true
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
