package pbrt

// Go commands to generate lexer/parser.
//
//go:generate nex -o pbrtlex.go pbrtlex.nn
//go:generate go tool yacc -o pbrtparser.go pbrtparser.y

type pbrtPointer interface{}

type ParamSet struct {
	tokens []string
	params []pbrtPointer
}

func extractFloatParam(parm []pbrtPointer) (float64, bool) {
        v, ok := parm[0].(float64)
        return float64(v), ok
}
func extractColorParam(parm []pbrtPointer) (*Spectrum, bool) {
        if len(parm) != 3 {
                return nil, false
        }
        if r, ok := parm[0].(float64); ok {
                if g, ok := parm[1].(float64); ok {
                        if b, ok := parm[2].(float64); ok {
                                return CreateSpectrumRGB(float64(r), float64(g), float64(b)), true
                        }
                }
        }
        return nil, false
}
func extractPointParam(parm []pbrtPointer) (*Vector, bool) {
        if len(parm) != 3 {
                return nil, false
        }
        if x, ok := parm[0].(float64); ok {
                if y, ok := parm[1].(float64); ok {
                        if z, ok := parm[2].(float64); ok {
                                return &Vector{float64(x), float64(y), float64(z)}, true
                        }
                }
        }
        return nil, false
}

func ParseFile(filename string) bool {

	return false	
}