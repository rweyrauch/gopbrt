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
	"fmt"
	"strings"
	"unicode"
	"unicode/utf8"
)

const EOF int = -1
const ILLEGAL int = -2
const COMMENT int = -3

var keywords map[string]int = map[string]int{
	"Accelerator":        ACCELERATOR,
	"ActiveTransform":    ACTIVETRANSFORM,
	"All":                ALL,
	"AreaLightSource":    AREALIGHTSOURCE,
	"AttributeBegin":     ATTRIBUTEBEGIN,
	"AttributeEnd":       ATTRIBUTEEND,
	"Camera":             CAMERA,
	"ConcatTransform":    CONCATTRANSFORM,
	"CoordinateSystem":   COORDINATESYSTEM,
	"CoordSysTransform":  COORDSYSTRANSFORM,
	"EndTime":            ENDTIME,
	"Film":               FILM,
	"Identity":           IDENTITY,
	"Include":            INCLUDE,
	"LightSource":        LIGHTSOURCE,
	"LookAt":             LOOKAT,
	"MakeNamedMaterial":  MAKENAMEDMATERIAL,
	"Material":           MATERIAL,
	"NamedMaterial":      NAMEDMATERIAL,
	"ObjectBegin":        OBJECTBEGIN,
	"ObjectEnd":          OBJECTEND,
	"ObjectInstance":     OBJECTINSTANCE,
	"PixelFilter":        PIXELFILTER,
	"Renderer":           RENDERER,
	"ReverseOrientation": REVERSEORIENTATION,
	"Rotate":             ROTATE,
	"Sampler":            SAMPLER,
	"Scale":              SCALE,
	"Shape":              SHAPE,
	"StartTime":          STARTTIME,
	"SurfaceIntegrator":  SURFACEINTEGRATOR,
	"Texture":            TEXTURE,
	"TransformBegin":     TRANSFORMBEGIN,
	"TransformEnd":       TRANSFORMEND,
	"TransformTimes":     TRANSFORMTIMES,
	"Transform":          TRANSFORM,
	"Translate":          TRANSLATE,
	"Volume":             VOLUME,
	"VolumeIntegrator":   VOLUMEINTEGRATOR,
	"WorldBegin":         WORLDBEGIN,
	"WorldEnd":           WORLDEND,
}

type Position struct {
	Filename string // filename, if any
	Offset   int    // offset, starting at 0
	Line     int    // line number, starting at 1
	Column   int    // column number, starting at 1 (byte count)
}

func (pos *Position) IsValid() bool { return pos.Line > 0 }
func (pos Position) String() string {
	s := pos.Filename
	if pos.IsValid() {
		if s != "" {
			s += ":"
		}
		s += fmt.Sprintf("%d:%d", pos.Line, pos.Column)
	}
	if s == "" {
		s = "-"
	}
	return s
}

type Pos int

const NoPos Pos = 0

type Scanner struct {
	filename string

	src []byte       // source
	err ErrorHandler // error reporting; or nil

	ch         rune
	offset     int // character offset
	rdOffset   int // reading offset (position after current character)
	lineOffset int // current line offset

	ErrorCount int
}

func (s *Scanner) next() {
	if s.rdOffset < len(s.src) {
		s.offset = s.rdOffset
		if s.ch == '\n' {
			s.lineOffset = s.offset
			//s.file.AddLine(s.offset)
		}
		r, w := rune(s.src[s.rdOffset]), 1
		switch {
		case r == 0:
			s.error(s.offset, "illegal character NUL")
		case r >= 0x80:
			// not ASCII
			r, w = utf8.DecodeRune(s.src[s.rdOffset:])
			if r == utf8.RuneError && w == 1 {
				s.error(s.offset, "illegal UTF-8 encoding")
			}
		}
		s.rdOffset += w
		s.ch = r
	} else {
		s.offset = len(s.src)
		if s.ch == '\n' {
			s.lineOffset = s.offset
			//s.file.AddLine(s.offset)
		}
		s.ch = -1 // eof
	}
}

type ErrorHandler func(pos Position, msg string)

func (s *Scanner) error(offs int, msg string) {
	if s.err != nil {
		//s.err(s.file.Position(s.file.Pos(offs)), msg)
	}
	s.ErrorCount++
}

func (s *Scanner) Init(filename string, src []byte, err ErrorHandler) {
	// Explicitly initialize all fields since a scanner may be reused.
	s.filename = filename
	s.src = src
	s.err = err
	s.ch = ' '
	s.offset = 0
	s.rdOffset = 0
	s.lineOffset = 0
	s.ErrorCount = 0
	s.next()
}

func (s *Scanner) skipWhitespace() {
	for s.ch == ' ' || s.ch == '\t' || s.ch == '\n' || s.ch == '\r' {
		s.next()
	}
}
func isLetter(ch rune) bool {
	return 'a' <= ch && ch <= 'z' || 'A' <= ch && ch <= 'Z' || ch == '_' || ch >= 0x80 && unicode.IsLetter(ch)
}
func isDigit(ch rune) bool {
	return '0' <= ch && ch <= '9' || ch >= 0x80 && unicode.IsDigit(ch)
}

func (s *Scanner) scanIdentifier() string {
	offs := s.offset
	for isLetter(s.ch) || isDigit(s.ch) {
		s.next()
	}
	return string(s.src[offs:s.offset])
}

func digitVal(ch rune) int {
	switch {
	case '0' <= ch && ch <= '9':
		return int(ch - '0')
	case 'a' <= ch && ch <= 'f':
		return int(ch - 'a' + 10)
	case 'A' <= ch && ch <= 'F':
		return int(ch - 'A' + 10)
	}
	return 16 // larger than any legal digit val
}
func (s *Scanner) scanMantissa(base int) {
	for digitVal(s.ch) < base {
		s.next()
	}
}

func (s *Scanner) scanNumber(seenDecimalPoint bool) (int, string) {
	// digitVal(s.ch) < 10
	offs := s.offset
	tok := NUMBER
	if seenDecimalPoint {
		offs--
		s.scanMantissa(10)
		goto fraction
	}
	if s.ch == '-' || s.ch == '+' {
		s.next()
	}
	if s.ch == '0' {
		// int or float
		offs := s.offset
		s.next()
		if s.ch == 'x' || s.ch == 'X' {
			// hexadecimal int
			s.next()
			s.scanMantissa(16)
			if s.offset-offs <= 2 {
				// only scanned "0x" or "0X"
				s.error(offs, "illegal hexadecimal number")
			}
		} else {
			// octal int or float
			seenDecimalDigit := false
			s.scanMantissa(8)
			if s.ch == '8' || s.ch == '9' {
				// illegal octal int or float
				seenDecimalDigit = true
				s.scanMantissa(10)
			}
			if s.ch == '.' || s.ch == 'e' || s.ch == 'E' || s.ch == 'i' {
				goto fraction
			}
			// octal int
			if seenDecimalDigit {
				s.error(offs, "illegal octal number")
			}
		}
		goto exit
	}
	// decimal int or float
	s.scanMantissa(10)
fraction:
	if s.ch == '.' {
		s.next()
		s.scanMantissa(10)
	}
	if s.ch == 'e' || s.ch == 'E' {
		s.next()
		if s.ch == '-' || s.ch == '+' {
			s.next()
		}
		s.scanMantissa(10)
	}
exit:
	return tok, string(s.src[offs:s.offset])
}

func (s *Scanner) scanString() string {
	// '"' opening already consumed
	offs := s.offset - 1
	for {
		ch := s.ch
		if ch == '\n' || ch < 0 {
			s.error(offs, "string literal not terminated")
			break
		}
		s.next()
		if ch == '"' {
			break
		}
	}
	return string(s.src[offs:s.offset])
}

func (s *Scanner) scanComment() {
	s.next()
	for s.ch != '\n' && s.ch >= 0 {
		s.next()
	}
}

func LookupToken(name string) int {
	token := keywords[name]
	if token != 0 {
		return token
	}
	return ILLEGAL
}

func (s *Scanner) Scan() (pos Pos, tok int, lit string) {
	s.skipWhitespace()
	// current token start
	pos = 0 //s.file.Pos(s.offset)
	// determine token value
	switch ch := s.ch; {
	case isLetter(ch):
		lit = s.scanIdentifier()
		if len(lit) > 1 {
			// keywords are longer than one letter - avoid lookup otherwise
			tok = LookupToken(lit)
		} else {
			tok = IDENTIFIER
		}
	case '0' <= ch && ch <= '9' || ch == '-':
		tok, lit = s.scanNumber(false)
	default:
		s.next() // always make progress
		switch ch {
		case -1:
			tok = EOF
		case '"':
			tok = STRING
			lit = strings.TrimRight(strings.TrimLeft(s.scanString(), "\""), "\"")
		case '.':
			if '0' <= s.ch && s.ch <= '9' {
				tok, lit = s.scanNumber(true)
			}
		case '[':
			tok = LBRACK
			lit = "["
		case ']':
			tok = RBRACK
			lit = "]"
		case '#':
			s.scanComment()
			tok = COMMENT
			lit = ""
		default:
			// next reports unexpected BOMs - don't repeat
			//s.error(s.file.Offset(pos), fmt.Sprintf("illegal character %#U", ch))
			tok = ILLEGAL
			lit = string(ch)
		}
	}
	return pos, tok, lit
}
