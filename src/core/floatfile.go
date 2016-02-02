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
	"bufio"
	"os"
	"strconv"
)

func isSpace(ch rune) bool {
	if ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' {
		return true
	} else {
		return false
	}
}

func ReadFloatFile(filename string) (bool, []float64) {
	fi, err := os.Open(filename)
	if err != nil {
		Error("Unable to open file \"%s\"", filename)
		return false, nil
	}

	values := make([]float64, 0, 16)
	rdr := bufio.NewReader(fi)

	inNumber := false
	var curNumber string
	lineNumber := 1
	c, err := rdr.ReadByte()
	for err == nil {
		if c == '\n' {
			lineNumber++
		}
		if inNumber {
			if isDigit(rune(c)) || c == '.' || c == 'e' || c == '-' || c == '+' {
				curNumber = curNumber + string(c)
			} else {
				v, err := strconv.ParseFloat(curNumber, 64)
				if err == nil {
					values = append(values, v)
					inNumber = false
					curNumber = ""
				}
			}

		} else {
			if isDigit(rune(c)) || c == '.' || c == '-' || c == '+' {
				inNumber = true
				curNumber = curNumber + string(c)
			} else if c == '#' {
				_, err = rdr.ReadBytes('\n')
				lineNumber++
			} else if !isSpace(rune(c)) {
				Warning("Unexpected text found at line %d of float file \"%s\"", lineNumber, filename)
			}
		}
		c, err = rdr.ReadByte()
	}
	return true, values
}
