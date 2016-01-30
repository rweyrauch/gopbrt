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
