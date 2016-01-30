package core

import (
	"os"
)

func ReadFloatFile(filename string) (bool, []float64) {
    f, err := os.Open(filename)
    defer f.Close()
    
    if err != nil {
        Error("Unable to open file \"%s\"", filename)
        return false, nil
    }

	values := make([]float64, 0, 16)
	
    int c;
    bool inNumber = false;
    char curNumber[32];
    int curNumberPos = 0;
    int lineNumber = 1;
    while ((c = getc(f)) != EOF) {
        if (c == '\n') ++lineNumber;
        if (inNumber) {
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+') {
                curNumber[curNumberPos++] = c;
            } else {
                curNumber[curNumberPos++] = '\0';
                values->push_back(atof(curNumber));
                Assert(curNumberPos < (int)sizeof(curNumber));
                inNumber = false;
                curNumberPos = 0;
            }
        } else {
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {
                inNumber = true;
                curNumber[curNumberPos++] = c;
            } else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)
                    ;
                ++lineNumber;
            } else if (!isspace(c)) {
                Warning("Unexpected text found at line %d of float file \"%s\"",
                        lineNumber, filename);
            }
        }
    }
    return true, values
}

