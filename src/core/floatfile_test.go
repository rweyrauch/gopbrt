package core

import (
	"testing"
	"fmt"
	"path/filepath"
)

func TestFloatFileRead(t *testing.T) {
	var options *Options
	PbrtInit(options)
	
	ok, values := ReadFloatFile(filepath.Join("testdata", "test.spd"))
	if ok {
		for i, v := range values {
			fmt.Printf("V[%d]: %f\n", i, v)
		}
	} else {
		t.Fail()
	}
	
	PbrtCleanup()
}