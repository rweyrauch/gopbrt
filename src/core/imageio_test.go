package core

import (
	"testing"
	"fmt"
	"path/filepath"
)

func TestLoadExrImage(t *testing.T) {
	var options Options	
	options.Debug = true
	PbrtInit(&options)
	
	image, width, height := ReadImage(filepath.Join("testdata", "test.exr"))
	
	fmt.Printf("Image size: %d x %d.  Expect: 32 x 48.\n", width, height)
	
	if width != 32 {
		t.Fail()
	}
	if height != 48 {
		t.Fail()
	}
	if image == nil {
		t.Fail()
	}
	PbrtCleanup()
}
