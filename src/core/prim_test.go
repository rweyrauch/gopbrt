package core

import(
	"testing"
)

var nextVal int = 0

func appendToSlice(prims *[]int) {
	*prims = append(*prims, nextVal)
	nextVal++
}

func TestRefine(t *testing.T) {

	prims := make([]int, 0, 8)
	
	appendToSlice(&prims)
	if len(prims) != 1 {
		t.Fail()
	}	
	
	appendToSlice(&prims)
	if len(prims) != 2 {
		t.Fail()
	}
}
