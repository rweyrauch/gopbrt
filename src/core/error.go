package core

import (
	"fmt"
)

func Debug(format string, args ...interface{}) {
    if !options.Debug { 
    	return
    }
    fmt.Printf(format, args...)
    fmt.Println("")	
}

func Info(format string, args ...interface{}) {
    if !options.Verbose || options.Quiet { 
    	return
    }
    fmt.Printf(format, args...)
    fmt.Println("")	
}

func Warning(format string, args ...interface{}) {
    if options.Quiet {
    	return
  	}
    fmt.Printf(format, args...)
    fmt.Println("")	
}

func Error(format string, args ...interface{}) {
    fmt.Printf(format, args...)
    fmt.Println("")	
}

func Severe(format string, args ...interface{}) {
    panic(fmt.Errorf(format, args...))
}

func Assert(assertion bool) {
	if !assertion { panic("Assertion failed") }
}
func AssertMsg(assertion bool, msg string) {
	if !assertion { panic(fmt.Errorf("Assertion failed: %s", msg)) }
}
func Unimplemented() {
	panic(fmt.Errorf("Unimplemented."))
}