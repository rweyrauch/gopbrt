package pbrt

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
