package main

import (
    "github.com/rweyrauch/gopbrt/src/core"
    "flag"
    "fmt"
)


func main() {
	var options pbrt.Options
	var filenames []string
	
	flag.IntVar(&options.NumCores, "ncores", 1, "Number of cores to use.")
	flag.StringVar(&options.ImageFile, "outfile", "output.png", "Output image file.")
	flag.BoolVar(&options.QuickRender, "quick", false, "Quick render mode.")
	flag.BoolVar(&options.Quiet, "quiet", false, "Quiet mode.")
	flag.BoolVar(&options.Verbose, "verbose", false, "Verbose mode.")
	
	flag.Parse()
	
	for _, arg := range flag.Args() {
		filenames = append(filenames, arg)
	}
	
    // Print welcome banner
    if !options.Quiet {
        fmt.Printf("gopbrt version %s [Detected %d core(s)]\n",
               pbrt.PBRT_VERSION, pbrt.NumSystemCores())
        fmt.Printf("gopbrt based on pbrt 2.0.0 by Matt Pharr ang Grey Humphreys.\n")
        fmt.Printf("Copyright (c)1998-2014 Matt Pharr and Greg Humphreys.\n")
        fmt.Printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n")
        fmt.Printf("See the file LICENSE.txt for the conditions of the license.\n")
    }
	
    pbrt.PbrtInit(&options)
    // Process scene description
    //PBRT_STARTED_PARSING();
    if len(filenames) == 0 {
        // Parse scene from standard input
        pbrt.ParseFile("-")
    } else {
        // Parse scene from input files
        for _, f := range filenames {
            if !pbrt.ParseFile(f) {
                fmt.Printf("Couldn't open scene file \"%s\"", f)
            }
        }
    }
    pbrt.PbrtCleanup()
}
