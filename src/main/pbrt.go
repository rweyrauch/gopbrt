package main

import (
    "github.com/rweyrauch/pbrt"
    "flag"
)


func main() {
	var options pbrt.Options
	var filenames []string
	
	flag.IntVar(&options.nCores, "ncores", 1, "Number of cores to use.")
	flag.StringVar(&options.imageFile, "outfile", "output.png", "Output image file.")
	flag.BoolVar(&options.quickRender, "quick", false, "Quick render mode.")
	flag.BoolVar(&options.quiet, "quiet", false, "Quiet mode.")
	flag.BoolVar(&options.verbose, "verbose", false, "Verbose mode.")
	
	flag.Parse()
	
	for _, arg := range flag.Args() {
		filenames = append(filenames, arg)
	}
	
    // Print welcome banner
    if !options.quiet {
        fmt.Printf("gopbrt version %s [Detected %d core(s)]\n",
               pbrt.PBRT_VERSION, prbt.NumSystemCores())
        fmt.Printf("gopbrt based on pbrt 2.0.0 by Matt Pharr ang Grey Humphreys.\n")
        fmt.Printf("Copyright (c)1998-2014 Matt Pharr and Greg Humphreys.\n")
        fmt.Printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n")
        fmt.Printf("See the file LICENSE.txt for the conditions of the license.\n")
    }
	
   pbrt.Init(options)
    // Process scene description
    PBRT_STARTED_PARSING();
    if len(filenames) == 0 {
        // Parse scene from standard input
        ParseFile("-")
    } else {
        // Parse scene from input files
        for _, f := range filenames {
            if !ParseFile(f) {
                fmt.Printf("Couldn't open scene file \"%s\"", f)
            }
        }
    }
    pbrt.Cleanup()
}
