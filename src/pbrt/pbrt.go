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
package main

import (
	"flag"
	"fmt"
	"github.com/rweyrauch/gopbrt/src/core"
	"os"
	"runtime/pprof"
)

func main() {
	var options core.Options
	var filenames []string
	var profile bool
	var profileOutputFile string

	flag.IntVar(&options.NumCores, "ncores", 1, "Number of cores to use.")
	flag.StringVar(&options.ImageFile, "outfile", "", "Output image file.")
	flag.BoolVar(&options.QuickRender, "quick", false, "Quick render mode.")
	flag.BoolVar(&options.Quiet, "quiet", false, "Quiet mode.")
	flag.BoolVar(&options.Verbose, "verbose", false, "Verbose mode.")
	flag.BoolVar(&options.Debug, "debug", false, "Debug mode.")
	flag.BoolVar(&profile, "profile", false, "Enable go profiling.")
	flag.StringVar(&profileOutputFile, "profout", "pbrt.prof", "Profile samples output file.")
	flag.Parse()

	if profile {
		f, err := os.Create(profileOutputFile)
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	for _, arg := range flag.Args() {
		filenames = append(filenames, arg)
	}

	// Print welcome banner
	if !options.Quiet {
		fmt.Printf("gopbrt version %s\n", core.GOPBRT_VERSION)
		fmt.Printf("Copyright (c)2016 Rick Weyrauch.\n\n")
		fmt.Printf("gopbrt based on pbrt 2.0.0 (see http://pbrt.org)\n")
		fmt.Printf("Copyright (c)1998-2014 Matt Pharr and Greg Humphreys.\n")
		fmt.Printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n")
		fmt.Printf("See the file LICENSE.txt for the conditions of the license.\n")
	}

	if options.Debug {
		options.Quiet = false
		options.Verbose = true
	}

	core.PbrtInit(&options)
	// Process scene description
	if len(filenames) == 0 {
		// Parse scene from standard input
		core.ParseFile("-")
	} else {
		// Parse scene from input files
		for _, f := range filenames {
			if !core.ParseFile(f) {
				core.Error("Couldn't open scene file \"%s\"", f)
			}
		}
	}
	core.PbrtCleanup()
}
