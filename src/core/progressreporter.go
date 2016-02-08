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
package core

import (
	"fmt"
	"math"
	"time"
	"github.com/rweyrauch/gopbrt/src/os"
)

func TerminalWidth() int {
	return os.TerminalWidth()
}

type ProgressReporter struct {
	totalWork              int
	workDone, totalPlusses int
	title                  string
	startTime              time.Time
}

func NewProgressReporter(totalwork int, title string, barLength int) *ProgressReporter {
	pr := new(ProgressReporter)
	pr.totalWork = totalwork

	if barLength <= 0 {
		barLength = TerminalWidth() - 28
	}
	pr.title = title
	pr.totalPlusses = Maxi(2, barLength-len(title))
	pr.workDone = 0
	pr.startTime = time.Now()

	return pr
}

const (
	COMPLETED = "+"
	REMAINING = " "
)

func (pr *ProgressReporter) Update(num int) {
	if num == 0 || options.Quiet {
		return
	}
	pr.workDone += num
	percentDone := float64(pr.workDone) / float64(pr.totalWork)
	plussesNeeded := Round2Int(float64(pr.totalPlusses) * percentDone)
	if plussesNeeded > pr.totalPlusses {
		plussesNeeded = pr.totalPlusses
	}
	fmt.Printf("\r%s: [", pr.title)
	plussesPrinted := 0
	for plussesPrinted < plussesNeeded {
		fmt.Print(COMPLETED)
		plussesPrinted++
	}
	spacesRemaining := pr.totalPlusses - plussesPrinted
	for spacesRemaining > 0 {
		fmt.Printf(REMAINING)
		spacesRemaining--
	}
	fmt.Print("]")

	// Update elapsed time and estimated time to completion
	elapsedTime := time.Since(pr.startTime)
	seconds := elapsedTime.Seconds()
	estRemaining := seconds/percentDone - seconds
	if percentDone == 1.0 {
		fmt.Printf(" (%.1fs)       ", seconds)
	} else {
		fmt.Printf(" (%.1fs|%.1fs)  ", seconds, math.Max(0.0, estRemaining))
	}
}

func (pr *ProgressReporter) Done() {
	fmt.Println("")
}
