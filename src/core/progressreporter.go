package core

import (
	"fmt"
	"time"
	"math"
)

func TerminalWidth() int { return 80 }

type ProgressReporter struct{
    totalWork int
    workDone, totalPlusses int
    title string
    startTime time.Time
}

func NewProgressReporter(totalwork int, title string, barLength int) *ProgressReporter {
	pr := new(ProgressReporter)
	pr.totalWork = totalwork
	
    if barLength <= 0 {
        barLength = TerminalWidth() - 28
	}
	pr.title = title
    pr.totalPlusses = Maxi(2, barLength - len(title))
    pr.workDone = 0
    pr.startTime = time.Now()
    
	return pr
}
const (
	COMPLETED = "+"
	REMAINING = " "
)

func (pr *ProgressReporter) Update(num int) {
	if num == 0 || options.Quiet { return }
	pr.workDone += num
	percentDone := float64(pr.workDone) / float64(pr.totalWork)
	plussesNeeded := Round2Int(float64(pr.totalPlusses) * percentDone)
	if plussesNeeded > pr.totalPlusses { plussesNeeded = pr.totalPlusses }
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
    estRemaining := seconds / percentDone - seconds
    if percentDone == 1.0 {
        fmt.Printf(" (%.1fs)       ", seconds)
    } else {
        fmt.Printf(" (%.1fs|%.1fs)  ", seconds, math.Max(0.0, estRemaining))
    }    
}

func (pr *ProgressReporter) Done()   {
}
