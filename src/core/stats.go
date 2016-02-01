package core

import (
)

type (
	
	StatsAccumulator struct {
		counters map[string]int64
		timers map[string]int64
	}
)

var (
		
)

func (stats *StatsAccumulator) ReportCounter(name string, val int64) {
	stats.counters[name] += val
}

func (stats *StatsAccumulator) ReportTimer(name string, val int64) {
	stats.timers[name] += val
}