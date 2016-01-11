package pbrt

import (
	"math/rand"
)

type RNG struct {
	rng *rand.Rand
}

func CreateRNG(seed int64) {
	rng := new(RNG)
	rng.rand = rand.New(math.rand.NewSource(seed))
	return rng
}

func (r *RNG) Seed(seed int64) {
	r.rng.Seed(seed)
}

// generates a random number on [0,1)-real-interval
func (r *RNG) RandomFloat() float64 {
	return r.rng.Float64()
}

func (r *RNG) RandomUInt() uint32 {
	return r.rng.Uint32()
}