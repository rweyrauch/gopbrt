package pbrt

import (
	"math"
)

type Material interface {
	GetBSDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSDF
	GetBSSRDF(dg, dgs *DifferentialGeometry, arena *MemoryArena) *BSSRDF
}

func Bump(texture *TextureFloat, dg, dgs *DifferentialGeometry) *DifferentialGeometry {
	// Compute offset positions and evaluate displacement texture
	dgEval := *dgs

	// Shift _dgEval_ _du_ in the $u$ direction
	du := 0.5 * (math.Abs(dgs.dudx) + math.Abs(dgs.dudy))
	if du == 0.0 {
		du = 0.01
	}
	dgEval.p = dgs.p.Add(dgs.dpdu.Scale(du))
	dgEval.u = dgs.u + du
	dgEval.nn = NormalizeNormal(CreateNormalFromVector(CrossVector(dgs.dpdu, dgs.dpdv)).Add(dgs.dndu.Scale(du)))

	uDisplace := texture.Evaluate(&dgEval)

	// Shift _dgEval_ _dv_ in the $v$ direction
	dv := 0.5 * (math.Abs(dgs.dvdx) + math.Abs(dgs.dvdy))
	if dv == 0.0 {
		dv = 0.01
	}
	dgEval.p = dgs.p.Add(dgs.dpdv.Scale(dv))
	dgEval.u = dgs.u
	dgEval.v = dgs.v + dv
	dgEval.nn = NormalizeNormal(CreateNormalFromVector(CrossVector(dgs.dpdu, dgs.dpdv)).Add(dgs.dndv.Scale(dv)))
	vDisplace := texture.Evaluate(&dgEval)
	displace := texture.Evaluate(dgs)

	// Compute bump-mapped differential geometry
	dgBump := *dgs
	dgBump.dpdu = dgs.dpdu.Add(CreateVectorFromNormal(dgs.nn).Scale((uDisplace - displace) / du)).Add(CreateVectorFromNormal(dgs.dndu).Scale(displace))
	dgBump.dpdv = dgs.dpdv.Add(CreateVectorFromNormal(dgs.nn).Scale((vDisplace - displace) / dv)).Add(CreateVectorFromNormal(dgs.dndv).Scale(displace))
	dgBump.nn = CreateNormalFromVector(NormalizeVector(CrossVector(dgBump.dpdu, dgBump.dpdv)))
	if Xor(dgs.shape.ReverseOrientation(), dgs.shape.TransformSwapsHandedness()) {
		dgBump.nn = dgBump.nn.Negate()
	}
	// Orient shading normal to match geometric normal
	dgBump.nn = FaceforwardNormal(dgBump.nn, dg.nn)

	return &dgBump
}
