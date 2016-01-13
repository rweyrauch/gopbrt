package shapes

import (
    "github.com/rweyrauch/gopbrt/src/core"
)

type Cone struct {
	pbrt.ShapeData
}

func CreateConeShape(o2w, w2o *pbrt.Transform, reverseOrientation bool, params *pbrt.ParamSet) *Cone {
    return nil
}

func (c *Cone) ObjectBound() *pbrt.BBox{
    return nil
}

func (c *Cone) 	WorldBound() *pbrt.BBox {
    return nil
}
func (c *Cone) 	CanIntersect() bool {
    return false
}
func (c *Cone) 	Refine() (refined []*pbrt.Shape) {
    return nil
}
func (c *Cone) 	Intersect(ray *pbrt.Ray) (hit bool, tHit, rayEpsilon float64, dg *pbrt.DifferentialGeometry) {
    return false, 0.0, 0.0, nil
}
func (c *Cone) 	IntersectP(ray *pbrt.Ray) bool {
    return false
}
func (c *Cone) 	GetShadingGeometry(obj2world *pbrt.Transform, dg *pbrt.DifferentialGeometry) *pbrt.DifferentialGeometry {
    return nil
}
func (c *Cone) 	Area() float64 {
    return 0.0
}
func (c *Cone) 	Sample(u1, u2 float64) (*pbrt.Point, *pbrt.Normal) {
    return nil, nil
}
func (c *Cone) 	Pdf(pshape *pbrt.Point) float64 {
    return 0.0
}
func (c *Cone) 	SampleAt(p *pbrt.Point, u1, u2 float64) (*pbrt.Point, *pbrt.Normal) {
    return nil, nil
}
func (c *Cone) 	Pdf2(p *pbrt.Point, wi *pbrt.Vector) float64 { 
    return 0.0
    } 
func (c *Cone) 	ObjectToWorld() *pbrt.Transform { 
    return c.objectToWorld
     }
func (c *Cone) 	WorldToObject() *pbrt.Transform { 
    return c.worldToObject 
    }
func (c *Cone) 	ReverseOrientation() bool { 
    return c.reverseOrientation 
    }
func (c *Cone) 	TransformSwapsHandedness() bool { 
    return c.transformSwapsHandedness 
    }
func (c *Cone) 	ShapeId() uint32 { 
    return c.shapeId 
    }
