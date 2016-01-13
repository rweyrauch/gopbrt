package pbrt

import (
)

type Sphere struct {
    ShapeData	
}

func CreateSphereShape(o2w, w2o *Transform, reverseOrientation bool, params *ParamSet) *Sphere {
    return nil
}


func (c *Sphere) ObjectBound() *BBox{
    return nil
}

func (c *Sphere) 	WorldBound() *BBox {
    return nil
}
func (c *Sphere) 	CanIntersect() bool {
    return false
}
func (c *Sphere) 	Refine() (refined []*Shape) {
    return nil
}
func (c *Sphere) 	Intersect(ray *Ray) (hit bool, tHit, rayEpsilon float64, dg *DifferentialGeometry) {
    return false, 0.0, 0.0, nil
}
func (c *Sphere) 	IntersectP(ray *Ray) bool {
    return false
}
func (c *Sphere) 	GetShadingGeometry(obj2world *Transform, dg *DifferentialGeometry) *DifferentialGeometry {
    return nil
}
func (c *Sphere) 	Area() float64 {
    return 0.0
}
func (c *Sphere) 	Sample(u1, u2 float64) (*Point, *Normal) {
    return nil, nil
}
func (c *Sphere) 	Pdf(pshape *Point) float64 {
    return 0.0
}
func (c *Sphere) 	SampleAt(p *Point, u1, u2 float64) (*Point, *Normal) {
    return nil, nil
}
func (c *Sphere) 	Pdf2(p *Point, wi *Vector) float64 { 
    return 0.0
    } 
func (c *Sphere) 	ObjectToWorld() *Transform { 
    return c.objectToWorld
     }
func (c *Sphere) 	WorldToObject() *Transform { 
    return c.worldToObject 
    }
func (c *Sphere) 	ReverseOrientation() bool { 
    return c.reverseOrientation 
    }
func (c *Sphere) 	TransformSwapsHandedness() bool { 
    return c.transformSwapsHandedness 
    }
func (c *Sphere) 	ShapeId() uint32 { 
    return c.shapeId 
    }
