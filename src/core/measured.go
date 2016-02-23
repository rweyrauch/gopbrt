package core

import (
	"os"
	"math"
	"path/filepath"
	"encoding/binary"
	"strings"
)

/*
  File format descriptions:

  -- Irregularly Sampled Isotropic BRDF --

  This is the format of the BRDFs in the scenes/brdfs/ folder of the pbrt
  distribution.  This is a simple text format of numbers in a particular
  format; the hash character # is used to denote a comment that continues
  to the end of the current line.

  The first number in the file gives the number of wavelengths at which the
  reflection data was measured, numWls.  This is followed by numWls values
  that give the frequency in nm of these wavelengths.  Each BRDF
  measurement is represented by 4+numWls values.  The first two give the
  (theta,phi) angles of the incident illumination direction, the next two
  give (theta,phi) for the measured reflection direction, and the following
  numWls give the spectral coefficients for the measurement at the
  wavelength specified at the start of the file.


  -- Regular Half-Angle BRDF --
  This is the file format used in the MERL BRDF database; see http://merl.com/brdf.

  This file format is a binary format, with numbers encoded in low-endian
  form.  It represents a regular 3D tabularization of BRDF samples in RGB
  color where the dimensions indexed over are (delta phi, delta theta,
  sqrt(theta_h)).  Here, theta_h is the angle between the halfangle vector
  and the normal, and delta theta and delta phi are the offset in theta and
  phi of one of the two directions.  (Note that the offset would be the
  same for the other direction, since it's from the half-angle vector.)

  The starts with three 32-bit integers, giving the resolution of the
  overall table.  It then containes a number of samples equal to the
  product of those three integers, times 3 (for RGB).  Samples are laid out
  with delta phi the minor index, then delta theta, then sqrt(theta_h) as
  the major index.

  In the file each sample should be scaled by RGB(1500,1500,1500/1.6) of
  the original measurement.  (In order words, the sample values are scaled
  by the inverse of that as they are read in.  
*/

func CreateMeasuredMaterial(xform *Transform, mp *TextureParams) *MeasuredMaterial {
    bumpMap := mp.GetFloatTextureOrNil("bumpmap")
    return NewMeasuredMaterial(mp.FindFilename("filename", ""), bumpMap)
}

var loadedRegularHalfangle map[string][]float64
var loadedThetaPhi map[string]*KdTree

func init() {
	loadedThetaPhi = make(map[string]*KdTree)
	loadedRegularHalfangle = make(map[string][]float64)
}

func (s *IrregIsotropicBRDFSample) Location() *Point {
	return &s.p
}

func NewMeasuredMaterial(filename string, bump TextureFloat) *MeasuredMaterial {
	material := new(MeasuredMaterial)
	material.bumpMap = bump
	material.regularHalfangleData = nil
	material.thetaPhiData = nil
	
	suffix := strings.ToLower(filepath.Ext(filename))
	if len(suffix) == 0 {
        Error("No suffix in measured BRDF filename \"%s\". Can't determine file type (.brdf / .merl)", filename)	
	} else if strings.Compare(suffix, ".brdf") == 0 {
        // Load $(\theta, \phi)$ measured BRDF data
		phiData := loadedThetaPhi[filename]
		if phiData != nil {
			material.thetaPhiData = phiData
			return material
		}
			
		ok, values := ReadFloatFile(filename)
		if !ok {
            Error("Unable to read BRDF data from file \"%s\"", filename)
            return material
		}
		pos := 0
		numWls := int(values[pos]); pos++
		if (len(values)-1-numWls) % (4 + numWls) != 0 {
           Error("Excess or insufficient data in theta, phi BRDF file \"%s\"", filename)
           return material
		}
		
		wls := make([]float64, 0, numWls)
		for i := 0; i < numWls; i++ {
			wls = append(wls, values[pos])
			pos++
		}
		
        //var bbox BBox
        samples := make([]NodeData, 0, len(values)-numWls)
        for pos < len(values) {
            thetai := values[pos]; pos++
            phii := values[pos]; pos++
            thetao := values[pos]; pos++
            phio := values[pos]; pos++
            wo := SphericalDirection(math.Sin(thetao), math.Cos(thetao), phio)
            wi := SphericalDirection(math.Sin(thetai), math.Cos(thetai), phii)
            s := SpectrumFromSampled(wls, values[pos:pos+numWls])
            pos += numWls
            p := BRDFRemap(wo, wi)
            samples = append(samples, &IrregIsotropicBRDFSample{*p, *s})
            //bbox = *UnionBBoxPoint(&bbox, p)
        }
        material.thetaPhiData = NewKdTree(samples)
        loadedThetaPhi[filename] = material.thetaPhiData
	} else {
        // Load RegularHalfangle BRDF Data
        material.nThetaH = 90
        material.nThetaD = 90
        material.nPhiD = 180
        
        angleData := loadedRegularHalfangle[filename]
        if angleData != nil {
            material.regularHalfangleData = angleData
            return material
        }
        
        fi, err := os.Open(filename)
        defer fi.Close()
        if err != nil {
            Error("Unable to open BRDF data file \"%s\"", filename)
            return material
        }
        
        var dims [3]int
        err = binary.Read(fi, binary.LittleEndian, dims)        
        if err != nil {
            Error("Premature end-of-file in measured BRDF data file \"%s\"", filename)
            return material
        }
        n := dims[0] * dims[1] * dims[2]
        if n != material.nThetaH * material.nThetaD * material.nPhiD {
            Error("Dimensions don't match.")
            return material
        }
        
        material.regularHalfangleData = make([]float64, 3*n, 3*n)
        
        chunkSize := 2*material.nPhiD
        tmp := make([]float64, chunkSize, chunkSize)
        nChunks := n / chunkSize
        Assert((n % chunkSize) == 0)
        scales := [3]float64{ 1.0/1500.0, 1.15/1500.0, 1.66/1500.0 }
        for c := 0; c < 3; c++ {
            offset := 0
            for i := 0; i < nChunks; i++ {
                err = binary.Read(fi, binary.LittleEndian, tmp)
            	if err != nil {
                    Error("Premature end-of-file in measured BRDF data file \"%s\"", filename)
                    material.regularHalfangleData = nil
                    return material
                }
                for j := 0; j < chunkSize; j++ {
                    material.regularHalfangleData[3 * offset + c] = math.Max(0.0, tmp[j] * scales[c])
                    offset++
				}                    
            }
        }
        
        loadedRegularHalfangle[filename] = material.regularHalfangleData		
	}
	return material	
}

func (m *MeasuredMaterial) GetBSDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSDF {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    var dgs *DifferentialGeometry
    if m.bumpMap != nil {
        dgs = Bump(m.bumpMap, dgGeom, dgShading)
    } else {
        dgs = dgShading
    }    
	bsdf := NewBSDF(dgs, dgGeom.nn, 1.0)
    if m.regularHalfangleData != nil {
        bsdf.Add(NewRegularHalfangleBRDF(m.regularHalfangleData, m.nThetaH, m.nThetaD, m.nPhiD))
    } else if m.thetaPhiData != nil {
        bsdf.Add(NewIrregIsotropicBRDF(m.thetaPhiData))
    }    
    return bsdf
}
func (m *MeasuredMaterial) GetBSSRDF(dgGeom, dgShading *DifferentialGeometry, arena *MemoryArena) *BSSRDF {
	return nil
}
