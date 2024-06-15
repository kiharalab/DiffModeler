#!/usr/bin/env python
#Geo3DUtils.py

import string
import sys
import os
import copy
import math
import gzip

def vecteur(x1,y1,z1,x2,y2,z2):
	return x2-x1, y2-y1, z2-z1

def distance(x1,y1,z1,x2,y2,z2):
	dx = x2-x1
	dy = y2-y1
	dz = z2-z1
	d = dx*dx+dy*dy+dz*dz
	return math.sqrt(d)

def mixtproduct(x1,y1,z1,x2,y2,z2,x3,y3,z3):
	x = y1 * z2 - z1 * y2
	y = z1 * x2 - x1 * z2
	z = x1 * y2 - y1 * x2
	n = math.sqrt(x*x + y*y + z*z)
	x = x / n
	y = y / n
	z = z / n
	
	return x * x3 + y * y3 + z * z3

def RTOD(x):
	return x * 180.0 / 3.14159265358979323846
	
def DTOR(x):
	return x * 3.14159265358979323846 / 180.0
	
def dihedral(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4):

#	print x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
	ab_x = (x2 - x1)
	ab_y = (y2 - y1)
	ab_z = (z2 - z1)
	bc_x = (x3 - x2)
	bc_y = (y3 - y2)
	bc_z = (z3 - z2)
	cd_x = (x4 - x3)
	cd_y = (y4 - y3)
	cd_z = (z4 - z3)
	
	d012  = ab_x * bc_x + ab_y * bc_y + ab_z * bc_z
	d123  = cd_x * bc_x + cd_y * bc_y + cd_z * bc_z
	d0123 = ab_x * cd_x + ab_y * cd_y + ab_z * cd_z
	
	d01   = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z
	d12   = bc_x * bc_x + bc_y * bc_y + bc_z * bc_z
	d23   = cd_x * cd_x + cd_y * cd_y + cd_z * cd_z
	
	num = d012 * d123 - d12 * d0123
	den = (d01*d12 - d012*d012)*(d12*d23 - d123*d123)
	arccos = num / math.sqrt(den)
	
	if arccos > 1.:
		arccos = 1.
			
	if arccos < -1.:
		arccos = -1.

	RS = math.acos(arccos)
  
	RS1 = cd_x * (ab_y * bc_z - ab_z * bc_y) + \
	      cd_y * (bc_x * ab_z - ab_x * bc_z) + \
	      cd_z * (ab_x * bc_y - ab_y * bc_x)
  
	if RS1 > 0.:
		return RTOD(RS)
	else:
		return - RTOD(RS)

## def oneHMMGeo(theCAs, aCA):
## 	CA1x, CA1y, CA1z = atmCrds(theCAs[aCA])
## 	CA2x, CA2y, CA2z = atmCrds(theCAs[aCA+1])
## 	CA3x, CA3y, CA3z = atmCrds(theCAs[aCA+2])
## 	CA4x, CA4y, CA4z = atmCrds(theCAs[aCA+3])
## 	d1 = distance(CA1x, CA1y, CA1z, CA3x, CA3y, CA3z)
## 	d2 = distance(CA1x, CA1y, CA1z, CA4x, CA4y, CA4z)
## 	d3 = distance(CA2x, CA2y, CA2z, CA4x, CA4y, CA4z)
## 	x1, y1, z1 = vecteur(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
## 	x2, y2, z2 = vecteur(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
## 	x3, y3, z3 = vecteur(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
## 	d4 = mixtproduct(x1, y1, z1, x2, y2, z2, x3, y3, z3)
## 	d5 = distance(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
## 	d6 = distance(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
## 	d7 = distance(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
## 	return d1,d2,d3,d4,d5,d6,d7

## def HMMGeo(theCAs, theId):
## 	for aCA in range(0,len(theCAs)-3):
## 		d1,d2,d3,d4,d5,d6,d7 = oneHMMGeo(theCAs, aCA)
## 		print "%s %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %3d %s" % (resNum(theCAs[aCA]), d1,d2,d3,d4,d5,d6,d7, len(theCAs)-3, theId)

## PDBDATADIR = "/raid5/HMM/data/"

## if __name__=='__main__':
## 	useChain = ' '
## 	if len(sys.argv[1]) == 5:
## 		useChain = sys.argv[1][4]
## 	theTrace = pdbTrace(PDBDATADIR+sys.argv[1]+".pdb", useChain, 1)
## 	HMMGeo(theTrace, sys.argv[1])

