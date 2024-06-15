#!/usr/bin/env python

#
# A parser for PDB files
#
# Written (2001-2003) by P. Tuffery, INSERM, France
# Contributions by R. Gautier, J. Maupetit, J. Herisson
# Version: 6.8 (2008 september)
#
# No warranty of any kind is provided
# This is free software. You can use it, modify it, distribute it
# but its origin (i.e. this present text) must remain clearly
# stated.
#
# Thanks for any feedback related to any bug fix, or any improvement.
# 
#
# Main classes:
#   PDB (for one PDB file)
#   PDBList (or a collection of PDBFiles)
#
"""
Classe PDB pour parser un PDB

Developpements par P. Tuffery (2004-) et
aussi R. Gautier, J. Maupetit, J. Herisson

L'idee est d'avoir une gestion aisee des PDBs.

Ex:

x = PDB("1tim")
y = PDB("1timA")

# 1er residu de x
x[0]

# Chaine A de x
x["A"]

#1er atome du 2eme residu de x
x[0][1]


Le schema des classes est:

PDBLine (une ligne d'un fichier PDB)
atmLine (une ligne ATOM/HETATM d'un fichier PDB)
atmList (une serie de lignes ATOM/HETATM d'un fichier PDB (ex: une serie
de ligne atomes))
residu (un residu d'un fichier (information atomique))
PDB (un PDB)

PDBList (pas une classe pour l'instant) gere un collection de PDBs.


protein : Classe pour gerer un PDB de type proteine
(n'est pas au niveau actuellement (adaptations en cours), mais
est partiellement fonctionelle)

NOTE:
(Le developpement a ete incremental, et un lifting semble
desormais necessaire.)

"""

import string
import sys
import os
import copy
import math
import gzip
import types
import popen2
import urllib
import urllib2

## sys.path.append('/home/raid5/PyTools/Classes/')
## sys.path.append('/home/tuffery/proteineDBTools/PyScripts/')
sys.path.append("/data/PyTools/Classes/")

from FileBasics import *
## from Lines import *
from Geo3DUtils import *
from html2text import *


GDFLTPDBDIR = "/home/raid5/pdb/data/structures/"
GDFLTPDBDIR = "/data/pdb/data/structures/"
GDFLTSCOPDIR = "/data/banks/Astral/current/"
GDFLTCATHDIR = "/data/banks/CATH/current/pdb"

AA1 = "ACDEFGHIKLMNPQRSTVWY"
# HSE, HSP, HSD rajoute pour les histidines protonees par charmm (stef)
AA3 = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HSE", "HSP", "HSD"]
AA1seq = "ACDEFGHIKLMNPQRSTVWYXXXSXWMCXWYMDPECXXYAMMSCMA"
AA3STRICT = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

# HSE, HSP, HSD rajoute pour les histidines protonees par charmm (stef)
AA3new = ['PAQ', 'AGM', 'PR3', 'DOH', 'CCS', 'GSC', 'GHG', 'OAS', 'MIS', 'SIN', 'TPL', 'SAC', '4HT', 'FGP', 'HSO', 'LYZ', 'FGL', 'PRS', 'DCY', 'LYM', 'GPL', 'PYX', 'PCC', 'EHP', 'CHG', 'TPO', 'DAS', 'AYA', 'TYN', 'SVA', 'SCY', 'BNN', '5HP', 'HAR', 'IAS', 'SNC', 'AHB', 'PTR', 'PHI', 'NPH', 'PHL', 'SNN', 'A66', 'TYB', 'PHD', 'MAA', 'APN', 'TYY', 'TYT', 'TIH', 'TRG', 'CXM', 'DIV', 'TYS', 'DTH', 'MLE', 'CME', 'SHR', 'OCY', 'DTY', '2AS', 'AEI', 'DTR', 'OCS', 'CMT', 'BET', 'NLP', 'LLY', 'SCH', 'CEA', 'LLP', 'TRF', 'HMR', 'TYI', 'TRO', 'NLE', 'BMT', 'BUC', 'PEC', 'BUG', 'SCS', 'NLN', 'MHO', 'CSO', 'FTR', 'DLE', 'TRN', 'CSE', 'CSD', 'OMT', 'CSA', 'DSP', 'CSB', 'DSN', 'SHC', 'CSX', 'YCM', 'CSZ', 'TRQ', 'CSW', 'EFC', 'CSP', 'CSS', 'CSR', 'CZZ', 'MSO', 'BTR', 'HLU', 'MGN', 'HTI', 'TYQ', '4IN', 'M3L', 'C5C', 'HTR', 'MPQ', 'KCX', 'GLH', 'DIL', 'ACA', 'NEM', '5CS', 'LYX', 'DVA', 'ACL', 'GLX', 'MLZ', 'GLZ', 'SME', 'SMC', 'DLY', 'NEP', 'BCS', 'ASQ', 'SET', 'SEP', 'ASX', 'DGN', 'DGL', 'MHS', 'SEG', 'ASB', 'ASA', 'SEC', 'SEB', 'ASK', 'GGL', 'ASI', 'SEL', 'CGU', 'C6C', 'ASL', 'LTR', 'CLD', 'CLE', 'GMA', '1LU', 'CLB', 'MVA', 'S1H', 'DNP', 'SAR', 'FME', 'ALO', 'ALM', 'LEF', 'MEN', 'TPQ', 'NMC', 'SBD', 'ALY', 'MME', 'GL3', 'ALS', 'SBL', '2MR', 'CAY', '3AH', 'DPR', 'CAS', 'NC1', 'HYP', 'FLA', 'LCX', 'MSE', 'IYR', 'DPN', 'BAL', 'CAF', 'MSA', 'AIB', 'HIP', 'CYQ', 'PCA', 'DAL', 'BFD', 'DAH', 'HIC', 'CYG', 'DAR', 'CYD', 'IIL', 'CYM', 'CYL', 'CY3', 'CY1', 'HAC', '143', 'DHI', 'CY4', 'YOF', 'HPQ', 'SOC', 'DHA', '2LU', 'MLY', 'TRW', 'STY', 'MCL', 'BHD', 'NRQ', 'ARM', 'PRR', 'ARO', "5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HSE", "HSP", "HSD"]

# HSE, HSP, HSD rajoute pour les histidines protonees par charmm (stef)
dico_AA = {
 "ALA": 'A',
 "CYS": 'C',
 "ASP": "D",
 "GLU": 'E',
 "PHE": 'F',
 "GLY": 'G',
 "HIS": 'H',
 "ILE": 'I',
 "LYS": 'K',
 "LEU": 'L',
 "MET": 'M',
 "ASN": 'N',
 "PRO": 'P',
 "GLN": 'Q',
 "ARG": 'R',
 "SER": 'S',
 "THR": 'T',
 "VAL": 'V',
 "TRP": 'W',
 "TYR": 'Y',
 'PAQ': 'Y', 
 'AGM': 'R', 
 'PR3': 'C', 
 'DOH': 'D', 
 'CCS': 'C', 
 'GSC': 'G', 
 'GHG': 'Q', 
 'OAS': 'S', 
 'MIS': 'S', 
 'SIN': 'D', 
 'TPL': 'W', 
 'SAC': 'S', 
 '4HT': 'W', 
 'FGP': 'C', 
 'HSO': 'H', 
 'LYZ': 'K', 
 'FGL': 'S', 
 'PRS': 'P', 
 'DCY': 'C', 
 'LYM': 'K', 
 'GPL': 'K', 
 'PYX': 'C', 
 'PCC': 'P', 
 'EHP': 'F', 
 'CHG': 'A', 
 'TPO': 'T', 
 'DAS': 'D', 
 'AYA': 'A', 'TYN': 'Y', 'SVA': 'S', 'SCY': 'C', 'BNN': 'A', '5HP': 'E', 'HAR': 'R', 'IAS': 'D', 'SNC': 'C', 'AHB': 'N', 'PTR': 'Y', 'PHI': 'F', 'NPH': 'C', 'PHL': 'F', 'SNN': 'D', 'A66': 'A', 'TYB': 'Y', 'PHD': 'D', 'MAA': 'A', 'APN': 'A', 'TYY': 'Y', 'TYT': 'Y', 'TIH': 'A', 'TRG': 'K', 'CXM': 'M', 'DIV': 'V', 'TYS': 'Y', 'DTH': 'T', 'MLE': 'L', 'CME': 'C', 'SHR': 'K', 'OCY': 'C', 'DTY': 'Y', '2AS': 'D', 'AEI': 'T', 'DTR': 'W', 'OCS': 'C', 'CMT': 'C', 'BET': 'G', 'NLP': 'L', 'LLY': 'K', 'SCH': 'C', 'CEA': 'C', 'LLP': 'K', 'TRF': 'W', 'HMR': 'R', 'TYI': 'Y', 'TRO': 'W', 'NLE': 'L', 'BMT': 'T', 'BUC': 'C', 'PEC': 'C', 'BUG': 'L', 'SCS': 'C', 'NLN': 'L', 'MHO': 'M', 'CSO': 'C', 'FTR': 'W', 'DLE': 'L', 'TRN': 'W', 'CSE': 'C', 'CSD': 'A', 'OMT': 'M', 'CSA': 'C', 'DSP': 'D', 'CSB': 'C', 'DSN': 'S', 'SHC': 'C', 'CSX': 'C', 'YCM': 'C', 'CSZ': 'C', 'TRQ': 'W', 'CSW': 'C', 'EFC': 'C', 'CSP': 'C', 'CSS': 'C', 'CSR': 'C', 'CZZ': 'C', 'MSO': 'M', 'BTR': 'W', 'HLU': 'L', 'MGN': 'Q', 'HTI': 'C', 'TYQ': 'Y', '4IN': 'W', 'M3L': 'K', 'C5C': 'C', 'HTR': 'W', 'MPQ': 'G', 'KCX': 'K', 'GLH': 'E', 'DIL': 'I', 'ACA': 'A', 'NEM': 'H', '5CS': 'C', 'LYX': 'K', 'DVA': 'V', 'ACL': 'R', 'GLX': 'Z', 'MLZ': 'K', 'GLZ': 'G', 'SME': 'M', 'SMC': 'C', 'DLY': 'K', 'NEP': 'H', 'BCS': 'C', 'ASQ': 'D', 'SET': 'S', 'SEP': 'S', 'ASX': 'B', 'DGN': 'Q', 'DGL': 'E', 'MHS': 'H', 'SEG': 'A', 'ASB': 'D', 'ASA': 'D', 'SEC': 'C', 'SEB': 'S', 'ASK': 'D', 'GGL': 'E', 'ASI': 'N', 'SEL': 'S', 'CGU': 'E', 'C6C': 'C', 'ASL': 'D', 'LTR': 'W', 'CLD': 'S', 'CLE': 'L', 'GMA': 'E', '1LU': 'L', 'CLB': 'S', 'MVA': 'V', 'S1H': 'S', 'DNP': 'A', 'SAR': 'G', 'FME': 'M', 'ALO': 'T', 'ALM': 'A', 'LEF': 'L', 'MEN': 'N', 'TPQ': 'Y', 'NMC': 'G', 'SBD': 'S', 'ALY': 'K', 'MME': 'M', 'GL3': 'G', 'ALS': 'C', 'SBL': 'S', '2MR': 'R', 'CAY': 'C', '3AH': 'H', 'DPR': 'P', 'CAS': 'C', 'NC1': 'S', 'HYP': 'P', 'FLA': 'A', 'LCX': 'K', 'MSE': 'M', 'IYR': 'Y', 'DPN': 'F', 'BAL': 'A', 'CAF': 'C', 'MSA': 'G', 'AIB': 'A', 'HIP': 'H', 'CYQ': 'C', 'PCA': 'E', 'DAL': 'A', 'BFD': 'D', 'DAH': 'F', 'HIC': 'H', 'CYG': 'C', 'DAR': 'R', 'CYD': 'C', 'IIL': 'I', 'CYM': 'C', 'CYL': 'C', 'CY3': 'C', 'CY1': 'C', 'HAC': 'A', '143': 'C', 'DHI': 'H', 'CY4': 'C', 'YOF': 'Y', 'HPQ': 'F', 'SOC': 'C', 'DHA': 'A', '2LU': 'L', 'MLY': 'K', 'TRW': 'W', 'STY': 'Y', 'MCL': 'K', 'BHD': 'D', 'NRQ': 'Y', 'ARM': 'R', 'PRR': 'A', 'ARO': 'R', 'HSE': 'H', 'HSP': 'H', 'HSD': 'H'}

RNA3 = ["U"]
DNA3 = ["A","T","G","C"]
SOLV = ["HOH","H2O","WAT","DOD"]

# BBATMS = ["N","CA","C","O","OXT"]
BBATMS = ["N","CA","C","O","OXT"]
SCATMS = ["-","N","CA","C","O","OXT"]
NCHIS  = [0,1,2,3,2,0,2,2,4,2,3,2,0,3,5,1,1,1,2,2]

CHIATMS = [ \
	[], \
	[["N","CA","CB","SG"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","ND1"]], \
	[["N","CA","CB","CG1"],["CA","CB","CG1","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","CE"],["CG","CD","CE","NZ"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","SD"], \
	 ["CB","CG","SD","CE"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","NE"],["CG","CD","NE","CZ"], \
	 ["CD","NE","CZ","NH1"]], \
	[["N","CA","CB","OG"]], \
	[["N","CA","CB","OG1"]], \
	[["N","CA","CB","CG1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]]] 

AASC=[["CB"], \
      ["CB","SG"], \
      ["CB","CG","OD1","OD2"], \
      ["CB","CG","CD","OE1","OE2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ"], \
      [],["CB","CG","ND1","CD2","CE1","NE2"], \
      ["CB","CG1","CG2","CD1"], \
      ["CB","CG","CD","CE","NZ"], \
      ["CD","CG","CD1","CD2"], \
      ["CB","CG","SD","CE"], \
      ["CB","CG","OD1","ND2"], \
      ["CB","CG","CD"], \
      ["CB","CG","CD","OE1","NE2"], \
      ["CB","CG","CD","NE","CZ","NH1","NH2"], \
      ["CB","OG"], \
      ["CB","OG","OG1","CG2"], \
      ["CB","CG1","CG2"], \
      ["CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]]

AABB=["N","CA","C","O"]

#GBINPATH="/data/bin/"
GBINPATH="/data/bin/"
GBINPATH="/home/tintin/tuffery/bin/"
#GHMMPATH="/data/HMM/models/HMM1/"
GHMMPATH="/data/HMM/models/HMM1/"

HNames = {
	"ALA" : ["HB1","HB2","HB3"],
	"CYS" : ["HB1","HB2","HG"],
	"ASP" : ["HB1","HB2"],
	"GLU" : ["HB1","HB2","HG1","HG2"],
	"PHE" : ["HB1","HB2","HD1","HE1","HZ","HE2","HD2"],
	"GLY" : ["HA2"],
	"HIS" : ["HB1","HB2","HD2","HE1","HD1"],
	"ILE" : ["HB","HG11","HG12","HD11","HD12","HD13","HG21","HG22","HG23"],
	"LYS" : ["HB1","HB2","HG1","HG2","HD1","HD2","HE1","HE2","HZ1","HZ2","HZ3"],
	"LEU" : ["HB1","HB2","HG","HD11","HD12","HD13","HD21","HD22","HD23"],
	"MET" : ["HB1","HB2","HG1","HG2","HE1","HE2","HE3"],
	"ASN" : ["HB1","HB2","HD21","HD22"],
	"PRO" : ["HB1","HB2","HG1","HG2","HD1","HD2"],
	"GLN" : ["HB1","HB2","HG1","HG2","HE21","HE22"],
	"ARG" : ["NH1","NH2","HB1","HB2","HG1","HG2","HD1","HD2","HE","HH11","HH12","HH21","HH22"],
	"SER" : ["HB1","HB2","HG"],
	"THR" : ["HB","HG1","HG21","HG22","HG23"],
	"VAL" : ["HB","HG11","HG12","HG13","HG21","HG22","HG23"],
	"TRP" : ["HB1","HB2","HD1","HE1","HZ2","HH2","HZ3","HE3"],
	"TYR" : ["HB1","HB2","HD1","HE1","HE2","HD2","HH"],
	"BCK" : ["HA","HN","HN1","HN2","HN3"]
	}

PDBHNames = {
	"ALA" : ["1HB","2HB","3HB"],
	"CYS" : ["1HB","2HB"," HG"],
	"ASP" : ["1HB","2HB"],
	"GLU" : ["1HB","2HB","1HG","2HG"],
	"PHE" : ["1HB","2HB"," HD1"," HE1"," HZ"," HE2"," HD2"],
	"GLY" : ["2HA"],
	"HIS" : ["1HB","2HB"," HD2"," HE1"," HD1"],
	"ILE" : [" HB","1HG1","2HG1","1HD1","2HD1","3HD1","1HG2","2HG2","3HG2"],
	"LYS" : ["1HB","2HB","1HG","2HG","1HD","2HD","1HE","2HE","1HZ","2HZ","3HZ"],
	"LEU" : ["1HB","2HB"," HG","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2"],
	"MET" : ["1HB","2HB","1HG","2HG","1HE","2HE","3HE"],
	"ASN" : ["1HB","2HB","1HD2","2HD2"],
	"PRO" : ["1HB","2HB","1HG","2HG","1HD","2HD"],
	"GLN" : ["1HB","2HB","1HG","2HG","1HE2","2HE2"],
	"ARG" : ["1HB","2HB","1HG","2HG","1HD","2HD"," HE","1HH1","2HH1","1HH2","2HH2"],
	"SER" : ["1HB","2HB"," HG"],
	"THR" : [" HB"," HG1","1HG2","2HG2","3HG2"],
	"VAL" : [" HB","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2"],
	"TRP" : ["1HB","2HB"," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3"],
	"TYR" : ["1HB","2HB"," HD1"," HE1"," HE2"," HD2"," HH"],
	"BCK" : [" HA"," H","1H","2H","3H"]
	}


def clean_PDB(pdbpath):
    os.system("grep '^ATOM' %s > %s_temp"%(pdbpath,pdbpath))
    os.system("mv %s_temp %s"%(pdbpath,pdbpath))


def resType(aName):
	if aName == "":
		return "Error"
	if string.count(AA3, aName) > 0:
		return string.index(AA3, aName)
	else:
		return "Error"

def aa3Type(aName):
	if aName == "":
		return "Error"
	if string.count(AA3, aName) > 0 :
		return string.index(AA3, aName)
	else:
		return "Error"

def aa1Type(aName):
	if aName == "":
		return "Error"
	if string.count(AA1, aName) > 0:
		return string.index(AA1, aName)
	else:
		return "Error"

#
# a series of AA3 separated by blanks into aa1 string
#
def SEQREStoAA1(seqres, verbose = 0):

	seq = ""
	aList = string.split(seqres)
	for aRes in aList:
		if verbose:
			sys.stderr.write("SEQREStoAA1: %s\n" % aRes)
# 		if AA3.count(aRes) != 0:
# 			if verbose:
# 				sys.stderr.write("Found as AA3 %d\n" % AA3.index(aRes))
# 			seq = seq + AA1seq[AA3.index(aRes)]
		if dico_AA.has_key(aRes):
			if verbose:
				sys.stderr.write("Found as AA3 %d\n" % AA3.index(aRes))
			seq = seq + dico_AA[aRes]
		else:
			seq = seq + "X"
		# print seq
	return seq

#
# This will convert an alignement into a selection mask.
# s1 and s2 must be 2 strings of identical lengths.
# gaps (indels) must be represented by '-'
#
def aln2mask(s1,s2):
	res = ""
	if len(s1) != len(s2):
		return res
	for i in range(0,len(s1)):
		if s1[i] == '-':
			continue
		if s2[i] == '-':
			res = res + '-'
		else:
			res = res + s1[i]
	return res

#
# any PDB line
#
class PDBLine:
	"""
	PDBLine : basic management (mostly to access line type (ATOM, REMARK, etc) for one text line of a PDB datafile.
	"""
	def __init__(self, aLine = ""):
		self.txt = aLine

 	def __getslice__(self,ffrom=0,tto=-1):
 		return self.txt[ffrom:tto]

	def __repr__(self):
		return str(self.txt)

	def __getitem__(self,aPos):
		return self.txt[aPos]

	def __len__(self):
		return len(self.txt)

	# back to list of lines
	def flat(self):
		return self.txt

	# header de la ligne
	def header(self):
		"""
		PDBLine header, i.e. its 6 first chars
		one of:
		HEADER
		REMARK
		ATOM
		CONECT
		etc
		"""
		try:
			return string.split(self.txt[0:6])[0]
		except:
			return ""


#
# a PDB ATOM (HETATM) line
#
class atmLine(PDBLine):
	"""
	class amtLine

	This models one PDB ATOM / HETATOM line
	and its accessors
	"""
	def __init__(self, aLine = ""):
		if isinstance(aLine,atmLine):
			## print "atmLine from atmLine"
			self.txt = aLine.txt
		elif isinstance(aLine,PDBLine):
			## print "atmLine from PDBLine"
			self.txt = aLine.txt
		elif isinstance(aLine,types.StringType):
			## print "atmLine from string"
			self.txt = aLine
		else:
			self.txt = aLine

	def header(self, hdr = ""):
		"""
		header()
		PDB line HEADER
		self.header() return ATOM or HETATM
		self.header(HETATM) sets HETATM as header
		header must be string
		"""
		if hdr != "":
			self.txt = "%-.6s%s" % (hdr,self.txt[6:])
			return hdr[:6]
		try:
			return string.split(self.txt[0:6])[0]
		except:
			return ""

	def atmNum(self, anum = ""):
		"""
		atmNum()
		PDB line atom number
		self.atmNum() returns atom number (string)
		self.atmNum(anum) sets anum as atom number
		anum may be int or string
		"""
		if anum != "":
			self.txt = "%s%5d%s" % (self.txt[:6],int(anum),self.txt[11:])
			return anum
		try:
			anum=string.split(self.txt[6:11])[0]
			return anum
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return "UNK"

	def atmName(self, aname = ""):
		"""
		atmName()
		PDB line atom name
		self.atmName() returns atom name (string)
		self.atmName(aname) sets aname as atom number
		anum must be string of size 4 at max, must begin by a blank if necessary !
		"""
		if aname != "":
			self.txt = "%s%-4s%s" % (self.txt[:12],aname,self.txt[16:])
			return aname
		try:
			rnum=string.split(self.txt[12:16])[0]
			return rnum
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return "UNK"

	def alt(self, acode = ""):
		"""
		alt()
		PDB line alternate code
		self.alt() returns alternate code (1 char string)
		self.alt(acode) sets acode as alternate code
		acode must be one character
		"""
		if acode != "":
			self.txt = "%s%c%s" % (self.txt[:16],acode,self.txt[17:])
			return acode
		try:
			alt=self.txt[16]
			return alt
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return " "
		
	def resName(self, rName = ""):
		"""
		resName()
		PDB line residue name
		self.resName() returns residue name (string)
		self.resName(rName) sets residue name
		rName must be string (3 chars)
		"""
		if rName != "":
			self.txt = "%s%-.3s%s" % (self.txt[:17],rName,self.txt[20:])
			return rName[:3]
		try:
			rname=string.split(self.txt[17:20])[0]
			return rname
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return "UNK"
		
	def chnLbl(self, lbl = ""):
		"""
		chnLbl()
		PDB line chain label
		self.chnLbl() returns atom chain label (1 character string)
		self.chnLbl(lbl) sets atom chain label
		lbl must be 1 character
		"""
		if lbl != "":
			self.txt = "%s%c%s" % (self.txt[:21],lbl[0],self.txt[22:])
			return lbl
		try:
			lbl=self.txt[21]
			return lbl
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return "UNK"

	def resNum(self, rnum = ""):
		"""
		resNum()
		PDB line residue number
		self.resNum() returns residue number (1 character string)
		self.resNum(rnum) sets atom residue number
		rnum may be string or int
		"""
		if rnum != "":
			self.txt = "%s%4d%s" % (self.txt[:22],int(rnum),self.txt[26:])
			return rnum
		try:
			rnum=string.split(self.txt[22:26])[0]
			return rnum
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return "UNK"

	def icode(self, thecode = ""):
		if thecode != "":
			self.txt = "%s%c%s" % (self.txt[:26],thecode,self.txt[27:])
			return thecode
		try:
			icode=self.txt[26]
			return icode
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return " "
		
	def xyz(self):
		try:
			x=string.split(self.txt[30:38])[0]
			y=string.split(self.txt[38:46])[0]
			z=string.split(self.txt[46:54])[0]
			return float(x), float(y), float(z)
		except ValueError:
			print "Incorrect ATOM line format for:", self.txt
			return 0., 0., 0.

	def crds(self):
		return self.txt[30:54]
	
	def setcrds(self,x,y,z):
		self.txt = "%s%8.3lf%8.3lf%8.3lf%s" % (self.txt[:30], x, y, z, self.txt[54:])
		return 

	def fpt(self, aQ=""):
		if aQ != "":
			self.txt = string.replace(self.txt,"\n","")
			self.txt = "%-60s %7.3f\n" % (self.txt[:61],aQ)
			return aQ
		try:
			occ=self.txt[54:61]
			return occ
		except ValueError:
			return "       "
		
	def q(self, aQ=""):
		if aQ != "":
			self.txt = "%s%7.3f%s" % (self.txt[:54],aQ,self.txt[61:])
			return aQ
		try:
			occ=self.txt[54:61]
			return occ
		except ValueError:
			return "       "
		
	def r(self, aR=""):
		if aR != "":
			self.txt = "%s%7.3f%s" % (self.txt[:61],aR,self.txt[68:])
			return aR
		try:
			occ=self.txt[61:68]
			return occ
		except ValueError:
			return "       "
		
	def occ(self, aOcc=""):
		if aOcc != "":
			self.txt = "%s%6.2f%s" % (self.txt[:54],aOcc,self.txt[60:])
			return aOcc
		try:
			occ=self.txt[54:60]
			return occ
		except ValueError:
			return "      "
		
	def tfac(self, tFac = ""):
                if tFac != "":
                        self.txt = "%s%6.2f%s" % (self.txt[:60],float(tFac),self.txt[66:])
                        return tFac
		try:
			tfac=self.txt[60:66]
			return tfac
		except ValueError:
			return "      "
		
	def segId(self):
		try:
			segId=self.txt[72:76]
			return segId
		except ValueError:
			return "    "
		
	def ele(self):
		try:
			ele=self.txt[76:78]
			return ele
		except ValueError:
			return "  "
		
	def chrg(self):
		try:
			chrg=self.txt[78:80]
			return chrg
		except ValueError:
			return "  "


## ========================================
## A series of PDB ATOM (HETATM) lines
## Considered as a set of residues
## Tabulation of residues is achieved
##
## atmList always return atmList,
## EXCEPT for __getitem__ when requesting in 1 residue
## where it is desirable to return atmLine
##
## atom lines accessible as: x.atms
##
## With this class, we are simply manipulating text
## No semantics associated
## ========================================
class atmList(atmLine):

	# instanciate
	def __init__(self, data = "", chId = "", hetSkip = 0, verbose = 0):
		#
		# Order of parsing is important (inheritance)
		#
		
		# from PDB: just retain PDB.data field
		if isinstance(data,PDB):
			if verbose == 2:
				print "atmList from PDB"
			self.list = data.data
		# from residue: just retain residue.list field
		elif isinstance(data,residue):
			if verbose == 2:
				print "atmList from residue"
			self.list = data.atms
		# from atmList: just propagate
		elif isinstance(data,atmList):
			if verbose == 2:
				print "atmList from atmList"
			self.list = data.list
		# from atmLine: just wrap
		elif isinstance(data,atmLine):
			## We force one line as a residue
			if verbose == 2:
				print "atmList  from atmLine"
			## print data
			self.list = []
			self.list.append(data)
		# from list: suppose a list of atomic lines
		elif isinstance(data,types.ListType):
			if verbose == 2:
				print "atmList from ListType"
			## print len(data)
			self.list = []
			for aLine in data:
				self.list.append(atmLine(aLine))
			## print self.__class__
			## self.resTab(verbose)
		else:
			if verbose == 2:
				print "atmList from unknown"
			self.list  = []
## 		print "len is",len(self.list)
## 		print "len res is",len(self)

	def __len__(self):
		return len(self.list)

	# return atmList
	def __add__(self,new):
		print "__add__.atmList"
		return atmList(self.list[:] + new.list[:])

	# return sub atmList
	def __getslice__(self,ffrom,tto):
		return atmList(self.list[ffrom:tto])

	# return one atmLine
	def __getitem__(self,aPos):
		return self.list[aPos]

	# del x[i] : ready for deletions !!
	def __delitem__(self,aPos):
		# delete old series
		aDex = self.rt[aPos][0]
		##print "Removing atoms ",self.rt[aPos][0]," to ",self.rt[aPos+1][0]
		for aAtm in range(self.rt[aPos][0],self.rt[aPos+1][0]):
			del self.atms[aDex]
		self.resTab(0)

	# Managing x[i] = y : ready for mutations !!
	def __setitem__(self,aPos, new):
		del self[aPos]

		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.atms.insert(aDex,aAtm)
			aDex = aDex + 1
		self.resTab(0)

	# return string for visualization
	def __repr__(self, altLbl = "", OXTSkip = 0, HSkip = 0):
		res = ""
		## print "atmList repr"
		for aAtm in self.list:
			if altLbl != "":
				alt = aAtm.alt()
				if alt != " " and alt != altLbl:
					continue
			if OXTSkip != 0:
				if aAtm.atmName() == "OXT":
					continue
			if HSkip != 0:
				atmName = aAtm.atmName()
				if atmName[0] == "H" or (atmName[0] in  "1234" and atmName[1] == "H"):
					continue
			res = res + str(aAtm)
		return res

	# back to list of lines
	def flat(self, altLbl = "", OXTSkip = 0, PDBMac = 0, keepH = 1):
		res = []
		for aAtm in self.list:
			if altLbl != "":
				alt = aAtm.alt()
				if alt != " " and alt != altLbl:
					continue
			if PDBMac and aAtm.atmName() == "O1":
				aAtm.atmName(" O")

			if PDBMac and aAtm.atmName() == "O2":
				aAtm.atmName(" OXT")

			if PDBMac and aAtm.atmName() == "OT1":
				aAtm.atmName(" O")

			if PDBMac and aAtm.atmName() == "OT2":
				aAtm.atmName(" OXT")

			if OXTSkip != 0:
				if aAtm.atmName() == "OXT":
					continue

			if not keepH:
				atmName = aAtm.atmName()
				if atmName[0] == "H":
					continue
				if atmName[0] in "1234":
					if atmName[1] == "H":
						continue
					if atmName[1] in "1234":
						if atmName[2] == "H":
							continue
			res.append(aAtm.flat())
		return res

	# Managing x.insert(i,y) : ready for insertions !!
	def insert(self,aPos, new):
		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.list.insert(aDex,aAtm)
			aDex = aDex + 1
		self.resTab(0)



	#
	# A series of coordinates of the range
	#
	def crds(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self.list)
		res = []
		for aAtm in range(ffrom, tto):
			res.append(atmLine(self.list[aAtm]).crds())
		return res

	#
	# A series of coordinates of the range
	#
	def xyz(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self.list)
		##print "atmList xyz", ffrom , tto
		if (ffrom == 0) and (len(self.list) == 1):
			return atmLine(self.list[ffrom]).xyz()
		else:
			res = []
			for aAtm in range(ffrom, tto):
				res.append(atmLine(self.list[aAtm]).xyz())
			return res

	#
	# center of geometry of a collection of atoms
	#
	def BC(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self)
		(x,y,z) = (0.,0.,0.)
		nAtm = 0.
		for aAtm in self[ffrom:tto].atms:
			(x1,y1,z1) = aAtm.xyz()
			x = x + x1
			y = y + y1
			z = z + z1
			nAtm = nAtm + 1.
		x = x / nAtm
		y = y / nAtm
		z = z / nAtm
		return (x,y,z)

	#
	# center of geometry of a collection of atoms
	#
	def radius(self, ffrom = 0, tto = -1):
		if tto == -1:
			tto = len(self)
		(x,y,z) = self.BC(ffrom, tto)
		rs = 0.		
		for aAtm in self[ffrom:tto].atms:
			(x1,y1,z1) = aAtm.xyz()
			r = distance(x,y,z,x1,y1,z1)
			if r > rs:
				rs = r
		return rs

	def oneChis(self):

		resTpe = resType(self.list[0].resName())
		if resTpe == "Error":
			return

		res = [AA3[resTpe]]
		for aChi in CHIATMS[resTpe]:
			aAtm = self.theAtm(aChi[0])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[1])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[2])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[3])
			if aAtm == []:
				return res
			a = self.theAtm(aChi[0]).xyz()
			b = self.theAtm(aChi[1]).xyz()
			c = self.theAtm(aChi[2]).xyz()
			d = self.theAtm(aChi[3]).xyz()
			res.append(apply(dihedral,a+b+c+d))
 		return res

	def chis(self):

		res = []
		if len(self) == 1:
			res.append(self.oneChis())
			return res
		for aRes in range(0,len(self)):
			res.append(self[aRes].oneChis())
 		return res

	def outChis(self):
		chis = self.chis()
		for i in chis:
			print i[0],
			for j in i[1:]:
				print j,
			print
			
	def atmPos(self, aName):
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == aName:
				return aPos
		return None
		
	def Npos(self):
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "N":
				# print str(self[aPos])
				return aPos
		return None

	def CApos(self):
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "CA":
				# print str(self[aPos])
				return aPos
		return None

	def Cpos(self):
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "C":
				# print str(self[aPos])
				return aPos
		return None

	def Opos(self):
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "O":
				# print str(self[aPos])
				return aPos
		return None

	def PDBHNames(self):
		"""
		Setup standard PDB names for hydrogens.
		Only works for standard amino-acids.
		"""
		for aPos in range(0,len(self.list)):
			aName = self.list[aPos].atmName()
			if aName[0] == "H":
				rName = self.list[aPos].resName()
				try:
					i = HNames[rName].index(aName)
					self.list[aPos].atmName(PDBHNames[rName][i])
				except:
					try:
						i = HNames["BCK"].index(aName)
						self.list[aPos].atmName(PDBHNames["BCK"][i])
					except:
						pass

	def out(self):
		pass

	def resName(self):
## 		print self.list[0]
## 		print self.list[0].__class__
## 		print atmLine(self.list[0])
## 		print "tutu"
## 		print self.__class__, "resName",len(self.list)
## 		print self.list[0].__class__, "resName",len(self.list[0])
## 		return self.list[0].resName()
		return atmLine(self.list[0]).resName()
	
	def theAtm(self,atmName = ""):
		for aLine in self.list:
			if atmLine(aLine).atmName() == atmName:
				return atmLine(aLine)
		return []

	def isPDB(self):
		return 1
	#
	# write PDB or PDB chain(s) to file
	#
	def write(self, outName = "", label="", hetSkip = 0,verbose = 0):
		if outName == "":
			f = sys.stdout
		else:
			f = open(outName,"w")

		f.write("HEADER %s (%d residues)\n" % (label, len(self)))
		for aAtm in self.list:
			f.write("%s" % aAtm)
	# from PDB import *
	# x = protein("/home/raid5/PDB/pdb1acc.ent.gz",hetSkip=1)
	# x.frg(0).write()

	def oneHMMGeo(self, aCA):
		CA1x, CA1y, CA1z = self[aCA].xyz()
		CA2x, CA2y, CA2z = self[aCA+1].xyz()
		CA3x, CA3y, CA3z = self[aCA+2].xyz()
		CA4x, CA4y, CA4z = self[aCA+3].xyz()
		d1 = distance(CA1x, CA1y, CA1z, CA3x, CA3y, CA3z)
		d2 = distance(CA1x, CA1y, CA1z, CA4x, CA4y, CA4z)
		d3 = distance(CA2x, CA2y, CA2z, CA4x, CA4y, CA4z)
		x1, y1, z1 = vecteur(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		x2, y2, z2 = vecteur(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		x3, y3, z3 = vecteur(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		d4 = mixtproduct(x1, y1, z1, x2, y2, z2, x3, y3, z3)
		d5 = distance(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		d6 = distance(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		d7 = distance(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		return d1,d2,d3,d4,d5,d6,d7

## ========================================
## The ONE residue class
## ========================================
class residue(atmList):
	def __init__(self,data="",verbose=0):
		if data == "":
			self.atms = []
			self.type = None
			self.name = None
		else:
			if isinstance(data,residue): # residue instance
				## print "residue from residue"
				self.atms = data.atms
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data,atmList): # atmList instance
				## print "residue from atmList"
				self.atms = data
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data,atmLine): # atmLine instance
				## print "residue from atmLine"
				self.atms = atmList(data)
				## self.type = self.rType()
				## self.name = self.rName()
			else:
				## print "residue from unknown",data.__class__
				self.atms = atmList(data)
				self.type = self.rType()
				self.name = self.rName()

	def __len__(self):
		return len(self.atms)

	def __repr__(self, altCare = 0, altLbl = "", OXTCare = 0, HSkip = 0):
		OXTSkip = 0
		if OXTCare != 0:
			if self.atms.atmPos("O") != None and self.atms.atmPos("OXT") != None :
				OXTSkip = 1
		if altCare != 0:
			altLbls = self.altLbls()
			if altLbl == "":
				if altLbls != "":
					altLbl = altLbls[0]
			else:
				if string.count(altLbls, altLbl):
					pass
				else:
					altLbl = altLbls[0]	
		return self.atms.__repr__(altLbl, OXTSkip = OXTSkip, HSkip = HSkip)

	def __getslice__(self,ffrom,tto):
		## print "__getslice__"
		## return atmList(self.atms[ffrom:tto])
		if len(self.atms) == 1:
			return residue(self.atms)
		if tto > len(self.atms):
			tto = len(self.atms)
		return self.atms[ffrom:tto]

	# managing x[i]
	def __getitem__(self,aPos):
		if isinstance(aPos,types.IntType):
			## print "residue.__getitem__[",aPos,"]"
			if aPos > len(self.atms):
				return None
			elif aPos < 0:
				if aPos + len(self.atms) < 0:
					return None
				else:
					return self.atms[aPos]
			## else, we return atmList
			return self.atms[aPos]

		elif isinstance(aPos,types.StringType):
			for iAtm in range(0,len(self.atms)):
				if self.atms[iAtm].atmName() == aPos:
					return self.atms[iAtm]

	# back to list of lines
	def flat(self, altCare = 0, altLbl = "", OXTCare = 0, PDBMac = 0, keepH = 1):
		OXTSkip = 0
		if OXTCare != 0:
			if self.atms.atmPos("O") != None and self.atms.atmPos("OXT") != None :
				OXTSkip = 1
		if altCare != 0:
			altLbls = self.altLbls()
			if altLbl == "":
				if altLbls != "":
					altLbl = altLbls[0]
			else:
				if string.count(altLbls, altLbl):
					pass
				else:
					altLbl = altLbls[0]	
		return self.atms.flat(altLbl, OXTSkip = OXTSkip, PDBMac = PDBMac, keepH = keepH)
	
	def rName(self, name = "", verbose = 0):
		if name == "":
			return self.atms[0].resName()
		else:
			for atm in self.atms:
				atm.resName(name)

	def rNum(self,aNum = "", verbose = 0):
		if aNum == "":
			return self.atms[0].resNum()
		else:
			for atm in self.atms:
				atm.resNum(aNum)

	def riCode(self,icode = "",verbose = 0):
		if icode == "":
			return self.atms[0].icode()
		else:
			for atm in self.atms:
				atm.icode(icode)

	def rType(self,verbose = 0):
		aName = self.atms[0].resName()
		if AA3.count(aName) > 0:
			return "AMINO-ACID"
		elif RNA3.count(aName) > 0:
			return "RNA"
		elif DNA3.count(aName) > 0:
			return "DNA"
		elif SOLV.count(aName) > 0:
			return "SOLVENT"
		else:
			return "HETERO"

	def chnLbl(self,lbl = "", verbose = 0):
		if lbl == "":
			return self.atms[0].chnLbl()
		else:
			for atm in self.atms:
				atm.chnLbl(lbl)

 	def atmPos(self, aName):
		return self.atms.atmPos(aName)

	def hasAltAtms(self,verbose = 0):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
##			print aAtm
			
			alt = aAtm.alt()

			if alt != ' ':
				isAlt = 1
				if string.count(string.digits,aAtm.txt[12]):
					isAlt = 0
				if aAtm.txt[12] == ' ' and aAtm.txt[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = aAtm.atmName()
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "O":
					BBAltAtm = "Yes"
				else:
					SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm

	def altLbls(self,verbose = 0):
		rs = ""
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
##			print aAtm
			
			alt = aAtm.alt()

			if alt != ' ':
				if string.count(rs,alt) == 0:
					rs += alt
		return rs

	def PDBHNames(self):
		self.atms.PDBHNames()

	#
	# return a selection of atoms
	# an atmList
	#
	def select(self,awhat=[""]):

		res = atmList()
		for iAtm in range(0,len(self.atms)):
			if awhat == [""]:
				res.list.append(atmLine(self.atms[iAtm].txt))
			else:
				if awhat[0] !=  "-":
					if awhat.count(self.atms[iAtm].atmName()) > 0:
						res.list.append(atmLine(self.atms[iAtm].txt))
				else:
					if awhat.count(self.atms[iAtm].atmName()) == 0:
						res.list.append(atmLine(self.atms[iAtm].txt))
		return res

	def delete(self,awhat=None):
		"""
		residue.delete([list of atom names])
		This will remove atoms from the residue
		based on their names
		"""
		if awhat == None:
			return
		for iAtm in range(len(self.atms)-1, -1, -1):
			# delete atoms
			if self.atms[iAtm].atmName() in awhat:
				self.atms.list.remove(self.atms.list[iAtm])

	def BBAtmMiss(self, verbose = 0):
		missp = []
		for atms in AABB:
			if self.atms.atmPos(atms) == None:
				missp.append(atms)
				break
		if verbose:
			print missp
		return missp

	def findAtm(self, atmName = "CA", chId = None, rName = None, rNum = None, icode = None, verbose = 0):
		"""
		To identify an atom given residue chain Id, name, PDB number, insertion code,
		and atom Name
		Return:
		either the atom instance
		"""
		if chId != "" and chId != None:
			if self.chnLbl() != chId:
				return None
		if rName != "" and rName != None:
			if self.rName() != rName:
				return None
		if rNum != "" and rNum != None:
			if self.rNum() != rNum:
				return None
		if icode != "" and icode != None:
			if self.riCode() != icode:
				return None
		for aAtm in self.atms:
			if verbose:
				print "aAtm loop", aAtm.atmName()
			if aAtm.atmName() == atmName:
				return aAtm
		if verbose:
			print "findAtm: ", atmName,": None"
		return None


## ========================================
## The PDB file parser
## p.chn("A") does not work !!
## ========================================
class PDB(PDBLine,residue):
	"""
	PDB(PDBLine,residue)
	Classe gerant un PDB.

	x = PDB(fname=\"1tim\", chId = \"\", model = 1, hetSkip = 0, altCare = 0, OXTCare = 0, keepH = 1, id = None, verbose = 0)

	
	"""
	def __init__(self, fname = "", chId = "", model = 1, hetSkip = 0, altCare = 0, OXTCare = 0, PDBMac = 0, keepH = 1, id = None, verbose = 0):

		if fname != "":
			if fname == None:
				return None
			elif isinstance(fname,PDB):                 # already a PDB instance
				if verbose:
					print "PDB from PDBType. hetSkip : ", hetSkip
				self.info  = fname.info
				self.id    = fname.id
				self.mdls  = fname.mdls
				# self.atms  = fname.atms
				iPDB = PDB(fname.flat(altCare = altCare, OXTCare = OXTCare, PDBMac = PDBMac, keepH = keepH), model = model, hetSkip = hetSkip, id = fname.id, verbose = verbose)
				self.atms  = iPDB.atms
				# self.data  = fname.data
				self.data  = iPDB.data
				self.mdls  = iPDB.mdls
				self.seq   = iPDB.seq
				self.seq3D = iPDB.seq3D
				self.ss    = iPDB.ss
				self.s2    = iPDB.s2
				self.nModel = iPDB.nModel
				self.dbref = iPDB.dbref
				self.chns  = iPDB.chns
				self.setModel(model, verbose)
				self.resTab(verbose)
				if altCare :
					self.atms = PDB(self.flat(altCare = 1), verbose = verbose).atms
					self.resTab(verbose)
				if OXTCare :
					# print "altCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)

				if PDBMac :
					# print "altCare set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)

			# a flat series of text lines
			elif isinstance(fname,types.ListType):    # a list of atoms
				if verbose:
					print "PDB from ListType. hetSkip : ", hetSkip
				self.id = "unkwn"
				self = self.parse(fname, "", chId, hetSkip, verbose)
				if id != None:
					self.id = id
				self.setModel(model, verbose)
				self.resTab(verbose)
				if altCare :
					# print "altCare set"
					self.atms = PDB(self.flat(altCare = 1)).atms
					self.resTab(verbose)
				if PDBMac :
					# print "altCare set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)
				if OXTCare :
					# print "altCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)
				if not keepH :
					# print "altCare set"
					self.atms = PDB(self.flat(keepH = 0)).atms
					self.resTab(verbose)
			# from disk file
			elif isinstance(fname,types.StringType):  # read file from disk
				if verbose:
					print "PDB from StringType : "
				self.load(fname, chId, hetSkip, PDBDIR=GDFLTPDBDIR, verbose = verbose)
				# if fname[0] in "0123456789" and len(fname < 6):
				# 	self.id = fname
				if id != None:
					self.id = id
				self.setModel(model, verbose)
				self.resTab(verbose)
				if PDBMac :
					# print "altCare set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)
				if OXTCare :
					# print "altCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)
				if altCare :
					# print "altCare set"
					self.atms = PDB(self.flat(altCare = 1)).atms
					self.resTab(verbose)
				if not keepH :
					# print "altCare set"
					self.atms = PDB(self.flat(keepH = 0)).atms
					self.resTab(verbose)

	# return PDB
	def __getslice__(self, ffrom = 0, tto = None):
		"""
		PDB.__getslice__(ffrom = 0, tto = None)
		PDB[ffrom:tto] (tto excluded, as in python)
		return a PDB instance, a slice of residues
		"""
		# print "__getslice__", ffrom, tto, len(self), tto.__class__
		if tto == None:
			tto = len(self)
		# print "__getslice__", ffrom, tto
		res = self[ffrom].flat()
		for i in range(ffrom+1,tto):
			res= res + self[i].flat()
		return PDB(res)

	# return residue  or PDB (chains)
	def __getitem__(self,aPos):
		"""
		PDB.__getitem__(aPos)
		PDB[aPos]
		PDB["CHAINS"]
		return one residue or PDB instance of chains matching CHAINS (e.g. \"AB\")
		"""
		if isinstance(aPos,types.IntType):
			return self.rt[aPos]
		elif isinstance(aPos,types.StringType):
			res = []
			for i in self:
				if (aPos[0] != "-") and (aPos.count(i.chnLbl())):
					res = res + i.flat()
				elif (aPos[0] == "-") and (not aPos.count(i.chnLbl())):
					res = res + i.flat()
			return PDB(res)
		
	# return number of residue 
	def __len__(self):
		"""
		PDB.__len__()
		number of residues of the PDB instance
		"""
		return len(self.rt)

	# merge two PDB
	def __add__(self,new):
		"""
		PDB.__add__()
		To concatenate 2 PDBs as one
		x = PDB()
		y = PDB()
		z = x + y
		"""

		return PDB(self.flat() + new.flat())
		#return protein(atmList(self.atms[:] + new.atms[:]))

	# del x[i] : ready for deletions !!
	def __delitem__(self,aPos):
		"""
		del x[i]
		preserves anything in comments:
		modify atms
		then performs resTab()

		MIGHT BE BUGGY (P. Tuffery, 2007)
		"""
		try:
			rName  = self[aPos].rName()
			rNum   = self[aPos].rNum()
			riCode = self[aPos].riCode()
			rChn   = self[aPos].chnLbl()
		except:
			return
		
		for iAtm in range(len(self.atms)-1, -1, -1):
			aAtm = self.atms[iAtm]
##			print aAtm
		
			resName = aAtm.resName()
			resNum  = aAtm.resNum()
			iCode   = aAtm.icode()
			chn     = aAtm.chnLbl()
			if (resNum == rNum) and (resName == rName) and (iCode == riCode) and (chn == rChn):
				del self.atms[iAtm]
		self.resTab(0)

	def __repr__(self, altCare = 0, altLbl = "", OXTCare = 0, HSkip = 0):
		"""
		PDB.__repr__
		show atomic information of PDB
		(print PDB ATOM lines)
		"""

		res = ""
## 		for aRes in self:
## 			res = res+aRes.rName()+"_"+aRes.riCode()+"_"+str(aRes.rNum())+"_"+aRes.chnLbl()+"  "
## 		res = res+"\n"
## 		return res
		i = 0
		curChn = ""
		for aRes in self:
			if (i > 0) and (curChn != aRes.chnLbl()):
				res = res + "TER \n"
			curChn = aRes.chnLbl()
			res = res + aRes.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip)
			i += 1
		return res

	# back to a string
	# an internal vital function to transit from PDB
	# to atmList etc
	def flat(self, altCare = 0, altLbl = "", OXTCare = 0, PDBMac = 0, keepH = 1):
		res = []
		for i in self:
			res = res + i.flat(altCare, altLbl, OXTCare = OXTCare, PDBMac = PDBMac, keepH = keepH)
		return res

	#
	# read PDB or PDB chain(s) from disk file
	# chainId: may be a string of several accepted Ids,
	#          or a string starting with - to indicate a list of rejected Ids.
	# hetSkip : 1 to avoid all non peptidic residuesn 0 else
	# model : the number of the model to install (from 1)
	# if header == 0: do not print header line
	# if ter == 0: do not print ter line
	#
	def out(self, outName = "", chainId = "", altCare = 0, altLbl = "", OXTCare = 0, hetSkip = 0, fmode = "w", header = 1, ter = 1, model = 0, end = 0, info = 0, HSkip = 0, allModels = 0, verbose = 0):
		"""
		PDB.out(outName = \"\", chainId = \"\", altCare = 0, altLbl = \"\", OXTCare = 0, hetSkip = 0, fmode = \"w\", header = 1, ter = 1, model = 0, end = 0, info = 0, verbose = 0)
		output PDB content to file
		"""
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName,fmode)
			except:
				sys.stderr.write("Failed to write to %s\n" % outName)
				return

		# print "Opened :",outName,fmode
		# f = sys.stdout
		if header:
			f.write("HEADER                                                        %s\n" % self.id)
		if info:
			for aLine in self.info:
				f.write("%s" % aLine)
		if allModels:
			for aModel in range(1, self.nModel + 1):
				self.setModel(model = aModel)
				f.write("MODEL %d \n" % (aModel))
				res = self.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip )
				f.write("%s" % res)
				f.write("ENDMDL\n")
		else:
			if model:
				f.write("MODEL %d \n" % (model))
			res = self.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip )
			f.write("%s" % res)
			if model:
				f.write("ENDMDL\n")
		if end:
			f.write("END\n")
		f.flush()
		if f != sys.stdout:
			f.close()
		else:
			sys.stdout.flush()

	def xyzout(self, outName = "", chainId = "", hetSkip = 0, fmode = "w", verbose = 0):
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
			
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName,fmode)
			except:
				print "Failed to write to ",outName

		for aCrd in res:
			f.write("%s\n" % aCrd)
		if f != sys.stdout:
			f.close()

	def xyz(self, outName = "", chainId = "", hetSkip = 0, fmode = "w", verbose = 0):
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
			
		return res

	def load(self,fname, chainId = "", hetSkip = 0, PDBDIR = GDFLTPDBDIR, CATHDIR = GDFLTCATHDIR, SCOPDIR = GDFLTSCOPDIR, verbose = 0, model = 1):

		try:
			if verbose:
				print "Trying: ",fname
			allPDB=gsimpleload(fname, 0)
		except IOError, UnboundLocalError:
			if fname[0] in "0123456789" and (len(fname) == 7) and (fname[-1] in "0123456789") and (fname[-2] in "0123456789"):
				# try CATH
				
				try:
					if verbose > 1:
						sys.stderr.write("Trying local CATH entry: %s/%s\n" % (CATHDIR,fname ))
					allPDB=gsimpleload("%s/%s1" % (CATHDIR,fname),0)
				except :
					try: 
						if verbose > 1:
							sys.stderr.write("Trying remote CATH entry: http://data.cathdb.info/v3_2_0/pdb/%s\n" % (fname ))

						from urllib import urlretrieve
						file, log = urlretrieve("http://data.cathdb.info/v3_2_0/pdb/%s" % fname)
						allPDB=simpleload(file,0)
						if verbose > 2:
							print "urlretrieve at %s" % file
							del urlretrieve
						# print len(allPDB)
						if len(allPDB) < 1:
							raise IOError
					except:
						if verbose:
							print 'Sorry: PDB entry ',pdbEntry,'not found'
						raise UnboundLocalError

			elif fname[0] in "0123456789":
				pdbEntry = fname[:4]
				if chainId == "":
					chainId = fname[4:]
				try:
				## Experimental structure
				## Experimental structure
					if verbose > 1:
						sys.stderr.write("Trying local PDB copy: %s/all/pdb/pdb%s.ent.Z\n" % (PDBDIR, pdbEntry))
					try:
						"""
						Old PDB is .Z
						"""
						allPDB=gsimpleload(PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.Z",0)
					except:
						"""
						wwPDB is .gz
						"""
						allPDB=gsimpleload(PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.gz",0)
				## print InfoChain
				except IOError:
					if verbose:
						print "Failed"
				## Model structure
					try:
						if verbose > 1:
							sys.stderr.write("Trying local PDB model: %s/models/current/pdb/%s/pdb%s.ent.Z\n" % (PDBDIR,pdbEntry[1:3],pdbEntry ))
						try:
							allPDB=gsimpleload(PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.Z",0)
						except:
							allPDB=gsimpleload(PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.gz",0)
					except IOError:
						try:
							if verbose > 1:
								sys.stderr.write("Attempting: PDB from wwpdb.org for %s\n" % pdbEntry)
							from urllib import urlretrieve
							# file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(pdbEntry,pdbEntry))
							# We query the RCSB
							# file, log = urlretrieve("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" %(pdbEntry))
							# We query the wwpdb
							file, log = urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz" % pdbEntry)
							allPDB=simpleload(file,0)
							if verbose > 2:
								print "urlretrieve at %s" % file
							del urlretrieve
							# return
						except:
							if verbose:
								print 'Sorry: PDB entry ',pdbEntry,'not found'
							raise UnboundLocalError
			elif fname[0] in "deg":
				# Astral / scop 
				try:
					subdir  = fname[2:4]
					lname = fname
					if fname[0] != "d":
						lname = "d"+fname[1:]
					if verbose > 1:
						sys.stderr.write("Trying Astral/Scop local copy: %s/%s/%s.ent\n" % (SCOPDIR, subdir, lname))
					allPDB=gsimpleload("%s/%s/%s.ent" % (SCOPDIR, subdir, lname),0)
				except:
					try:
						"""
						We try for astral on the net.
						"""
						if verbose > 1:
							sys.stderr.write("Attempting: Astral/Scop entry from astral.berkeley.edu %s\n" % lname)

						from urllib import urlretrieve
						file, log = urlretrieve("http://astral.berkeley.edu/pdbstyle.cgi?id=%s&output=text" % lname)
						allPDB=simpleload(file,0)
						del urlretrieve
						if verbose:
							sys.stderr.write("Astral entry %s retrieved via the net\n" % lname)
					except:
						print 'Sorry: Astral/SCOP entry '+pdbEntry+'not found'
						# raise UnboundLocalError
						return self
			else:
				print "Sorry: %s does not sound as PDB or Astral/SCOP entry!" % fname
				return self
					
			
##		print "PDB_load: ",len(InfoChain)," lines"
		## print InfoChain

	
		# Organize series of lines
		idName = fname
		if string.find(fname,"/") != -1:
			idName = fname[string.rindex(fname,"/")+1:]
		(prefix, tail) = os.path.split(fname)
		(Id, tail) =  os.path.splitext(tail)

		## self.parse(allPDB, idName[:40]+"_"+chainId, chainId, hetSkip, verbose, model)
		self.parse(allPDB, Id[:40], chainId, hetSkip, verbose, model)
		return self
	#
	# Flat line format to PDB format
	#
	def parse(self, allPDB, id="", chainId = "",  hetSkip = 0, verbose = 0, model = 1):

		# print "parse  chainId ","\""+chainId+"\"","hetSkip ",hetSkip
		if id == "":
			id = "unkwn"
		self.info  = []
		self.id    = string.replace(id," ","")
		self.data  = []   # All ATOM DATA, N MODELS
		self.mdls  = []
		self.atms  = []
		self.seq   = []
		self.seq3D = []
		self.ss    = []
		self.s2    = []
		self.nModel = 0
		self.mdls.append(0)
		self.dbref = ""
		self.chns  = ""

		for curLine in allPDB:

			aLine = PDBLine(curLine)

			#print items
			header = aLine.header()
			if header == "ATOM" or header == "HETATM":
				aLine = atmLine(aLine)
				OK = 0
				if chainId == "":
					OK = 1
				elif chainId[0] != '-':
					if string.count(chainId, aLine.chnLbl()):
						OK = 1
				else:
					if string.count(chainId, aLine.chnLbl()) == 0:
						OK = 1
					
				if OK:
					if hetSkip:
						if AA3.count(aLine.resName()) > 0:
							self.data.append(aLine)
						elif hetSkip == 2:
							if SOLV.count(aLine.resName()) == 0:
								self.data.append(aLine)
					else:
						self.data.append(aLine)
			elif header == "TER":
				## self.data.append(curLine)
				pass
			elif header == "HEADER":
				self.info.append(curLine)
				try:
					if curLine[62:66] != "    ":
						self.id = string.split(curLine[62:])[0]
				except:
					pass
			elif header == "COMPND":
				self.info.append(curLine)
			elif header == "SOURCE":
				self.info.append(curLine)
			elif header == "REMARK":
				self.info.append(curLine)
			elif header == "SEQRES":
				self.seq.append(curLine)
			elif header == "HELIX" or header == "SHEET" or header == "TURN":
				## self.s2.append(allPDB[aLine])
				self.s2.append(curLine)
			elif header == "SSBOND":
				## self.ss.append(curLine)
				self.ss.append(curLine)
			elif header == "DBREF":
				## self.dbref = allPDB[aLine]
				self.dbref = aLine
			elif header == "ENDMDL":
				## self.mdls.append(len(self.data))
## 				## self.nModel = self.nModel+1
				self.mdls.append(len(self.data))
				self.nModel = self.nModel+1
			else:
				if header not in ["MASTER", "END"]:
					self.info.append(curLine)
				
				#return self.atms
		if self.nModel == 0:
			self.nModel = self.nModel+1
		self.mdls.append(len(self.data))
		return self


	# tabulate residues
	def resTab(self, verbose):
		"PDB.resTab"

		start   = 1
		self.rt = []
		curResNum = "-1000"
		curResName = "XXX"
		curICode  = ""
		curChn = ""
		atmFrom = 0

		if len(self.atms) == 0:
		       if verbose:
		              sys.stderr.write("Empty PDB instance\n")
		       return
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
##			print aAtm
		
			resName = aAtm.resName()
			resNum  = aAtm.resNum()
			iCode   = aAtm.icode()
			chn     = aAtm.chnLbl()
			if resNum != curResNum or resName != curResName or iCode != curICode or chn != curChn:
				curResNum = resNum
				curResName = resName
				curICode = iCode
				curChn = chn

				if start:
					start = 0
				else:
					self.rt.append(residue(self.atms[atmFrom:iAtm]))
					atmFrom = iAtm
		self.rt.append(residue(self.atms[atmFrom:iAtm+1]))
		if verbose:
			print " Found ",len(self.rt),' residues'


	#
	# How many models in the PDB ?
	#
	def nModels(self):
		"""
		PDB.nmodels()
		return number of models of PDB file (as defined by MODEL / ENDMDL lines)
		"""
		return self.nModel


	#
	# Install current model
	# x = PDB("/home/raid5/PDB/pdb1g25.ent.gz", hetSkip = 1)
	def setModel(self,model = 1,verbose = 0):
		"""
		PDB.setModel(model = 1,verbose = 0)
		set model \# model as working model
		"""
		if model > self.nModels():
			print "Sorry: no model number ",model," (Total of ",self.nModels(),")"
			return
		self.atms = []
		#print "before: self.atoms",self.atms
		if verbose:
			print "Installing model ",model," (atoms ",self.mdls[model-1]," - ",self.mdls[model],")"
		for aLine in range(self.mdls[model-1],self.mdls[model]):
			self.atms.append(atmLine(self.data[aLine]))
		self.curMdl = model
		#print "after: self.atoms",self.atms
		self.resTab(0)
		return
	
	#
	# what chains in the PDB ?
	#
	def chnList(self):
		"""
		PDB.chnList()
		return a string, the concatenated chain identifiers present in the PDB (including blank).
		"""
		curChn = ""
		self.chns = ""
		for aLine in range(0,len(self.atms)):
			if string.count(self.chns,self.atms[aLine][21]) == 0:
				curChn = self.atms[aLine][21]
				self.chns = self.chns + curChn
		return self.chns

	#
	# what chains in the PDB ?
	#
	def nChn(self):
		"""
		PDB.nChn()
		return a number, the number of chain identifiers present in the PDB (including blank).
		"""
		if self.chns == "":
			return len(self.chnList())
		else:
			return len(self.chns)

	#
	# is there such a chain in the PDB file ?
	#
	def hasChn(self, chnId):
		"""
		PDB.hasChn(chnId)
		check if chain of id chainId is present in the PDB instance.
		"""
		if self.chns == "":
			return string.count(self.chnList(),chnId)
		else:
			return string.count(self.chns,chnId)


	#
	# extract particular chain(s) passed in string chainId
	# the default is to return all the chains
	#
	def chn(self,chainId="", hetSkip = 0):
		"""
		PDB.chn(chnId, hetSkip = 0)
		return a PDB instance of chains of the PDB.
		chainId might contain several chain Ids (e.g. AB)
		if chainId starts by \"-\" (minus) returns all but the chains specified.
		"""
		res = []
		for i in self:
			if chainId == "":
				res = res + i.flat()
			elif (chainId[0] != '-') and ( chainId.count(i.chnLbl()) ):
				res = res + i.flat()
			elif (chainId[0] == '-') and ( not chainId.count(i.chnLbl())):
				res = res + i.flat()
		return PDB(res, hetSkip=hetSkip)

	def PDBHNames(self):
		"""
		PDB.PDBHNames()
		Convert H IUPAC names into H PDB names for 20 standard amino-acids
		"""
		for i in self:
			i.PDBHNames()


	#
	# the molecular type of chain(s) in string chainId
	#
	def chnType(self, chainId = "", verbose = 0):
		"""
		PDB.chnType(chainId = \"\", verbose = 0)
		Heuristic detection if the chain type is one of:
		Protein, RNA, DNA, SOLVENT, HETERO
		"""
		if chainId == "":
			chainId = self.chnList()
#			print chainId
		res = []
		unres = []
		for aChain in chainId:
#			print aChain
			theChain = self.chn(aChain)
			nAA  = 0
			nRNA  = 0
			nDNA  = 0
			nHET = 0
			nH2O = 0
			# print "Chain ",aChain," len: ",len(theChain)
			for i in range(0,len(theChain)):
				# resName = atmList(theChain[i]).resName()
				resName = theChain[i].rName()
				#print "\'"+resName+"\'"
				if AA3.count(resName) > 0:
					nAA = nAA +1
				elif RNA3.count(string.split(resName)[0]) > 0:
					nRNA = nRNA +1
				elif DNA3.count(string.split(resName)[0]) > 0:
					nDNA = nDNA +1
				elif SOLV.count(string.split(resName)[0]) > 0:
					nH2O = nH2O +1
				else:
					nHET = nHET + 1
					if verbose:
						if unres.count(resName) == 0:
							unres.append(resName)
							print unres
							print "Unknown residue type (1)",resName

			if verbose:
				print "nAA : ",nAA," nNA : ",nDNA + nRNA," nHET : ",nHET
			nOTHER = nHET + nDNA + nRNA
			if nOTHER < nAA:
				res.append("Protein")
			elif nAA > nDNA + nRNA:
				res.append("Protein")
			# elif nRNA + nDNA > nHET:
			elif nRNA + nDNA > 0:
				if nRNA > 0:
					res.append("RNA")
				else:
					res.append("DNA")
			else:
				if nH2O > nHET:
					res.append("SOLVENT")
				else:
					res.append("HETERO")
		## return res, nAA, nDNA, nRNA, nHET, nH2O
		if len(chainId) == 1:
			return res[0]
		return res

	#
	# return a list of the residue names
	#
	def resTypes(self, what = "all", types = "all", solvent = True):
		"""
		resTypes(self, what, solvent)
		return: a list of residue names (3 characters), combining mask values of what and solvent

		what ("all","aminoacid","nucleicacid")
		presently not considered
		
		types ("all","standard","largestandard","heteros","strictheteros"
		types mask is for all residues but solvent.
		
		"all" : (default) all residue types.
		"none": no residue types (only solvent mask is effective).
		"std" : all standard residue types (standard amino acids, standard nucleotides).
		"lstd": all amino acid types in addition to std.
		"het" : all non standard residues (includes non standard amino-acids)
		"shet": true heteros groups (does not include non standard amino-acids)

		solvent (True/False): consider solvent
		solvent is set to True by default
		"""
		rs = []
		for i in self:
			rName = i.rName()
			if (solvent == False) and (rName in SOLV):
				continue
			if (types == "none") and (rName not in SOLV):
				continue
			if (types == "std") and ((rName not in AA3STRICT) and (rName not in DNA3) and (rName not in RNA3)):
				continue
			if (types == "lstd") and ((rName not in AA3) and (rName not in DNA3) and (rName not in RNA3)):
				continue
			if (types == "het") and ((rName in AA3STRICT) or (rName in DNA3) or (rName in RNA3)):
				continue
			if (types == "shet") and ((rName in AA3) or (rName in DNA3) or (rName in RNA3)):
				continue
			if rName not in rs:
				rs.append(rName)
		return rs
	#
	# return a selection of (sub) residues
	#
	def select(self,rwhat=[""],awhat=[""]):
		res = []
		for i in self:
			if rwhat == [""]:
				res = res + i.select(awhat).flat()
			elif rwhat[0] !=  "-":
				if rwhat.count(i.rName()) > 0:
					res = res + i.select(awhat).flat()
			else:
				if rwhat.count(i.rName()) == 0:
					res = res + i.select(awhat).flat()
		## print res
		if res == []:
			return None
		return PDB(res, id = self.id)

	#
	# return a selection of (sub) residues for a structure
	# the mask (if specified) is a string of length to-from
	# positions corresponding to '-' will be discarded
	#
	def mask(self,ffrom=0,tto=-1,mask=""):
		res = []
		aPos = 0
		if tto == -1:
			tto = len(self)
		if (mask != "") and (len(mask) < tto-ffrom):
			tto = ffrom + len(mask)
			
		for i in range(ffrom,tto):
			if mask == "" or ((mask != "") and (mask[aPos] != '-')):
				res = res + self[i].flat()
			aPos = aPos + 1
		## print res
		if res == []:
			return None
		return PDB(res, id = self.id)

	# Titre du fichier
 	def header(self):
		title=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items[1:]:
					if string.count(aItem,"-") == 2:
						break
					if title != '':
						title = title + " "
					title = title + aItem
		return title

	# Nature du fichier
 	def compound(self):
		title=''
		for Line in self.info:
			if Line[:6]=='COMPND':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title
				

	# Provenance de la molecule
 	def source(self):
		title=''
		for Line in self.info:
			if Line[:6]=='SOURCE':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title

			
	# Auteur
 	def author(self):
		title=''
		for Line in self.info:
			if Line[:6]=='AUTHOR':
				print Line
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title
				
	# KEYWDS lines
	def  keywords(self):
		keylist = ''
		for Line in self.info:
			if Line[:6]=='KEYWDS':
				keylist=keylist+Line[10:-1]

		aPos = 0
		OK = 1
		while string.find(keylist,'\'',aPos) != -1:
			aPos = string.find(keylist,'\'',aPos)
			afunc = keylist[0:aPos]+"\\"+keylist[aPos:]
			keylist = afunc
			aPos = aPos + 1
		return keylist

	# Date de creation du fichier
 	def date(self):
		date=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items:
					if string.count(aItem,"-") == 2:
						date = aItem
				break
		if date != '':
			return date

		# If no creation date, try revision date
		return self.revdate()

	# Revision date (supposes last revision is first REVDAT
 	def revdate(self):
		date=''
		for Line in self.info:
			if Line[:6]=='REVDAT':
				date=string.split(Line[13:22])[0]
				break
	       ## 	print 'creation fichier',date
		return date

	# method by which crds were generated
	def expmethod(self, verbose = 0):
		for Line in self.info:
			if Line[:6]=='EXPDTA':
				if string.find(Line,'X-RAY DIFFRACTION')!=-1:
					return 'X-RAY DIFFRACTION'
				if string.find(Line,'X-RAY POWDER DIFFRACTION')!=-1:
					return 'X-RAY POWDER DIFFRACTION'
				elif string.find(Line,'NMR')!=-1:
					return 'NMR'
				elif string.find(Line,'ELECTRON DIFFRACTION')!=-1:
					return 'ELECTRON DIFFRACTION'

				elif string.find(Line,'FIBER DIFFRACTION')!=-1:
					return 'FIBER DIFFRACTION'

				elif string.find(Line,'FLUORESCENCE TRANSFER')!=-1:
					return 'FLUORESCENCE TRANSFER'

				elif string.find(Line,'NEUTRON DIFFRACTION')!=-1:
					return 'NEUTRON DIFFRACTION'


				elif string.find(Line,'THEORETICAL MODEL')!=-1:
					return 'THEORETICAL MODEL'

				elif string.find(Line,'SYNCHROTRON')!=-1:
					return 'SYNCHROTRON'

				elif string.find(Line,'ELECTRON MICROSCOPY')!=-1:
					return 'ELECTRON MICROSCOPY'

				else:
					return ''
		# Suppose if resolution set: Xray
		if self.resolution() != -1.:
			return 'X-RAY DIFFRACTION'

	# Coordinates resolution
	def resolution(self, verbose = 0):
		"""
		PDB.resolution()
		returns the Resolution of the file, if specified somewhere.
		(data mining in the REMARK lines)
		"""
		resol = -1.
		for Line in self.info:
			if string.find(Line,'REMARK   2 RESOLUTION')!=-1:
				posMax=string.find(Line,'ANGSTROM')-1
				posMin=string.find(Line,'RESOLUTION')+11
				if posMax!=-1:
					try:
						resol=float(Line[posMin:posMax])
					except ValueError:
						pass
		return resol

	# R Value
	def rvalue(self, verbose = 0):
		"""
		PDB.rvalue()
		returns the RValue of the file, if specified somewhere.
		(data mining in the REMARK lines)
		"""
		R_VALUE = "NULL"
		checkRValue = 0

		for Line in self.info:
			if string.find(Line,'REMARK   3') != -1:
				# Case where it is on the next line !!
				if R_VALUE == "NULL" and ((checkRValue == 1) or (checkRValue == 2)):
					if checkRValue == 1:
						if string.find(Line,'.') != -1:
							# print Line
							pos=string.find(Line,'.')-1
							checkRValue == 0
							try:
								R_VALUE=float(Line[pos:pos+5])
								#print R_VALUE
							except ValueError:
								R_VALUE   = "NULL"
					elif checkRValue == 2:
						startPos = string.find(Line,'VALUE')
						if string.find(Line,'.', startPos) != -1:
							pos=string.find(Line,'.', startPos)-1
							toPos = pos+5
							# check for cases such as: 0.20.
							if string.count(Line,'.', pos,toPos) > 1:
								toPos = string.find(Line,'.', pos+2)
								# print Line[pos:pos+5]
							try:
								R_VALUE=float(Line[pos:toPos])
								#print R_VALUE
							except ValueError:
								R_VALUE   = "NULL"
								#print R_VALUE
					checkRValue = 0
					
				# On one line ?
				if R_VALUE == "NULL" and (string.find(Line,' R ') != -1 or string.find(Line,'R VALUE') != -1 or string.find(Line,'R-VALUE') != -1 or string.find(Line,'R-FACTOR') != -1) and string.find(Line,'TEST') == -1 and string.find(Line,'FREE') == -1 and string.find(Line,'ESTIMATE') == -1 and string.find(Line,'BIN') == -1 and  string.find(Line,'ERROR') == -1:
					#print Line
					startPos = string.find(Line,'R VALUE')
					if startPos == -1:
						startPos = string.find(Line,'R-VALUE')
					if startPos == -1:
						startPos = string.find(Line,'R-FACTOR')
					if startPos == -1:
						if string.find(Line,' R '):
							checkRValue = 2
					if verbose:
						print Line[:-1]
						print Line[startPos:-1]
					if string.find(Line,'.', startPos) != -1:
						pos=string.find(Line,'.', startPos)-1
						toPos = pos+5
						# check for cases such as: 0.20.
						if string.count(Line,'.', pos,toPos) > 1:
							toPos = string.find(Line,'.', pos+2)
						#print Line[pos:pos+5]
						#print Line[pos:toPos]
						try:
							R_VALUE=float(Line[pos:toPos])
							# print R_VALUE
						except ValueError:
							if Line[pos] == 'O':
								try:
									R_VALUE=float(Line[pos+1:toPos])
									#print R_VALUE," O error"
								except ValueError:
									R_VALUE   = "NULL"
							else:
								R_VALUE   = "NULL"
							#print R_VALUE

					else:
						#print "checkRValue = 1"
						checkRValue = 1

		return R_VALUE

	def freervalue(self):

		FREE_R_VALUE   = "NULL"
		for Line in self.info:

			if string.find(Line,'FREE R VALUE') != -1 and string.find(Line,'TEST') == -1 and string.find(Line,'ESTIMATE') == -1 and string.find(Line,'BIN') == -1 and  string.find(Line,'ERROR') == -1:
				if string.find(Line,'.') != -1:
					pos=string.find(Line,'.')-1
					try:
						FREE_R_VALUE=float(Line[pos:pos+5])
					except ValueError:
						FREE_R_VALUE   = "NULL"
		return FREE_R_VALUE

	# def seqres(self,chIds='NOTSPECIFIED', verbose = 0):
	def seqresaa3(self, chIds=None, verbose = 0):
		"""
		PDB.seqres()
		returns the sequence of the PDB as specified in the SEQRES lines.
		"""
		
		if chIds == None:
			chIds = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES':
					if string.count(chIds,Line[11]) == 0:
						chIds = chIds + Line[11]

		if len(chIds) > 1:
			rs = []
		for chId in chIds:
			aseqres = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES' and Line[11] == chId:
					aseqres = aseqres + Line[19:70]+' '

			if len(chIds) > 1:
				rs.append(aseqres.split())
			else:
				rs = aseqres.split()
		return rs
		
	def seqres(self,chIds=None, verbose = 0):
		"""
		PDB.seqres()
		returns the sequence of the PDB as specified in the SEQRES lines.
		"""
		
		if chIds == None:
			chIds = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES':
					if string.count(chIds,Line[11]) == 0:
						chIds = chIds + Line[11]

		if len(chIds) > 1:
			rs = []
		for chId in chIds:
			aseqres = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES' and Line[11] == chId:
					aseqres = aseqres + Line[19:70]+' '

			type = self.chnType(chId)

			if type == 'Protein':
				if verbose:
					sys.stderr.write("seqres for %s\n" % aseqres)
				aa1seq=SEQREStoAA1(aseqres, verbose = verbose)
			elif type == "DNA" or type == "RNA":
				curseqres = string.split(aseqres)
				aa1seq = ""
				for i in curseqres:
					aa1seq = aa1seq + i[0]
			else:
				aa1seq = aseqres

			if len(chIds) > 1:
				rs.append(aa1seq)
			else:
				rs = aa1seq
		return rs

	#
	# Does the file contain only CAs ?
	#
	def CAonly(self,verbose=0):
		"""
		PDB.CAonly()
		Does the file contain only CAs ? (Yes/No)
		"""
		res="Yes"
		for aLine in self.data:
			#print aLine[12:15]
			if string.find(aLine[12:15],"CA")==-1:
				res="No"
				if verbose:
					print 'pas uniquement les CA'
				return res
				break
		if verbose:
			print 'uniquement les CA'
		return res
		

	def SCatmMiss(self, verbose = 0):
		"""
		PDB.SCAtmMiss()
		Does the file has side chains with missing atoms ? (Yes/No)
		(For residues having at least one atomic coordinate present)
		"""
		SCatmMiss=""
		status ="No"
		nSCMiss = 0
## 		if NChaine=='_':
## 			theChain=PDB(infochain,hetSkip=1)
## 		else:	
## 			theChain=PDB(infochain,NChaine,hetSkip=1)

		for i in range(0,len(self)):
			resName = self[i].rName()
			#print "\'"+resName+"\'"
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue

			aaTpe = AA3STRICT.index(resName)
			chaine = ""
			for atm in self[i].atms:
				## chaine=chaine+string.split(str(atm))[2]+' '
				chaine=chaine+atm.atmName()+' '
			if verbose:
				print chaine
			missp = 0
			for atms in AASC[aaTpe]:
				if string.find(chaine,atms)==-1:
					missp = 1
					break
			if missp:
				#print res, missp
				status ="Yes"
				nSCMiss = nSCMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				#print resName, resNum
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				SCatmMiss = SCatmMiss+Res
	
		return 	nSCMiss, SCatmMiss		

	def BBatmMiss(self, verbose = 0):
		"""
		PDB.BBAtmMiss()
		Does the file has backbone with missing atoms ? (Yes/Ext/No)
		(For residues having at least one atomic coordinate present)
		"""
		BBatmMiss=""
		status ="No"
		nBBMiss = 0
## 		if NChaine=='_':
## 			theChain=PDB(infochain,hetSkip=1)
## 		else:	
## 			theChain=PDB(infochain,NChaine,hetSkip=1)

		for i in range(0,len(self)):
			resName = self[i].rName()
			#print "\'"+resName+"\'"
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue

			aaTpe = AA3STRICT.index(resName)
			theCheck = self[i].BBAtmMiss()
			missp = 0
			if theCheck != []:
				missp = 1
			if (missp == 1) and (theCheck[0] == "O") and (self[i].atmPos("OXT") != None):
				missp = 0
## 			chaine = ""
## 			for atm in self[i].atms:
## 				## chaine=chaine+string.split(str(atm))[2]+' '
## 				chaine=chaine+atm.atmName()+' '
## 			if verbose:
## 				print chaine
## 			missp = 0
## 			for atms in AABB:
## 				if string.find(chaine,atms)==-1:
## 					missp = 1
## 					break
			if missp:
				#print res, missp
				status ="Yes"
				if i == 0:
					status = "Ext"
				if i == len(self) -1 and status != "Yes":
					status = "Ext"
				if i > 0 and i < len(self) -1:
					status = "Yes"
				nBBMiss = nBBMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				#print resName, resNum
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				BBatmMiss = BBatmMiss+Res
	
		return 	nBBMiss, BBatmMiss		

	def hasAltAtms(self,verbose = 0):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for i in self:
			BB, SC = i.hasAltAtms()
			if BB == "Yes":
				BBAltAtm = "Yes"
			if SC == "Yes":
				SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm
	

	def altAtmsResList(self,verbose = 0):
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		for i in self:
			BB, SC = i.hasAltAtms()
			if BB == "No" and SC == "No":
				continue
			resName = i.rName()
			resNum = i.rNum()
			#print resName, resNum
			icode  = i.riCode()
			lbl  = i.chnLbl()
			if icode == ' ':
				icode = ''
			if lbl == ' ':
				lbl = ''
			resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "

			if BB == "Yes":
				nBBAltAtm = nBBAltAtm + 1
				BBAltAtm = BBAltAtm + resLabel
			if SC == "Yes":
				nSCAltAtm = nSCAltAtm + 1
				SCAltAtm = SCAltAtm + resLabel

		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm


	# 
	# Check if BB peptidic geometry is correct (distance)
	# THIS WILL NOT DETECT FRAGMENTS. IF MANY, THE GAPS ARE IGNORED
	# AND DO NOT RESULT IN "Bad" RETURN.
	# This allows to scan that all the fragments are correct at once.
	#
	def geomCheck(self,verbose=0):

		aN = None
		aC = None
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		
		for aRes in self:
			# aRes = atmList(self[aPos])
			# skip heteros
			if AA3.count(aRes.rName()) == 0:
				continue
			aN = aRes.atmPos("N")
			if aN != None:
				# Nx, Ny, Nz = atmLine.atmCrds(aRes[aN])
				Nx, Ny, Nz = aRes[aN].xyz()
				theN = aRes[aN]
                        if aC != None:
				if theN.chnLbl() == theC.chnLbl():
					aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
				## print Nx, Ny, Nz, Cx, Cy, Cz,aDist,aRes.rName(), aRes.rNum()
					if aDist > 1.50 and aDist < 3.:
						if verbose:
							print "Poor peptidic bond of ",aDist," for ", theC.resName(), theC.resNum(), theN.resName(), theN.resNum()
						if BBGeoOK == "Ok":
							BBGeoOK = "Poor"
					elif aDist > 3.:
						if verbose:
							print "Bad peptidic bond  of ",aDist," for :", theC.resName(), theC.resNum(), theN.resName(), theN.resNum()
						BBGeoOK = "Bad"
			aC  = aRes.atmPos("C")
			if aC != None:
				# Cx, Cy, Cz =atmLine.atmCrds(aRes[aC])
				Cx, Cy, Cz = aRes[aC].xyz()
				theC = aRes[aC]

		return BBGeoOK

	# 
	# Check if BB peptidic geometry is correct (distance)
	# 
	def traceCheck(self,hetSkip = 0, maxCADist = 4.2, verbose = 0):
		theTrace = self.select(awhat=["CA"])
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Ok"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		tracePB = ""

		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1][0].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes - 1]
                                return CisWarning,"No"

			try:
				x2, y2, z2 = theTrace[aRes][0].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes]
                                return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < 3.60: # CIS peptide
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if resName != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + resLabel
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + resLabel

			if aDist > maxCADist: # mauvaise geometrie
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					print "Bad Trace for ",theTrace[aRes-1]," dist = ",aDist," / ",maxCADist

		return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep

	# 
	# Check if BB peptidic geometry is correct (distance)
	# 
	def traceCheck2(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		theTrace = self.select(awhat=["CA"])
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Ok"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		tracePB = ""

		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1][0].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes - 1]
                                return CisWarning,"No"

			try:
				x2, y2, z2 = theTrace[aRes][0].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes]
                                return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < 3.3: # CIS peptide
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if resName != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + resLabel
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + resLabel
				if verbose:
					sys.stderr.write( "%s : CIS Trace for %s dist = %f \n" % (self.id, theTrace[aRes-1],aDist))

			if (aDist < minCADist) and (aDist >= 3.3): # mauvaise geometrie
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					sys.stderr.write( "%s : Bad Trace for %s dist = %f / %f\n" % (self.id, theTrace[aRes-1],aDist,minCADist))

			if aDist > maxCADist: # mauvaise geometrie
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					sys.stderr.write( "%s : Bad Trace for %s dist = %f / %f\n" % (self.id, theTrace[aRes-1],aDist,maxCADist))

		return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep

	def resLabel(self,aRes):
		resName = self[aRes].rName()
		resNum = self[aRes].rNum()
		icode  = self[aRes].riCode()
		lbl  = self[aRes].chnLbl()
		if icode == ' ':
			icode = ''
		if lbl == ' ':
			lbl = ''
		resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
		return resLabel

	def traceCheck3(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		traceOK   = "Ok"
		tracePB   = ""
		nCISPRO   = 0
		nCISPEP   = 0

		for aRes in range(0,len(self)-1):
			CAPos = self[aRes].atms.CApos()
			if CAPos == None:
				traceOK = "No"
				tracePB = tracePB + self.resLabel(aRes)
				continue
			try:
				x1, y1, z1 = self[aRes][CAPos].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", self[aRes].rName(),self[aRes].rNum()
				traceOK = "No"
				tracePB = tracePB + self.resLabel(aRes)
				continue

			CAPos1 = self[aRes+1].atms.CApos()
			if CAPos1 == None:
				traceOK = "No"
				continue
			try:
				x2, y2, z2 = self[aRes+1][CAPos1].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", self[aRes+1].rName(),self[aRes+1].rNum()
				traceOK = "No"
				continue
			
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist > maxCADist: # mauvaise geometrie
				resLabel = self.resLabel(aRes)
				tracePB  = tracePB + resLabel
				if traceOK != "No":
					traceOK = "Bad"
				if verbose:
					sys.stderr.write( "%s : Bad Trace for %s dist = %f / %f\n" % (self.id, resLabel,aDist,maxCADist))
			if aDist < minCADist: # mauvaise geometrie
				# We check the case for CIS
				try: 
					CA  = self[aRes].atms.theAtm("CA")
					C   = self[aRes].atms.theAtm("C")
					x1,y1,z1 = CA.xyz()
					x2,y2,z2 = C.xyz()
				except:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
				try: 
					N   = self[aRes+1].atms.theAtm("N")
					CA2 = self[aRes+1].atms.theAtm("CA")
					x3,y4,z3 = N.xyz()
					x4,y3,z4 = CA2.xyz()
				except:
					resLabel = self.resLabel(aRes+1)
					tracePB  = tracePB + resLabel
				ome = dihedral(x1,y1,z1,x2,y2,z2,x3,y4,z3,x4,y3,z4)
				if abs(ome) < 20:
					resName = self[aRes+1].rName()
					resLabel = self.resLabel(aRes)
					if resName != "PRO":
						nCISPEP += 1
						if verbose:
							sys.stderr.write( "%s : CISPEP conformation for %s ome = %f\n" % (self.id, resLabel, ome))
					else:
						sys.stderr.write( "%s : CISPRO conformation for %s ome = %f\n" % (self.id, resLabel, ome))
						nCISPRO += 1
				else:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
					if traceOK != "No":
						traceOK = "Bad"
					if verbose:
						sys.stderr.write( "%s : Bad Trace for %s dist = %f / %f\n" % (self.id, resLabel,aDist,maxCADist))
					
		return traceOK, tracePB, nCISPEP, nCISPRO

	def CISSeq(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		oSeq = ""
		
		for aRes in range(0,len(self)-1):
			CAPos = self[aRes].atms.CApos()
			if CAPos == None:
				return "PROBLEM"
			try:
				x1, y1, z1 = self[aRes][CAPos].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", self[aRes].rName(),self[aRes].rNum()
				return "PROBLEM"

			CAPos1 = self[aRes+1].atms.CApos()
			if CAPos1 == None:
				return "PROBLEM"
			try:
				x2, y2, z2 = self[aRes+1][CAPos1].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", self[aRes+1].rName(),self[aRes+1].rNum()
				return "PROBLEM"
			
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < minCADist: # mauvaise geometrie
				# We check the case for CIS
				try: 
					CA  = self[aRes].atms.theAtm("CA")
					C   = self[aRes].atms.theAtm("C")
					x1,y1,z1 = CA.xyz()
					x2,y2,z2 = C.xyz()
				except:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
				try: 
					N   = self[aRes+1].atms.theAtm("N")
					CA2 = self[aRes+1].atms.theAtm("CA")
					x3,y4,z3 = N.xyz()
					x4,y3,z4 = CA2.xyz()
				except:
					resLabel = self.resLabel(aRes+1)
					tracePB  = tracePB + resLabel
				ome = dihedral(x1,y1,z1,x2,y2,z2,x3,y4,z3,x4,y3,z4)
				if abs(ome) < 20:
					resName = self[aRes+1].rName()
					resLabel = self.resLabel(aRes)
					if resName != "PRO":
						oSeq += "2"
						if verbose:
							sys.stderr.write( "%s : CISPEP conformation for %s ome = %f\n" % (self.id, resLabel, ome))
					else:
						oSeq += "1"
				else:
					oSeq += "0"
			else:
				oSeq += "0"	
					
		return oSeq


	#
	# determine fragments based on alpha carbon inter-atomic distance alone
	# 4.10 is default threshold
	#
	def chnCAFrgList(self, chId = "", maxDist = 4.10): #x = PDB("12as")

		if chId == "" and len(self.chnList()) > 1:
			print "PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\""
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0

		for aRes in range(1,len(self)):

			# print aRes , "/", len(self)
			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
			# skip non amino acid
				continue

			aC = self[aRes-1].atmPos("CA")
			if aC == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1
			# print Cx,Cy,Cz

			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			
			# print self[aRes-1].rName(),self[aRes-1].rNum()," / ",self[aRes].rName()
			aN = self[aRes].atmPos("CA")
			
			# print aN, self[aRes][aN]
			if aN == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue

			# print aN
			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		# lRes.append(len(self) - 1)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
			
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def asOneChn(self,chnId = ' '):
		for aRes in range(0,len(self)):
			self[aRes].chnLbl(chnId)
			self[aRes].rNum(aRes+1)
		return PDB(self.flat())
			

	#
	# determine fragments based on inter-atomic distance C'-N
	# 1.70 is default threshold
	#
	def chnFrgList(self, chId = "", maxDist = 1.7): #x = PDB("12as")

		if chId == "" and len(self.chnList()) > 1:
			print "PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\""
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0

		for aRes in range(1,len(self)):

			# print aRes , "/", len(self)
			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
			# skip non amino acid
				continue

			aC = self[aRes-1].atmPos("C")
			# print aC, self[aRes-1][aC]
			if aC == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1

			# print Cx,Cy,Cz

			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			
			# print self[aRes-1].rName(),self[aRes-1].rNum()," / ",self[aRes].rName()
			aN = self[aRes].atmPos("N")
			
			# print aN, self[aRes][aN]
			if aN == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue

			# print aN
			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		# lRes.append(len(self) - 1)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
			
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	#
	# This will return the number of fragments
	# and their boundaries
	# This will manage different chains
	#
	def frgList(self, maxNCDist = 1.7, maxCADist = 4.1, verbose = 0): #x = PDB("12as")
		"""
		PDB.frgList(maxNCDist = 1.7, maxCADist = 4.1, verbose = 0)

		return a list of the fragments of the PDB if some geometric inconsistencies are
		detected.
		"""
		# (3, [[0, 326], [327, 655], [656, 860]])
		res = []
		oriRes = 0
		nFrg = 0

		chnIds = self.chnList()

		curDex = 0
		for chId in chnIds:
			curChn = self.chn(chId)

			if self.chnType(chId) != "Protein":
				curDex = curDex + len(curChn)
				continue
			# print len(curChn), len(theBB)

			CAonly = curChn.CAonly()
			if CAonly=="No":
				curNFrg, curFrgList = curChn.chnFrgList(maxDist=maxNCDist)
			else:
				curNFrg, curFrgList = curChn.chnCAFrgList(maxDist=maxCADist)
				# curNFrg = 1
				# curFrgList = [[0,len(curChn)-1]]

			for i in range(0,len(curFrgList)):
				curFrgList[i][0] = curFrgList[i][0] + curDex
				curFrgList[i][1] = curFrgList[i][1] + curDex
				res.append(curFrgList[i])

			# print curFrgList
			nFrg = nFrg + curNFrg
			
			curDex = curDex + len(curChn)
		return nFrg, res

	def nFrg(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)
		return nFrg

	#
	# This will not check for fragments
	#
	def aaseq_ori(self, verbose = 0):
		"""
		PDB.aaseq()

		return the sequence of residues present in the PDB file, having coordinates.
		Converts non standard amino-acids to equivalent standard amino-acid.
		"""
		res = ""
		unres = []
		for aRes in self:
			if AA3STRICT.count(aRes.rName()):
				res = res + AA1[AA3STRICT.index(aRes.rName())]
			elif AA3.count(aRes.rName()):
				rName = aRes.rName()
				if verbose:
					print "Unfrequent residue type: ",
				if rName == "MSE":   # seleno MET
					res = res+"M"
				elif rName == "CSE": # seleno CYS
					res = res+"C"
 				elif rName == "FGL": # amino propane dioique
 					res = res+"S"
				elif rName == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif rName == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif rName == "TRO": # HYDROXY TRYPTOPHANE
					res = res+"W"
				elif rName == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif rName == "MHO": # Hydroxy-MET
					res = res+"M"
				elif rName == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif rName == "HYP": # HYDROXY PRO
					res = res+"P"
				elif rName == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				elif rName == "AYA": # acetyl ALA
					res = res+"A"
				elif rName == "SEG": # hydroxy ALA
					res = res+"A"
				elif rName == "MME": # N methyl MET
					res = res+"M"
## 				elif rName == "BET": # 3methyl GLY
## 					res = res+"G"
## 				elif rName == "DAR": # D ARG
## 					res = res+"R"
				elif rName == "FME": # formyl MET
					res = res+"M"
				elif rName == "CXM": # carboxy MET
					res = res+"M"
				elif rName == "SAC": # acetyl SER
					res = res+"S"
				elif rName == "CSO": # -HYDROXYCYSTEINE
					res = res+"C"
				elif rName == "HTR": # BETA-HYDROXYTRYPTOPHANE
					res = res+"W"
				else:
					res = res+'X'
			else:
				if unres.count(aRes.rName()) == 0:
					unres.append(aRes.rName())
					if verbose:
						print "Unknown residue type (2): ",aRes.rName()
						print unres
				# res = res+'X'
		return res


	#
	# This will not check for fragments
	#
	def aaseq(self, matchAtms = None, verbose = 0):
		"""
		PDB.aaseq()

		return the sequence of residues present in the PDB file, having coordinates.
		Converts non standard amino-acids to equivalent standard amino-acid.
		"""
		res = ""
		unres = []
		for aRes in self:
			if matchAtms != None:
				lOK = 1
				for aAtm in matchAtms:
					if verbose:
						sys.stderr.write("aaseq checking atmName %s" % aAtm)
					if aRes.findAtm(atmName = aAtm, verbose = verbose) == None:
						lOK = 0
						break
				if not lOK:
					continue
			if AA3STRICT.count(aRes.rName()):
				res = res + AA1[AA3STRICT.index(aRes.rName())]
			# elif AA3.count(aRes.rName()):
			elif AA3new.count(aRes.rName()):
				rName = aRes.rName()
				if verbose:
					sys.stderr.write( "%s\n" % "Unfrequent residue type: ")
					
				# au lieu de mettre 10^999999 elif on utilise le dictionnaire
				if (dico_AA.has_key(rName)):
					res = res+dico_AA[rName]
					#print "coucou",rName,"FIN1",dico_AA[rName],"FIN1"
					
				else:
					res = res+'X'
					
			else:
				if unres.count(aRes.rName()) == 0:
					unres.append(aRes.rName())
					if verbose:
						sys.stderr.write( "%s\n" % "Unknown residue type (2): ",aRes.rName())
						sys.stderr.write( "%s\n" %  unres)
			
		return res



	def frgseq(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):

		res = []
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)

		if verbose:
			print "frgSeq:",nFrg, frgList

		for i in frgList:
			res.append( self[i[0]:i[1]+1].aaseq())
		return res

	def SGList(self):
		SGList = []
		for aRes in self:
			if aRes.rName() == "CYS":
				lSGList = []
				for aAtm in aRes.atms:
					if aAtm.atmName() == "SG":
						# print str(aRes[aAtm])
						lSGList.append(aAtm.xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList
	
	def nSSIntra(self):
		nSSBond = 0
		aSGList = self.SGList()
		# print aSGList, len(aSGList)
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				# print aSGList[aRes1][aSG1]
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < 2.35:
							nSSBond = nSSBond + 1
							break


		return nSSBond

	def SSIntra(self):
		SSBonds = []
		for aRes1 in range(0,len(self)):
			if self[aRes1].rName() != "CYS":
				continue
			for aAtm in self[aRes1].atms:
				if aAtm.atmName() == "SG":
					# print str(aRes[aAtm])
					xyz1 = aAtm.xyz()

			for aRes2 in range(aRes1+1,len(self)):
				if self[aRes2].rName() != "CYS":
					continue

				for aAtm in self[aRes2].atms:
					if aAtm.atmName() == "SG":
						# print str(aRes[aAtm])
						xyz2 = aAtm.xyz()
				if apply(distance, xyz1+xyz2) < 2.35:
					SSBonds.append([aRes1, aRes2])
		return SSBonds

		


	def isHalfCys(self, aRes):	
		if self[aRes].rName() != "CYS":
			return 0,0,0
		x = 0.
		y = 0.
		z = 0.
		isSet = 0
		for aAtm in range(0,len(self[aRes].atms)):
			if self[aRes].atms[aAtm].atmName() == "SG":
				x,y,z = self[aRes].atms[aAtm].xyz()
				isSet = 1
		if isSet == 0:
			return 0,0,0
		for aPos in range(0,len(self)):
			if self[aPos].rName() != "CYS":
				continue
			if aPos == aRes:
				continue
			for aAtm in range(0,len(self[aPos].atms)):
				if self[aPos].atms[aAtm].atmName() == "SG":
					x1,y1,z1 = self[aPos].atms[aAtm].xyz()
					if distance(x,y,z,x1,y1,z1) < 2.35:
						return 1, aPos, distance(x,y,z,x1,y1,z1)
		return 0,0,0	
	
	def findRes(self,chId,rName,rNum, icode, what = None, verbose = 0):
		"""
		PDB.findRes(chId,rName,rNum, icode, what = None, verbose = 0)
		
		To identify a residue given its chain Id, name, PDB number, insertion code
		Return:
		either the residue if what == None
		or the residue rank (from 0) if what != None
		"""
		aPos = -1
		for aRes in self:
			aPos = aPos + 1
			if verbose:
				print aRes.chnLbl(),aRes.rName(), aRes.rNum()
			if chId != "" and chId != None:
				if aRes.chnLbl() != chId:
					continue
			if rName != "" and rName != None:
				if aRes.rName() != rName:
					continue
			if rNum != "" and rNum != None:
				if aRes.rNum() != rNum:
					continue
			if icode != "" and icode != None:
				if aRes.riCode() != icode:
					continue
			if what != None:
				return aPos
			return aRes
		return None

	def findAtm(self,chId,rName,rNum, icode, atmName = "CA", verbose = 0):
		"""
		PDB.findAtm(chId,rName,rNum, icode, atmName, verbose = 0)

		To identify an atom given residue chain Id, name, PDB number, insertion code,
		and atom Name
		Return:
		either the atom instance
		"""

		res = self.findRes(chId,rName,rNum, icode)
		if res != None:
			for aAtm in res.atms:
				if aAtm.atmName() == atmName:
					return aAtm
		return None

        ##########################################################################
	#                                                                        #
	# Modification de la structure PDB :                                     #
	#                                                                        #
	#     + suppression des residus PCA                                      #
	#     + transformation des residus MSE en MET                            #
	#          - transformation des atomes SE en S                           #
	#     + transformation des residus CSE en CYS                            #
	#          - transformation des atomes SE en S                           #
	#     + transformation des residus CEA en CYS                            #
	#          - suppression des atomes O1 et HO1                            #
	#     + transformation des residus CGU en GLU                            #
	#          - suppression des atomes CD2, OE3, OE4, HE4                   #
	#     + transformation des residus HTR en TRP                            #
	#          - suppression des atomes O et OH                              #
	#     + transformation des residus TPQ en PHE                            #
	#          - suppression des atomes O2, O4, O5 et HO4                    #
	#                                                                        #
        ##########################################################################
	def clean_ori(self, whatRes = ["PCA", "5HP", "FGL", "MSE", "CSE", "CEA", "CGU", "HTR", "MHO", "IAS", "TPQ", "TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HYP", "TRO"], verbose = 0 ):
		"""
		clean(whatRes = ["5HP", "PCA","MSE", "CSE", "CEA", "CGU", "HTR", "TPQ"])
		cleanup PDB files by converting some non standard residue into standard ones.
		whatRes is the list of residues that may be affected.
		the default is all. whatRes could be only part of the default list.

		+ transformation des residus PCA en GLU
		+ transformation des residus MHO en MET   
		+ transformation des residus IAS en ASP   
		+ transformation des residus HYP en PRO
		+ transformation des residus TPQ en TYR
		+ transformation des residus TRO en TRP
		+ transformation des residus TYS en TYR
		+ transformation des residus MSE en MET   
		   - transformation des atomes SE en S  
		+ transformation des residus CSE en CYS   
		   - transformation des atomes SE en S  
		+ transformation des residus CEA en CYS   
		   - suppression des atomes O1 et HO1   
		+ transformation des residus CGU en GLU   
		   - suppression des atomes CD2, OE3, OE4, HE4   
		+ transformation des residus HTR en TRP            
		   - suppression des atomes O et OH              
		+ transformation des residus TPQ en PHE           
		   - suppression des atomes O2, O4, O5 et HO4
		+ transformation des residus FGL en SER           
		   - suppression des atomes OG1, renomme OG1 en OG
	        + transformation des residus AYA (acetyl ALA) en ALA
		   - suppression des atomes CT, OT, CM
	        + transformation des residus FME (formyl MET) en MET
		   - suppression des atomes OF, CF, (also CN, O1, HCN)
	        + transformation des residus CXM (carboxy MET) en MET
		   - suppression des atomes CN, O1, O2, HO1, HO2
	        +  transformation des residus SAC (acetyl SER) en SER
		   - suppression des atomes C1A, C2A, OAC, 1H2A, 2H2A, 3H2A
	        +  transformation des residus CSO (s-hydroxycysteine) en CYS
		   - suppression des atomes OD, HD
	        +  transformation des residus BET (3methyl GLY) en GLY (NOT BY DEFAULT)
		   - suppression des atomes C1, C2, C3, 1H1, 1H2, 1H3, 2H1, 2H2, 2H3, 3H1, 3H2, 3H3
	        +  transformation des residus DAR (D ARG) en ARG (NOT BY DEFAULT)
	        +  transformation des residus MME (n-methyl MET) en MET
		   - suppression des atomes CM, 1HM, 2HM, 3HM
	        +  transformation des residus SEG (HYDROXYALANINE) en ALA
		   - suppression des atomes OH, HOD
		   

		   AA3 = ["5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS"]
		   AA1seq = "XXXXXWMCXXYMDXECXXY"

		   A faire:
		   
		   FGL->CYS ??, ERROR !
		   TPQ->TYR *done* (not PHE),
		   MHO->MET *done*,
		   IAS->ASP *done: BUT requires atom addition!
		   CGU->GLU *done*,
		   TYS->TYR *done*
		"""
		# for iR in range(len(self)-1,-1,-1):
		# residu = self[iR]
		for residu in self:
			if residu.rName() not in whatRes:
				continue
			if verbose:
				print residu.rName()
			if residu.rName()=="PCA":
				# self.__delitem__(iR)
				residu.rName("GLU")
				residu.delete(["CD","OE"])
			if residu.rName()=="5HP":
				residu.rName("GLU")
				residu.delete(["CD", "OD"])
			if residu.rName()=="FGL":
				residu.rName("SER")
				residu.delete(["OG2"])
				for atom in residu:
					if atom.atmName()=="OG1":
						atom.atmName("OG")
			if residu.rName()=="MSE":
				residu.rName("MET")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CSE":
				residu.rName("CYS")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CEA":
				residu.rName("CYS")
				residu.delete(["O1","HO1"])
			if residu.rName()=="CGU":
				residu.rName("GLU")
				residu.delete(["CD2","OE3","OE4","HE4"])
			if residu.rName()=="HTR":
				residu.rName("TRP")
				residu.delete(["O","HO"])
			if residu.rName()=="MHO":
				residu.rName("MET")
				residu.delete(["OD1"])
			if residu.rName()=="AYA":
				residu.rName("ALA")
				residu.delete(["CT", "OT", "CM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="HYP":
				residu.rName("PRO")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="SEG":
				residu.rName("ALA")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="DAR":
				residu.rName("ARG")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="FME":
				residu.rName("MET")
				residu.delete(["OF", "CF", "CN", "O1", "HCN"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="MME":
				residu.rName("MET")
				residu.delete(["CM", "1HM", "2HM", "3HM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CXM":
				residu.rName("MET")
				residu.delete(["CN", "O1", "O2", "HO1", "HO2"])
				for atom in residu:
					atom.header("ATOM  ")
					atom.header("ATOM  ")
			if residu.rName()=="SAC":
				residu.rName("SER")
				residu.delete(["C1A", "C2A", "OAC", "1H2A", "2H2A", "3H2A"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CSO":
				residu.rName("CYS")
				residu.delete(["OD", "HD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="BET":
				residu.rName("GLY")
				residu.delete(["C1", "C2", "C3", "1H1", "1H2", "1H3", "2H1", "2H2", "2H3", "3H1", "3H2", "3H3"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="TRO":
				residu.rName("TRP")
				residu.delete(["OD1", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="IAS":
				residu.rName("ASP")
			if residu.rName()=="TPQ":
				residu.rName("TYR")
				residu.delete(["O2","O5"])
				for atom in residu:
					if atom.atmName()=="O4":
						atom.atmName("OH")
					if atom.atmName()=="HO4":
						atom.atmName("HH")
			if residu.rName()=="TYS":
				residu.rName("TYR")
				residu.delete(["S","O1","O2","O3","HO3"])
			for atom in residu:
				atom.header("ATOM  ") # 6 chars for header
		return None

	def clean(self, whatRes = AA3new, verbose = 0 ):
		"""
		clean(whatRes = ["5HP", "PCA","MSE", "CSE", "CEA", "CGU", "HTR", "TPQ"])
		cleanup PDB files by converting some non standard residue into standard ones.
		whatRes is the list of residues that may be affected.
		the default is all. whatRes could be only part of the default list.

		+ transformation des residus PCA en GLU
		+ transformation des residus MHO en MET   
		+ transformation des residus IAS en ASP   
		+ transformation des residus HYP en PRO
		+ transformation des residus TPQ en TYR
		+ transformation des residus TRO en TRP
		+ transformation des residus TYS en TYR
		+ transformation des residus MSE en MET   
		   - transformation des atomes SE en S  
		+ transformation des residus CSE en CYS   
		   - transformation des atomes SE en S  
		+ transformation des residus CEA en CYS   
		   - suppression des atomes O1 et HO1   
		+ transformation des residus CGU en GLU   
		   - suppression des atomes CD2, OE3, OE4, HE4   
		+ transformation des residus HTR en TRP            
		   - suppression des atomes O et OH              
		+ transformation des residus TPQ en PHE           
		   - suppression des atomes O2, O4, O5 et HO4
		+ transformation des residus FGL en SER           
		   - suppression des atomes OG1, renomme OG1 en OG
	        + transformation des residus AYA (acetyl ALA) en ALA
		   - suppression des atomes CT, OT, CM
	        + transformation des residus FME (formyl MET) en MET
		   - suppression des atomes OF, CF, (also CN, O1, HCN)
	        + transformation des residus CXM (carboxy MET) en MET
		   - suppression des atomes CN, O1, O2, HO1, HO2
	        +  transformation des residus SAC (acetyl SER) en SER
		   - suppression des atomes C1A, C2A, OAC, 1H2A, 2H2A, 3H2A
	        +  transformation des residus CSO (s-hydroxycysteine) en CYS
		   - suppression des atomes OD, HD
	        +  transformation des residus BET (3methyl GLY) en GLY (NOT BY DEFAULT)
		   - suppression des atomes C1, C2, C3, 1H1, 1H2, 1H3, 2H1, 2H2, 2H3, 3H1, 3H2, 3H3
	        +  transformation des residus DAR (D ARG) en ARG (NOT BY DEFAULT)
	        +  transformation des residus MME (n-methyl MET) en MET
		   - suppression des atomes CM, 1HM, 2HM, 3HM
	        +  transformation des residus SEG (HYDROXYALANINE) en ALA
		   - suppression des atomes OH, HOD
		   

		   AA3 = ["5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS"]
		   AA1seq = "XXXXXWMCXXYMDXECXXY"

		   A faire:
		   
		   FGL->CYS ??, ERROR !
		   TPQ->TYR *done* (not PHE),
		   MHO->MET *done*,
		   IAS->ASP *done: BUT requires atom addition!
		   CGU->GLU *done*,
		   TYS->TYR *done*
		"""
		# for iR in range(len(self)-1,-1,-1):
		# residu = self[iR]
		#print whatRes
		for residu in self:
			if residu.rName() not in whatRes:
				continue
			if verbose:
				print residu.rName()
				
			if residu.rName()=="PCA":
				# self.__delitem__(iR)
				residu.rName("GLU")
				residu.delete(["CD","OE"])
			if residu.rName()=="5HP":
				residu.rName("GLU")
				residu.delete(["CD", "OD"])
			if residu.rName()=="FGL":
				residu.rName("SER")
				residu.delete(["OG2"])
				for atom in residu:
					if atom.atmName()=="OG1":
						atom.atmName("OG")
			if residu.rName()=="MSE":
				residu.rName("MET")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CSE":
				residu.rName("CYS")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CEA":
				residu.rName("CYS")
				residu.delete(["O1","HO1"])
			if residu.rName()=="CGU":
				residu.rName("GLU")
				residu.delete(["CD2","OE3","OE4","HE4"])
			if residu.rName()=="HTR":
				residu.rName("TRP")
				residu.delete(["O","HO"])
			if residu.rName()=="MHO":
				residu.rName("MET")
				residu.delete(["OD1"])
			if residu.rName()=="AYA":
				residu.rName("ALA")
				residu.delete(["CT", "OT", "CM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="HYP":
				residu.rName("PRO")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="SEG":
				residu.rName("ALA")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="DAR":
				residu.rName("ARG")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="FME":
				residu.rName("MET")
				residu.delete(["OF", "CF", "CN", "O1", "HCN"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="MME":
				residu.rName("MET")
				residu.delete(["CM", "1HM", "2HM", "3HM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CXM":
				residu.rName("MET")
				residu.delete(["CN", "O1", "O2", "HO1", "HO2"])
				for atom in residu:
					atom.header("ATOM  ")
					atom.header("ATOM  ")
			if residu.rName()=="SAC":
				residu.rName("SER")
				residu.delete(["C1A", "C2A", "OAC", "1H2A", "2H2A", "3H2A"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CSO":
				residu.rName("CYS")
				residu.delete(["OD", "HD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="BET":
				residu.rName("GLY")
				residu.delete(["C1", "C2", "C3", "1H1", "1H2", "1H3", "2H1", "2H2", "2H3", "3H1", "3H2", "3H3"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="TRO":
				residu.rName("TRP")
				residu.delete(["OD1", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="IAS":
				residu.rName("ASP")
			if residu.rName()=="TPQ":
				residu.rName("TYR")
				residu.delete(["O2","O5"])
				for atom in residu:
					if atom.atmName()=="O4":
						atom.atmName("OH")
					if atom.atmName()=="HO4":
						atom.atmName("HH")
			if residu.rName()=="TYS":
				residu.rName("TYR")
				residu.delete(["S","O1","O2","O3","HO3"])
				
			#ajout des autres residus ici
			if residu.rName()=="OMT":
				residu.rName("MET")
				residu.delete(["OD1","OD2"])
				for atom in residu:
					atom.header("ATOM  ")
					
			if residu.rName()=="ACL":
				residu.rName("ARG")
				residu.delete(["CM","1HM","2HM"])
				for atom in residu:
					if atom.atmName()=="CL":
						atom.atmName("3HM")
					atom.header("ATOM  ")
			
			if residu.rName()=="AGM":
				residu.rName("ARG")
				residu.delete(["1HE2","2HE2","3HE2"])
				for atom in residu:
					if atom.atmName()=="CE2":
						atom.atmName("HD")
					atom.header("ATOM  ")
					
			if residu.rName()=="ARM":
				residu.rName("ARG")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					if atom.atmName()=="CM": #remplace CM par un O
						atom.atmName("O")
					if atom.atmName()=="1HM": #remplace CL par un H
						atom.atmName("H")
					atom.header("ATOM  ")		
					
			if residu.rName()=="HAR":
				residu.rName("ARG")
				for atom in residu:
					if atom.atmName()=="OH1":
						atom.atmName("HH")
					atom.header("ATOM  ")
					
			if residu.rName()=="HMR":
				residu.rName("ARG")
				residu.delete(["2HC","OXT","HXT","2HC"])
				for atom in residu:
					if atom.atmName()=="1HC":
						atom.atmName("O")
					if atom.atmName()=="C":
						atom.atmName("OH")
					if atom.atmName()=="O":
						atom.atmName("HO")
					atom.header("ATOM  ")
				
			if residu.rName()=="AIB":
				residu.rName("ALA")
				residu.delete(["1HB2","2HB2","3HB2"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="CB2":
						atom.atmName("HA")
					
			if residu.rName()=="ALM":
				residu.rName("ALA")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="O":
						atom.atmName("O1")
					if atom.atmName()=="CM":
						atom.atmName("O2")
					if atom.atmName()=="1HM":
						atom.atmName("HO2")
				
			# HETNAM PHE mais suppose ALA --> a verifier 
			if residu.rName()=="BNN":
				residu.rName("ALA")
				residu.delete(["O1","CH","N16","N17","C1","C2","C3","C4","C5","C6",
				"C15","1HH1","2HH1","3HH1","1H16","2H16","H17","H2","H3","H5","H6"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C7":
						atom.atmName("CB")
					if atom.atmName()=="C1":
						atom.atmName("3HB")
					if atom.atmName()=="2H7":
						atom.atmName("1HB")
					if atom.atmName()=="1H7":
						atom.atmName("2HB")
					if atom.atmName()=="C11":
						atom.atmName("2HN")
					if atom.atmName()=="H":
						atom.atmName("1HN")
			
			# EN GLY plutot que ALA
			if residu.rName()=="CHG":
				residu.rName("GLY")
				residu.delete(["H1","C2","C3","C4","C5","C6",
				"1H2","2H2","1H3","2H3","1H4","2H4","1H5","2H5","1H6","2H6"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C7":
						atom.atmName("CA")
					if atom.atmName()=="H7":
						atom.atmName("1HA")
					if atom.atmName()=="N7":
						atom.atmName("N")
					if atom.atmName()=="1HN7":
						atom.atmName("1HN")
					if atom.atmName()=="2HN7":
						atom.atmName("2HN")
					if atom.atmName()=="C1":
						atom.atmName("2HA")
					
			if residu.rName()=="CSD":
				residu.rName("ALA")
				residu.delete(["OD1","OD2","HD2"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="SG":
						atom.atmName("3HB")
					
			if residu.rName()=="DAL":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
				
			# manque 2 H
			if residu.rName()=="DHA":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
				
			if residu.rName()=="DNP":
				residu.rName("ALA")
				residu.delete(["1HG","2HG","3HG"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="NG":
						atom.atmName("3HB")
			
			if residu.rName()=="FLA":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="F1":
						atom.atmName("1HB")
					if atom.atmName()=="F2": 
						atom.atmName("2HB")	
					if atom.atmName()=="F3": 
						atom.atmName("3HB")	
						
			if residu.rName()=="HAC":
				residu.rName("ALA")
				residu.delete(["C2","C3","C4","C5","C6","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","1H6","2H6","H1"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C1": 
						atom.atmName("3HB")
				
			if residu.rName()=="MAA":
				residu.rName("ALA")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="CM": 
						atom.atmName("2HN")
						
			if residu.rName()=="PRR":
				residu.rName("ALA")
				residu.delete(["C4","H4","C3","C9","H9","N1","1H10","2H10","3H10",
				"C10","C2","H2","1H3"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C9": 
						atom.atmName("3H5")
						
			#Manque un H sur OG1
			if residu.rName()=="ALO":
				residu.rName("THR")
				for atom in residu:
					atom.header("ATOM  ")
					
			
			if residu.rName()=="BMT":
				residu.rName("THR")
				residu.delete(["1HH","2HH","3HH","CH","CZ","HZ","HE","CE",
				"1HD2","2HD2","2HD1","1HD1","3HD1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD2": 
						atom.atmName("2HG2")
					if atom.atmName()=="CD1": 
						atom.atmName("3HG2")
					
			if residu.rName()=="DTH":
				residu.rName("THR")
				for atom in residue:
					atom.header("ATOM   ")
					
			if residu.rName()=="TPO":
				residu.rName("THR")
				for atom in residue:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HG1")
				
			if residu.rName()=="BCS":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","CE","CZ1","CZ2","HZ1","HZ2",
				"CT1","CT2","HT1","HT2","CH"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
					
			if residu.rName()=="BUC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3","1H4","2H4","3H4",
				"C4","c3","C2","C1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
						
			if residu.rName()=="C5C":
				residu.rName("CYS")
				residu.delete(["H1","1H2","2H2","1H3","2H3","1H4","2H4",
				"1H5","2H5","C1","C2","C3","C4","C5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
					
			if residu.rName()=="C6C":
				residu.rName("CYS")
				residu.delete(["H1","1H2","2H2","1H3","2H3","1H4","2H4",
				"1H5","2H5","1H6","2H6","C1","C2","C3","C4","C5","C6"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
						
			if residu.rName()=="CCS":
				residu.rName("CYS")
				residu.delete(["HOZ","OZ1","OZ2","CE","1HD","2HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
					
			if residu.rName()=="CME":
				residu.rName("CYS")
				residu.delete(["HO","OH","CZ","1HZ","2HZ","CE","1HE","2HE"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
						
			# Verifier la valence du residu csp
			if residu.rName()=="CSP":
				residu.rName("CYS")
				residu.delete(["O1P","O2P","O3P","PHO2","PHO3","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("H")
						
			if residu.rName()=="CSS":
				residu.rName("CYS")
				residu.delete(["HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
						
			if residu.rName()=="CSW":
				residu.rName("CYS")
				residu.delete(["OD1","OD2"])
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="CSX":
				residu.rName("CYS")
				residu.delete(["OD"])
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="CY3":
				residu.rName("CYS")
				residu.delete(["2H1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1H1":
						atom.atmName("HXT")
					if atom.atmName()=="N1":
						atom.atmName("OXT")
					
			if residu.rName()=="CYG":
				residu.rName("CYS")
				residu.delete(["OE2","CG1","1HG1","2HG1","CB1","1HB1","2HB1",
				"CA1","HA1","N1","1HN1","2HN1","C1","O1","O2","HO2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD1":
						atom.atmName("H")
					
			if residu.rName()=="CYM":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","3HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
						
			if residu.rName()=="DCY":
				residu.rName("CYS")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="EFC":
				residu.rName("CYS")
				residu.delete(["F2","1H2","2H2","C2","C1","1H1","2H1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
					
			if residu.rName()=="OCS":
				residu.rName("CYS")
				residu.delete(["HD2","OD1","OD3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OD2":
						atom.atmName("H")
						
			if residu.rName()=="PEC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","3H5","C1","C2","C3","C4","C5",])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
						
			#Manque un H sur le C-ter
			if residu.rName()=="PEC":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","1HZ","2HZ","1HH","2HH","3HH",
				"CE","CZ","CH"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="H":
						atom.atmName("O")
					if atom.atmName()=="SD":
						atom.atmName("HS")
						
			if residu.rName()=="SCH":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("HS")
						
			#Manque un H sur le C-ter
			if residu.rName()=="SCS":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE","1HZ","2HZ","CZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("HS")
					if atom.atmName()=="H":
						atom.atmName("O")
						
			if residu.rName()=="SCY":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE","OCD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("HS")
							
			if residu.rName()=="SHC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","1H6","2H6","3H6","C6","C2","C3",
				"C4","C5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1":
						atom.atmName("H")
							
			if residu.rName()=="SMC":
				residu.rName("CYS")
				residu.delete(["1HCS","2HCS","3HCS"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CS":
						atom.atmName("H")
						
			if residu.rName()=="SOC":
				residu.rName("CYS")
				residu.delete(["OD2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SE":
						atom.atmName("SG")
					if atom.atmName()=="OD1":
						atom.atmName("H")
						
			if residu.rName()=="ALY":
				residu.rName("LYS")
				residu.delete(["1HH3","2HH3","3HH3","CH3","HO"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("H")
						
			if residu.rName()=="DLY":
				residu.rName("LYS")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="LLP":
				residu.rName("LYS")
				residu.delete(["1H4A","2H4A","C4","C5","C5A","1H5A","2H5A",
				"04P","P","O1P","O2P","O3P","2HOP","3HOP","C6","N1","H6","C2",
				"C2A","1H2A","2H2A","3H2A","C3","O3","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C4A":
						atom.atmName("H")
						
			if residu.rName()=="LLY":
				residu.rName("LYS")
				residu.delete(["C1","O1","O2","HO2","C2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("H")
						
			if residu.rName()=="LYM":
				residu.rName("LYS")
				residu.delete(["3HM","2HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("OXT")
					if atom.atmName()=="1HM":
						atom.atmName("HXT")
						
			if residu.rName()=="LYZ":
				residu.rName("LYS")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OH":
						atom.atmName("2HD")
						
			# Possiblite d'un ASP ?
			if residu.rName()=="SHR":
				residu.rName("LYS")
				residu.delete(["C1","C2","C3","C5","O1","O2","HO1",
				"1H2","2H2","1H3","2H3","O3","O4","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C4":
						atom.atmName("2HN")
						
			if residu.rName()=="TRG":
				residu.rName("LYS")
				residu.delete(["1HH1","2HH1","3HH1","1HH2","2HH2","3HH2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH1":
						atom.atmName("1HZ")
					if atom.atmName()=="CH2":
						atom.atmName("2HZ")
						
			#Valine plutot que Leucine
			if residu.rName()=="BUG":
				residu.rName("VAL")
				residu.delete(["2HN2","1HG3","2HG3","3HG3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="NA":
						atom.atmName("OXT")
					if atom.atmName()=="1HN2":
						atom.atmName("OXT")
					if atom.atmName()=="CG3":
						atom.atmName("H")
						
			if residu.rName()=="CLE":
				residu.rName("LEU")
				residu.delete(["2H2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="N2":
						atom.atmName("OXT")
					if atom.atmName()=="1H2":
						atom.atmName("OXT")
						
			if residu.rName()=="DLE":
				residu.rName("LEU")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="MLE":
				residu.rName("LEU")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
					
			# Transforme en met plutot que leucine
			if residu.rName()=="NLE":
				residu.rName("MET")
				residu.delete(["1HD","2HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
						
			# Transforme en met plutot que leucine
			if residu.rName()=="NLN":
				residu.rName("MET")
				residu.delete(["1HD","2HD","2HH2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
					if atom.atmName()=="NH2":
						atom.atmName("OXT")
					if atom.atmName()=="1HH2":
						atom.atmName("HXT")
					
			# Transforme en met plutot que leucine
			if residu.rName()=="NLP":
				residu.rName("MET")
				residu.delete(["1HD","2HD","O3","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
					if atom.atmName()=="P":
						atom.atmName("C")
						
			if residu.rName()=="CYQ":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","P","O1P","O2P","O3P","2HOP","3HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
						
			if residu.rName()=="DVA":
				residu.rName("VAL")
				for atom in residu:
					atom.header("ATOM   ")

			# Trasnforme en ALA plutot que VAL
			if residu.rName()=="DIV":
				residu.rName("VAL")
				residu.delete(["1HB2","2HB2","3HB2","1HG1","2HG1","3HG1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CB2":
						atom.atmName("HA")
					if atom.atmName()=="CG1":
						atom.atmName("3HB1")
						
			if residu.rName()=="MVA":
				residu.rName("VAL")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
						
			if residu.rName()=="2AS":
				residu.rName("ASP")
				residu.delete(["1HBB","2HBB","3HBB"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CBB":
						atom.atmName("2HB")
						
			
			# Manque un H sur le C-ter
			if residu.rName()=="ASA":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HXT":
						atom.atmName("OXT")
					
			
			if residu.rName()=="ASB":
				residu.rName("ASP")
				residu.delete(["O1","O2","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1":
						atom.atmName("HOD")
						
			if residu.rName()=="ASK":
				residu.rName("ASP")
				residu.delete(["1HM","2HM","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("OXT")
					if atom.atmName()=="3HM":
						atom.atmName("HXT")
						
			if residu.rName()=="ASL":
				residu.rName("ASP")
				residu.delete(["1HC3","2HC3","3HC3","C3","HC2","C1","O1","O2","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C2":
						atom.atmName("HOD")
					
			if residu.rName()=="ASQ":
				residu.rName("ASP")
				residu.delete(["O1P","O2P","O3P","2HOP","3HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOD")
						
			if residu.rName()=="BHD":
				residu.rName("ASP")
				residu.delete(["HOB"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OB":
						atom.atmName("2HB")
						
			if residu.rName()=="DAS":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DSP":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DSN":
				residu.rName("SER")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="MIS":
				residu.rName("SER")
				residu.delete(["C1","C2","C3","1H2","2H2","3H2","H1",
				"1H3","2H3","3H3","1HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOG")
						
			if residu.rName()=="OAS":
				residu.rName("SER")
				residu.delete(["1HC2","2HC2","3HC2","C2A","OAC"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1A":
						atom.atmName("HOG")
						
			# Manque un H sur le C-ter
			if residu.rName()=="SEL":
				residu.rName("SER")
				residu.delete(["2HB2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1HB2":
						atom.atmName("OXT")
						
			if residu.rName()=="SEP":
				residu.rName("SER")
				residu.delete(["2HOP","3HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOG")
						
			if residu.rName()=="SET":
				residu.rName("SER")
				residu.delete(["2HNT"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1HNT":
						atom.atmName("HXT")
					if atom.atmName()=="NT":
						atom.atmName("OXT")
						
			if residu.rName()=="SVA":
				residu.rName("SER")
				residu.delete(["O1","O2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="V":
						atom.atmName("HOG")
					
			if residu.rName()=="DGL":
				residu.rName("GLU")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="GGL":
				residu.rName("GLU")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="GMA":
				residu.rName("GLU")
				residu.delete(["2HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="N2":
						atom.atmName("OXT")
					if atom.atmName()=="1HN":
						atom.atmName("HXT")
					
			if residu.rName()=="DIL":
				residu.rName("ILE")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="IIL":
				residu.rName("ILE")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="GL3":
				residu.rName("GLY")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("OXT")
					if atom.atmName()=="HS":
						atom.atmName("HXT")
					
			# Manque un H sur le C-ter
			if residu.rName()=="GLZ":
				residu.rName("GLY")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HXT":
						atom.atmName("OXT")
					
			if residu.rName()=="GSC":
				residu.rName("GLY")
				residu.delete(["1H1","2H1","1H2","2H2","3H2","C1","C2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("2HA")
						
			if residu.rName()=="MPQ":
				residu.rName("GLY")
				residu.delete(["1HM","2HM","3HM","CD1","CD2","CE1","CE2","CZ",
				"1HD1","1HE1","1HE2","1HD2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("2HN")
					if atom.atmName()=="CG":
						atom.atmName("2HA")
					
			if residu.rName()=="MSA":
				residu.rName("GLY")
				residu.delete(["1HN","3HN","1HG","2HG","3HG"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
					if atom.atmName()=="SB":
						atom.atmName("2HA")
					
			if residu.rName()=="NMC":
				residu.rName("GLY")
				residu.delete(["1HCN","2HCN","CX1","CX2","CX3","1HC2",
				"2HC2","1HC3","2HC3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
					
			if residu.rName()=="SAR":
				residu.rName("GLY")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
					
			if residu.rName()=="DGL":
				residu.rName("GLN")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DHI":
				residu.rName("HIS")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DPN":
				residu.rName("PHE")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DPR":
				residu.rName("PRO")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="DTR":
				residu.rName("TRP")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="DTY":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
						
			if residu.rName()=="TYY":
				residu.rName("TYR")
				residu.delete(["HN5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="O2":
						atom.atmName("HD1")
					if atom.atmName()=="N5":
						atom.atmName("HE2")
						
			if residu.rName()=="TYQ":
				residu.rName("TYR")
				residu.delete(["HN51","HN52","HOZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OZ":
						atom.atmName("HD1")
					if atom.atmName()=="N5":
						atom.atmName("HE2")
						
						
			# Manque un H sur le C-ter
			if residu.rName()=="TYB":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HC":
						atom.atmName("OXT")
						
			if residu.rName()=="STY":
				residu.rName("TYR")
				residu.delete(["O2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("HH")
					
			if residu.rName()=="PTR":
				residu.rName("TYR")
				residu.delete(["O1P","O2P","O3P","PHO2","PHO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HH")
						
			if residu.rName()=="PAQ":
				residu.rName("TYR")
				residu.delete(["HN1","N2","NH2","C1","C2","C3","C4","C5",
				"H2","H3","H4","H5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="O2":
						atom.atmName("HD1")
					if atom.atmName()=="N1":
						atom.atmName("HE2")	
						
			if residu.rName()=="IYR":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="IE":
						atom.atmName("HE")
						
			if residu.rName()=="PHL":
				residu.rName("PHE")
				reside.delete(["H2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="H1":
						atom.atmName("O2")
						
			if residu.rName()=="PHI":
				residu.rName("PHE")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="I":
						atom.atmName("HH")
						
			if residu.rName()=="MEN":
				residu.rName("ASN")
				residu.delete(["1HE2","2HE2","3HE2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CE2":
						atom.atmName("2HD2")
						
			if residu.rName()=="KCX":
				residu.rName("LYS")
				residu.delete(["HX2","OX2","OX1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CX":
						atom.atmName("2H2")
						
			# Manque un H sur l'atome ND1
			if residu.rName()=="3AH":
				residu.rName("HIS")
				residu.delete(["N1","HN1","N2","C3","N4","N3A","1HN3","2HN3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C5":
						atom.atmName("HNE2")
						
			if residu.rName()=="DAH":
				residu.rName("PHE")
				residu.delete(["HOE"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OE2":
						atom.atmName("HE2")
			
			# Manque un H sur l'atome ND1
			if residu.rName()=="HIC":
				residu.rName("HIS")
				residu.delete(["1HZ","2HZ","3HZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CZ":
						atom.atmName("HE2")
						
			if residu.rName()=="HIP":
				residu.rName("HIS")
				residu.delete(["2HOP","3HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HD1")
						
			#Nomenclature atom 31HN ?
			if residu.rName()=="HPQ":
				residu.rName("PHE")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("OXT")
					if atom.atmName()=="1HM":
						atom.atmName("HXT")
						
			if residu.rName()=="LTR":
				residu.rName("TRP")
				for atom in residu:
					atom.header("ATOM   ")
					
			if residu.rName()=="TPL":
				residu.rName("TRP")
				residu.delete(["1HC"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="2HC":
						atom.atmName("O")
					
			if residu.rName()=="MHS":
				residu.rName("HIS")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("HD1")
						
			if residu.rName()=="NEM":
				residu.rName("HIS")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("HE2")
						
			if residu.rName()=="NEP":
				residu.rName("HIS")
				residu.delete(["O1P","O2P","O3P","1HOP","2HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HE2")

			# P. Tuffery, september 2008
# 			if residu.rName()=="CSD":
# 				residu.rName("CYS")
# 				residu.delete(["O1P","O2P","O3P","1HOP","2HOP"])
# 				for atom in residu:
# 					atom.header("ATOM   ")
# 					if atom.atmName()=="P":
# 						atom.atmName("HE2")
			#for atom in residu:
				#atom.header("ATOM  ") # 6 chars for header
		return None

	

	def site(self, verbose = 0):
		"""
		Parse the info lines and check for a site description
		according to the PDB format
		"""
		rs = None
		for Line in self.info:
			if Line[:6]=='SITE  ':
				if rs == None:
					rs = {}
				try:
					siteName = Line[11:14]
				except:
					return rs
				try:
					siteNRes = Line[15:17]
				except:
					return rs
				try:
					resName = Line[18:21]
				except:
					return rs
				try:
					resChn = Line[22]
				except:
					return rs
				try:
					resNum = Line[23:27]
				except:
					return rs
				try:
					resIcode = Line[27]
				except:
					return rs
				if rs.has_key(siteName):
					rs[siteName].append([resName,resChn,resNum,resIcode])
				else:
					rs[siteName] = [[resName,resChn,resNum,resIcode]]
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[29:32]
				except:
					pass
				try:
					resChn = Line[33]
				except:
					pass
				try:
					resNum = Line[34:38]
				except:
					pass
				try:
					resIcode = Line[38]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[40:43]
				except:
					pass
				try:
					resChn = Line[44]
				except:
					pass
				try:
					resNum = Line[45:49]
				except:
					pass
				try:
					resIcode = Line[49]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[51:54]
				except:
					pass
				try:
					resChn = Line[55]
				except:
					pass
				try:
					resNum = Line[56:60]
				except:
					pass
				try:
					resIcode = Line[60]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
		return rs

	def CSAsite(self, id = None, verbose = 0):
		"""
		Attempt to retrieve site from Catalytic Site Atlas:
		at http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=2lzm
		returns either None or a dictionnary of the sites
		Id must be a PDB id. CSA does not consider chains. Hence, one must parse it later.
		Return a list of dictionnaries.
		Status: "Literature" or "PsiBLAST"
		Referer: "PsiBLAST" match
		Comment: Comment on psiBlast and EC.
		Site:    Atoms involved
		"""
		from urllib import urlretrieve
		chIds = None
		if id == None:
			cmd = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=%s" % self.id
		else:
			if len(id) > 4:
				chIds = id[4:]
				id    = id[:4]
				# print "Id: %s ChIds: %s" % (id, chIds)
			cmd = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=%s" % id
		file, log = urlretrieve(cmd)
		if verbose:
			print "urlretrieve at %s" % file
		CSA=simpleload(file,0)
		del urlretrieve

		# print "Id: %s ChIds: %s" % (id, chIds)
		grs = []
		rs = {}
		on = 0
		EC = ""
		Title = ""
		Compound = ""
		StatusOn = False
		TitleCount = 0
		CompoundCount = 0
		for i in CSA:
			# Status du site
			if StatusOn:
				rs["EC"] = EC
				rs["Title"] = Title
				rs["Compound"] = Compound
				if string.count(i,"Literature reference"):
					rs["Status"] = "Literature"
				if string.count(i,"PsiBLAST"):
					# PsiBLAST alignment on <a href="CSA_Site_Wrapper.pl?pdb=1ds2">1ds2</a><br><font color=red>1ds2 has EC code 3.4.21.81 whereas 1ssx has EC code 0....<br>The difference in function suggests that the transfer of annotation from 1ds2 to 1ssx may be incorrect.</font>
					rs["Status"]  = "PsiBLAST"
					s = string.split(i,">")[1]
					r = string.split(s,"<")[0]
					rs["Referer"] = r
					try:
						aPos = string.index(i,">")
						aPos = string.index(i,">", aPos+1)
						aPos = string.index(i,">", aPos+1)
						aPos = string.index(i,">", aPos+1)
						rs["Comment"] = i[aPos+1:]
						rs["Comment"] = string.replace(rs["Comment"],"<br>"," ")
						rs["Comment"] = string.replace(rs["Comment"],"</font>\n","")
					except:
						rs["Comment"] = None
					# print rs["Status"],rs["Referer"], rs["Comment"]
				StatusOn = False
			# New Site
			if string.count(i,"Found by:"):
				# print i
				if rs != {}:
					# print rs
					if rs["Site"] != []:
						grs.append(rs)
				rs = {}
				rs["Site"] = []
				StatusOn = True
				
			if string.count(i,"http://www.ebi.ac.uk/intenz/query?cmd=SearchEC"):
				# <pre><a href="http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&amp;ec=1.1.1.1" target="_top">1.1.1.1</a>
				if EC == "":
					it = string.split(i,">")[2]
					lEC = string.replace(it,"</a","")
					if lEC != "\n":
						EC = lEC

			if TitleCount > 0:
				TitleCount -= 1
				if TitleCount == 0:
					Title = string.replace(i,"</div>","")
					Title = string.replace(Title,"\n","")

			if CompoundCount > 0:
				CompoundCount -= 1
				if CompoundCount == 0:
					Compound = string.replace(i,"</div>","")
					Compound = string.replace(Compound,"\n","")

			if string.count(i,"Title:"):
				TitleCount = 4

			if string.count(i,"Compound:"):
				CompoundCount = 4

			if string.count(i,"Residue"):
				on = 1
				continue
			if on:
				if string.count(i,"<td>"):
					l = string.replace(i,"<td>"," ")
					l = string.replace(l,"</td>\n"," ")
					resName = string.split(l)[0]
					resNum  = string.split(l[6:])[0]
					resChn  = l[5]
					what    = string.split(l[6:])[-1]
					# print resName, resNum, resChn, chIds != None
					if (chIds != None):
						if (resChn in chIds):
							rs["Site"].append( [resName,resNum,resChn,what])
					else:
						rs["Site"].append( [resName,resNum,resChn,what])
				if string.count(i,"</table>"):
					on = 0
		if rs != {}:
			if rs["Site"] != []:
				grs.append(rs)
		return grs

	def exposedAminoAcids(self, rH2O = "1.4", ASALimit = "0.25", what = "E", verbose = 0):
		import ASA
		# y = PDB(self, hetSkip = 2)
		lines = ASA.ASA2(self, rH2O = rH2O, ASALimit = ASALimit, verbose = verbose)
		del ASA
		aaseq  = ""
		be     = ""
		toKeep = []
		for i in range(1,len(lines)-3):
			it = string.split(lines[i])
			# if it[-1] == "E":
			if it[-1] == what: # "E" or "B"
				toKeep.append([it[0],it[1],it[2]])
				be += it[-1]
			else:
				be += "-"
		res = []
		for i in self:
			if [i.rName(),i.chnLbl(),i.rNum()] in toKeep:
				res = res + i.flat()
			elif AA3STRICT.count(i.rName()) == 0:
				res = res + i.flat()
		return PDB(res)

	def BB(self):
		"""
		PDB.BB()
		return a PDB of the backbone only
		"""
		theBB = self.select(awhat=BBATMS)
		return theBB

	def SC(self):
		"""
		PDB.SC()
		return a PDB of the side-chain atoms only
		"""
		theSC = self.select(awhat=SCATMS)
		return theSC

	def around(self, aPos, dist = 3.):
		rs = []
		for aAtm in self[aPos]:
			x, y, z = aAtm.xyz()
			# print aAtm
			# continue
			for aRes in range(0,len(self)):
				if aRes == aPos:
					continue
				aCa = self[aRes].findAtm("CA")
				if aCa == None:
					continue

				CAx, CAy, CAz = aCa.xyz()
				d =  distance(x, y, z, CAx, CAy, CAz)
				if d > 15.:
					continue
				if d < dist:
				#print CAx, CAy, CAz
					if aRes not in rs:
						rs.append(aRes)
					continue
				for aAtm2 in self[aRes]:
					CAx, CAy, CAz = aAtm2.xyz()
					d =  distance(x, y, z, CAx, CAy, CAz)
					if d < dist:
						if aRes not in rs:
							rs.append(aRes)
						break
						
				continue
		rs.sort()
		for aRes in rs:
			print self[aRes].rNum(), self[aRes].rName()
			


def PDBBiologicalUnit(PDBid = None, verbose = 0):
	"""
	get the PDB entry biological unit at PQS server (EBI):
	http://pqs.ebi.ac.uk/pqs-doc/macmol/2lzm.mmol
	"""
	from urllib import urlretrieve
	if PDBid == None:
		return []
	cmd = "http://pqs.ebi.ac.uk/pqs-doc/macmol/%s.mmol" % (string.lower(PDBid))
	file, log = urlretrieve(cmd)
	if verbose:
		print "urlretrieve at %s" % file
	list=simpleload(file,0)
	del urlretrieve
	x  =PDB(list)
	return x
	
def PDBEntries(what = None, isauthor = "no", verbose = 0):
	"""
	To search at pdb.org for a list of entries matching a word
	http://www.pdb.org/pdb/navbarsearch.do?newSearch=yes&isAuthorSearch=no&radioset=All&inputQuickSearch=calpain&image.x=0&image.y=0&image=Search
	"""
	from urllib import urlretrieve
	if what == None:
		return []
	cmd = "http://www.pdb.org/pdb/navbarsearch.do?newSearch=yes&isAuthorSearch=%s&radioset=All&inputQuickSearch=%s&image.x=0&image.y=0&image=Search" % (isauthor, what)
	file, log = urlretrieve(cmd)
	if verbose:
		print "urlretrieve at %s" % file
	list=simpleload(file,0)
	del urlretrieve
## 	f = open("calpain","w")
## 	for i in list:
## 		f.write("%s" % i)
## 	f.close()
	rs = []
	for i in list:
		if string.count(i,"<input type=\"checkbox\" name="):
			j = string.split(string.replace(i,"<input type=\"checkbox\" name=",""))[0]
			k = string.replace(j,"\"","")
			rs.append(k)
	return rs
	
	
	
def CSASite2Escan(PDBid = None, patterns = "strict", purge = 0, verbose = 0):
	"""
	extract Catalytic Site atoms for a PDB.
	gets the information at CSA directly at CSA.
	returns a list of atoms involved.
	Will not take into account psiblast sites.
	patterns: one of "strict", "medium", "light"
	if "strict": exact atom name match
	if "medium": atom name match using atom class compatible pattern
	if "light" : atom name match using light maks (atomic type)
	"""
	"""
	extracted from Catalytic Site Atlas
	"""
	catalyticAtoms = {
		"ASP" : [["CG",  "OD1", "OD2"],["CG",  "OD1", "OD2"]],
		"GLU" : [["CD",  "OE1", "OE2"],["CG",  "CD",  "OE1", "OE2"]],
		"ASN" : [["CG",  "OD1", "ND2"],["CG",  "OD1", "ND2"]],
		"GLN" : [["CD",  "OE1", "NE2"],["CG",  "CD",  "OE1", "NE2"]],
		"HIS" : [["ND1", "NE2"],["CG",  "ND1", "NE2"]],
		"ARG" : [["NH1", "NH2"],["NE",  "NH1", "NH2"]],
		"LYS" : [["CE",  "NZ"],["CD",  "CE",  "NZ"]],
		"SER" : [["CB",  "OG"],["CB",  "OG"]],
		"THR" : [["CB",  "OG", "OG1"], ["CB",  "OG", "OG1"]],
		"CYS" : [["CB",  "SG"],["CB",  "SG"]],
		"ALA" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"GLY" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"LEU" : [["CG",  "CD1", "CD2"],["CG",  "CD1", "CD2"]],
		"PHE" : [["CE1", "CE2", "CZ"],["CE1", "CE2", "CZ"]],
		"TRP" : [["NE1"],["NE1", "CZ2", "CH2"]],
		"TYR" : [["CZ",  "OH"],["CE1", "CZ",  "OH"]],
		"MET" : [["SD"],["SD",  "CE"]],
		}
	"""
	Remark: Some hetero atomes such as ZN/NAD (1qlh) peuvent intervenir
	"""
	"""
	Patterns to match
	The first is for similar class, the second pattern is light
	"""
	matchPatterns = {
		"OD" : ["OD|OE", "O.*"],
		"OE" : ["OD|OE", "O.*"],
		"OG" : ["OG|OH", "SG|OG|OH"],
		"OH" : ["OG|OH", "SG|OG|OH"],
		"NE" : ["ND.*|NE.*", "N.*"],
		"ND2": ["ND2|NE2", "ND|NE"],
		"ND1": ["ND1|NE1", "ND|NE"],
		"NZ" : ["NZ", "N.*"],
		"NH" : ["NH.*", "N.*"],
		"SD" : ["SD", "S.*"],
		"SG" : ["SG", "SG|OG|OH"],
		"CG" : ["CG|CD|CE", "C.*"],
		"CD" : ["CG|CD|CE", "C.*"],
		"CE" : ["CG|CD|CE", "C.*"],
		"CZ" : ["CZ", "C.*"],
		"CH" : ["CH", "C.*"],
		"CB" : ["CB|CG|CD|CE", "C.*"],
		"N"  : ["N","N"],
		"CA" : ["CA","CA"],
		"C"  : ["C","C"],
		}

	rs = None
	if PDBid == None:
		return None
	x = PDB(PDBid)
	if (x == None) or (len(x) == 0):
		return rs
	sites = x.CSAsite()
	if purge:
		sites = purgeSites(sites)
	if verbose == 2:
		print sites
	if len(sites) == 0:
		return None
	chnList = x.chnList()
	if patterns == "light":
		rank = 1
	elif patterns == "medium":
		rank = 0
	grs = []
	if verbose:
		sys.stderr.write( "%s : Found %d site(s)\n" % (PDBid, len(sites)))
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				sys.stderr.write("%s: PsiBlast referer : %s\n" % (PDBid, site["Referer"]))
			continue
		rs = []
		# print site
		cmpLine = "COMPND    %-70s\n" % site["Compound"]
		ECLine  = "REMARK    EC: %-60s\n" % site["EC"]
		rs.append(cmpLine)
		rs.append(ECLine)
		# Check for Hetero Groups
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			if AA3.count(aItem[0]) == 0:
				HTLine = "REMARK    HETGRP: %s\n" % aItem[0]
				rs.append(HTLine)
		# if len(site["Site"]) == 1:
		# Check if monoresidue output !
		count = 0
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			count += 1
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					count -= 1
					break
		if count < 2:
			rs.append("%-70s\n" % "REMARK    MONORESIDUE SITE")
					
		# Here, we format the site
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					if verbose:

						sys.stderr.write("%s: unreferenced catalytic atom %s %s\n" % (PDBid, aAtm.resName(), aAtm.atmName()))
						sys.stderr.flush()
					continue
				if aAtm.atmName() not in catalyticAtoms[aItem[0]][0]:
					continue
				rs.append(aAtm.flat())
				if patterns != "strict":
					try:
						rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
					except:
						try:
							rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
						except:
							pass
				else:
					rs.append("REMARK     RESMATCH. %s\n" % aItem[0])
		if len(rs) > 2:
			grs.append(atmList(rs))
	return grs
	# return rs
	
def purgeSites(sites, hetCheck = 0, verbose = 0):
	rs = []
	for aSite in range(0,len(sites)):
		OK = 1
		for aSite2 in range(aSite+1, len(sites)):
			id = identicalp(sites[aSite],sites[aSite2], hetCheck = hetCheck, verbose =verbose)
			if id == 1:
				 OK = 0
				 break
		if OK:
			rs.append(sites[aSite])
	return rs
				
	
def identicalp(site1, site2, hetCheck = 0, verbose = 0):
	"""
	samep: check is site1 is compatible with site2
	on the basis of residue names, residue number and which part
	(sidechain, backbone)
	hetCheck: if 1, also check heteros
	"""
	if verbose:
		print site1
		print site2
	if site1["EC"] != site2["EC"]:
		return 0
	m = []
	h1 = []
	h2 = []
	for aRes2 in site1["Site"]:
		if AA3.count(aRes2[0]) == 0:
			h1.append(aRes2[0])
	for aRes2 in site2["Site"]:
		m.append(0)
		if AA3.count(aRes2[0]) == 0:
			h2.append(aRes2[0])
	if len(site1["Site"]) - len(h1) != len(site2["Site"]) - len(h2):
		return 0
	for aRes in  site1["Site"]:
		OK = 0
		count = -1
		for aRes2 in site2["Site"]:
			count += 1
			if m[count]: # this residue already assigned
				continue
			if (aRes[0] == aRes2[0]) and (aRes[1] == aRes2[1]) and (aRes[3] == aRes2[3]):
				OK = 1
				# we flag aRes2 as matched
				if verbose:
					print aRes,"matches",aRes2
				break
		if (OK == 0):
			if AA3.count(aRes[0]) == 0:
				if hetCheck == 1:
					return 0
			else:
				if verbose:
					print "No match for",aRes
				return 0
	if verbose:
		print "samep"

	# amino acids are OK
	# what about het groups
	if len(h1) and len(h2):
		if len(h1) == len(h2):
			OK = 1
			for i in h1:
				if i not in h2:
					OK = 0
					break
			return OK
		else:
			return 0
	elif len(h1) or len(h2):
		if hetCheck == 0:
			if len(h1) > len(h2):
				return 2
			else:
				return 1
	return 1

def samep(site1, site2, hetCheck = 0, verbose = 0):
	"""
	samep: check is site1 is compatible with site2
	on the basis of residue names and which part
	(sidechain, backbone)
	hetCheck: if 1, also check heteros
	"""
	if verbose:
		print site1
		print site2
## 	if site1["EC"] != site2["EC"]:
## 		return 0
	h1 = []
	h2 = []
	m = []
	for aRes2 in site1["Site"]:
		if AA3.count(aRes2[0]) == 0:
			h1.append(aRes2[0])
	for aRes2 in site2["Site"]:
		m.append(0)
		if AA3.count(aRes2[0]) == 0:
			h2.append(aRes2[0])
	if len(site1["Site"]) - len(h1) != len(site2["Site"]) - len(h2):
		return 0
	for aRes in  site1["Site"]:
		OK = 0
		count = -1
		for aRes2 in site2["Site"]:
			count += 1
			if m[count]: # this residue already assigned
				continue
			if (aRes[0] == aRes2[0]) and (aRes[3] == aRes2[3]):
				OK = 1
				# we flag aRes2 as matched
				if verbose:
					print aRes,"matches",aRes2
				break
		if (OK == 0):
			if AA3.count(aRes[0]) == 0:
				if hetCheck == 1:
					if verbose:
						print "Het inconsistency"
					return 0
			else:
				if verbose:
					print "No match for",aRes
				return 0
	if verbose:
		print "samep"
	return 1

def EscanCASSites(PDBid = None, patterns = "strict", purge = 0, verbose = 0):
	"""
	extract Catalytic Site atoms for a PDB.
	recurse to validated catalytic sites (follow referer if required).
	Then gets the information at CSA directly at CSA.
	returns a list of atoms involved.
	patterns: one of "strict", "medium", "light"
	if "strict": exact atom name match
	if "medium": atom name match using atom class compatible pattern
	if "light" : atom name match using light maks (atomic type)
	"""
	"""
	extracted from Catalytic Site Atlas
	"""
	catalyticAtoms = {
		"ASP" : [["CG",  "OD1", "OD2"],["CG",  "OD1", "OD2"]],
		"GLU" : [["CD",  "OE1", "OE2"],["CG",  "CD",  "OE1", "OE2"]],
		"ASN" : [["CG",  "OD1", "ND2"],["CG",  "OD1", "ND2"]],
		"GLN" : [["CD",  "OE1", "NE2"],["CG",  "CD",  "OE1", "NE2"]],
		"HIS" : [["ND1", "NE2"],["CG",  "ND1", "NE2"]],
		"ARG" : [["NH1", "NH2"],["NE",  "NH1", "NH2"]],
		"LYS" : [["CE",  "NZ"],["CD",  "CE",  "NZ"]],
		"SER" : [["CB",  "OG"],["CB",  "OG"]],
		"THR" : [["CB",  "OG", "OG1"], ["CB",  "OG", "OG1"]],
		"CYS" : [["CB",  "SG"],["CB",  "SG"]],
		"ALA" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"GLY" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"LEU" : [["CG",  "CD1", "CD2"],["CG",  "CD1", "CD2"]],
		"PHE" : [["CE1", "CE2", "CZ"],["CE1", "CE2", "CZ"]],
		"TRP" : [["NE1"],["NE1", "CZ2", "CH2"]],
		"TYR" : [["CZ",  "OH"],["CE1", "CZ",  "OH"]],
		"MET" : [["SD"],["SD",  "CE"]],
		}
	"""
	Remark: Some hetero atomes such as ZN/NAD (1qlh) peuvent intervenir
	"""
	"""
	Patterns to match
	The first is for similar class, the second pattern is light
	"""
	matchPatterns = {
		"OD" : ["OD|OE", "O.*"],
		"OE" : ["OD|OE", "O.*"],
		"OG" : ["OG|OH", "SG|OG|OH"],
		"OH" : ["OG|OH", "SG|OG|OH"],
		"NE" : ["ND.*|NE.*", "N.*"],
		"ND2": ["ND2|NE2", "ND|NE"],
		"ND1": ["ND1|NE1", "ND|NE"],
		"NZ" : ["NZ", "N.*"],
		"NH" : ["NH.*", "N.*"],
		"SD" : ["SD", "S.*"],
		"SG" : ["SG", "SG|OG|OH"],
		"CG" : ["CG|CD|CE", "C.*"],
		"CD" : ["CG|CD|CE", "C.*"],
		"CE" : ["CG|CD|CE", "C.*"],
		"CZ" : ["CZ", "C.*"],
		"CH" : ["CH", "C.*"],
		"CB" : ["CB|CG|CD|CE", "C.*"],
		"N"  : ["N","N"],
		"CA" : ["CA","CA"],
		"C"  : ["C","C"],
		}

	rs = None
	if PDBid == None:
		return None
	x = PDB(PDBid)
	if (x == None) or (len(x) == 0):
		return rs
	sites = x.CSAsite()
	if purge:
		sites = purgeSites(sites)
	if verbose == 2:
		print sites
	if len(sites) == 0:
		return None
	chnList = x.chnList()
	if patterns == "light":
		rank = 1
	elif patterns == "medium":
		rank = 0
	grs = []
	if verbose:
		sys.stderr.write( "%s : Found %d site(s)\n" % (PDBid, len(sites)))

	# 1. We check that these are validated sites
	tsites = []
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				sys.stderr.write("%s: PsiBlast referer : %s\n" % (PDBid, site["Referer"]))
			y = PDB(site["Referer"])
			ysites = y.CSAsite()
			if purge:
				ysites = purgeSites(ysites)
			for ysite in ysites:
				if ysite["Status"] != "Literature":
 					continue
				# Now we check on residue names and what part (sidechain, etc)
				if samep(site,ysite, verbose = 0):
					ysite["Id"] = site["Referer"]
					tsites.append(ysite)
## 					print "MATCH !!!"
					break
			continue
		else:
			site["Id"] = PDBid
			tsites.append(site)
	sites = tsites
	
	if verbose:
		print "Before purge:"
		print sites
	if purge:
		sites = purgeSites(sites)

	if verbose:
		print "After purge:"
		print sites

	orix = x
	oriChnList = chnList
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				sys.stderr.write("%s: PsiBlast referer : %s\n" % (PDBid, site["Referer"]))
			continue
		if site["Id"] == PDBid:
			x = orix
			chnList = oriChnList
		else:
			# print "Loading %s" % site["Id"]
			x = PDB(site["Id"])
			chnList = x.chnList()
		
		rs = []
		# print site
		if site["Id"] == PDBid:
			# PDB HEADER ID is 63-66
			# headLine = "HEADER   %s CSA literature site %s\n" % (x.id,x.id)
			headLine = "HEADER    %s CSA literature site                            %s\n" % (x.id,x.id)
		else:
			# headLine = "HEADER   %s CSA psiblast site  %s\n" % (x.id,x.id)
			headLine = "HEADER    %s CSA psiblast %s site                         %s\n" % (PDBid, x.id,x.id)
		cmpLine = "COMPND    %-70s\n" % site["Compound"]
		ECLine  = "REMARK    EC: %-60s\n" % site["EC"]
		rs.append(headLine)
		rs.append(cmpLine)
		rs.append(ECLine)
		
		# Check for Hetero Groups
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			if AA3.count(aItem[0]) == 0:
				HTLine = "REMARK    HETGRP: %s\n" % aItem[0]
				rs.append(HTLine)
		# if len(site["Site"]) == 1:
		# Check if monoresidue output !
		count = 0
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			count += 1
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					count -= 1
					break
		if count < 2:
			rs.append("%-70s\n" % "REMARK    MONORESIDUE SITE")
					
		# Here, we format the site
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					if verbose:

						sys.stderr.write("%s: unreferenced catalytic atom %s %s\n" % (PDBid, aAtm.resName(), aAtm.atmName()))
						sys.stderr.flush()
					continue
				if aAtm.atmName() not in catalyticAtoms[aItem[0]][0]:
					continue
				rs.append(aAtm.flat())
				if patterns != "strict":
					try:
						rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
					except:
						try:
							rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
						except:
							pass
				else:
					rs.append("REMARK     RESMATCH. %s\n" % aItem[0])
		if len(rs) > 2:
			rs.append("END    \n")
			grs.append(atmList(rs))
	return grs
	# return rs


## ========================================
## Protein specific tools
## y = protein(x.chn("A"))
## ========================================
	
class protein(PDB):

	def __init__(self, data, chId = "", model = 1, hetSkip = 0, verbose = 0):
		if data != "":
			if isinstance(data,PDB):
				## print "protein init from PDB"
				self.atms = data.atms
				self.rt   = data.rt
				self.nFrg    = 0
				self.resTypes(verbose)
				self.frgs = []
				self.chns  = data.chns

			elif isinstance(data,types.ListType):
				#print "protein init from listType"
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
## 				self.atms = []
## 				for aLine in data:
## 					self.atms.append(aLine)
				## self.atms = data
				self.resTab(verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = ""
			elif isinstance(data,atmList):
				## print "protein init from atmList"
				self.atms = []
				for aLine in data.atms:
					self.atms.append(aLine)
				##self.atms = data
				PDB.resTab(self,verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = ""
			elif isinstance(data,types.StringType):
				## print "protein init from string"
				self.atms  = []
				self.info  = []
				self.seq   = []
				self.seq3D = []
				self.ss    = []
				self.s2    = []
				self.id    = ""
				self.dbref = ""
				self.chns  = ""
				self.nFrg    = 0
				self.frgs = []
				self.nModel = 0
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
				## print self.nModel
				## self.load(data, chId, hetSkip, verbose)
				##self.setModel(model, verbose)
				self.resTab(verbose)
				#PDB.resTab(self,verbose)
				self.resTypes(verbose)
	def resTypes(self, verbose = 0):
		self.tpe = []
		unres = []
		for aRes in range(0,len(self.rt) -1):
			aAtm = self.rt[aRes][0]
			aLine = self.atms[aAtm]
			if AA3.count(atmLine(aLine).resName()) != 0:
				idex = AA3.index(atmLine(aLine).resName())
				self.tpe.append(idex)
			else:
				if unres.count(atmLine(aLine).resName()) == 0:
					print "Unknown residue type (3): ",atmLine(aLine).resName()
					unres.append(atmLine(aLine).resName())
				self.tpe.append(-1)
				

	def frgList(self):
		res = []
		theBB = self.BB()
		oriRes = 0
		nFrg = 0
		## print len(self), len(theBB)
		for aRes in range(1,len(theBB)):
			#print aRes
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def nFrgs(self):
		res = []
		theBB = self.BB()
		#print "PDB5.nFrgs: theBB len ",len(theBB)
		oriRes = 0
		nFrg = 0
		for aRes in range(1,len(theBB)):
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		return nFrg

	def trace(self,fname = "", chId = "", hetSkip = 0, altSel = " "):
		#print "trace : hetSkip = ", hetSkip
		self.resTab(0)
		res = []
		#print len(self.rt)
		for aRes in range(0,len(self.rt) -1):
			for aAtm in range(self.rt[aRes][0], self.rt[aRes+1][0]):
				aLine = self.atms[aAtm]
				# print atmLine(aLine).resName()
				if (hetSkip == 2) and (AA3STRICT.count(atmLine(aLine).resName()) == 0):
					#print "HETPEP : ", atmLine(aLine).resName()
					break
				if hetSkip and (AA3.count(atmLine(aLine).resName()) == 0):
					#print "HET : ", atmLine(aLine).resName()
					break
				#print atmLine(aLine).resName(), "Checking for CA"
				if atmLine(aLine).atmName() == "CA":
					res.append(aLine)
					break


## 		for aLine in self.atms:
## 			if hetSkip and AA3.count(atmLine(aLine).resName()) == 0:
## 				continue
## 			if atmLine(aLine).atmName() == "CA":
## 				res.append(aLine)
## 		print "trace :", res.__class__
		return atmList(res)

## 	def chis(self):
## 		res = []
## 		for aRes in range(0,len(self)):
## 			res.append(self[aRes].chis())
##  		return res


	def outSeq(self, label, hetSkip = 0, verbose = 0):
		#print self
		#print hetSkip
		theTrace = self.trace("","",hetSkip, verbose)
		#print theTrace
		#sys.exit(0)
		seq = protein(theTrace).aaseq()
		print "> ",label,len(seq)
		while len(seq) > 0:
			print seq[:80]
			seq = seq[80:]

	def outRawSeq(self, hetSkip = 0, verbose = 0):
		#print self
		#print hetSkip
		theTrace = self.trace("","",hetSkip, verbose)
		#print theTrace
		#sys.exit(0)
		seq = protein(theTrace).aaseq()
		print seq

	def aaseq(self, verbose = 0):
		res = ""
		unres = []
		for aRes in self:
			if AA3STRICT.count(aRes[0].resName()):
				res = res + AA1[AA3STRICT.index(aRes[0].resName())]
			elif AA3.count(aRes[0].resName()):
				if verbose:
					print "Unfrequent residue type: ",aRes[0].resName()
				if aRes[0].resName() == "MSE": # seleno MET
					res = res+"M"
				elif aRes[0].resName() == "CSE": # seleno CYS
					res = res+"C"
				elif aRes[0].resName() == "FGL": # Formyl GLY
					res = res+"C"
				elif aRes[0].resName() == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif aRes[0].resName() == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif aRes[0].resName() == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif aRes[0].resName() == "MHO": # Hydroxy-MET
					res = res+"M"
				elif aRes[0].resName() == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif aRes[0].resName() == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				else:
					res = res+'X'
			else:
				if unres.count(aRes[0].resName()) == 0:
					unres.append(aRes[0].resName())
					print "Unknown residue type (2): ",aRes[0].resName()
					print unres
		return res
	
	def frg(self,whatFrg, frgs = []):
		if frgs == [] and self.frgs == []:
			self.nFrg, self.frgs = self.frgList()
		return protein(self[self.frgs[whatFrg][0]:self.frgs[whatFrg][1]+1])
	
	def hasAltAtms(self,verbose):
		BBAltAtm = "No"
		SCAltAtm = "No"
		for aLine in self.atms:
			if aLine[16] != ' ':
				isAlt = 1
				if string.count(string.digits,aLine[12]):
					isAlt = 0
				if aLine[12] == ' ' and aLine[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = string.split(aLine[12:15])[0]
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
					BBAltAtm = "Yes"
				else:
					SCAltAtm = "Yes"
		return BBAltAtm, SCAltAtm
	
	def altAtmsResList(self,verbose):
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		# print "altAtmsResList"
		for aPos in range(0,len(self)):
			curRes = self[aPos]
			for aLine in curRes.atms:
				if aLine[16] != ' ':
					isAlt = 1
					if string.count(string.digits,aLine[12]):
						isAlt = 0
					if aLine[12] == ' ' and aLine[13] == 'H':
						isAlt = 0
					if isAlt == 0:
						continue
					theAtmTpe = string.split(aLine[12:15])[0]
					#print aLine
					res    = aLine.resName()
					resNum = aLine.resNum()
					#print i, resNum
					icode  = aLine.icode()
					lbl    = aLine.chnLbl()
					if icode == ' ':
						icode = ''
					if lbl == ' ':
						lbl = ''
					Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
					if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
						nBBAltAtm = nBBAltAtm + 1
						BBAltAtm = BBAltAtm + Res
						break
					else:
						nSCAltAtm = nSCAltAtm + 1
						SCAltAtm = SCAltAtm + Res
						break

		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm

	# Check if all BB atoms are present
	def hasAllBBAtms(self,verbose):
		CAWarning = 0
		CWarning  = 0
		OWarning  = 0
		NWarning  = 0
		residuNameManquant=[]
		cp=0
		for aPos in range(0,len(self)):
		
			#aRes = atmList(self[aPos])
			aRes = self[aPos]
			if aRes.Npos() == None:
				if aPos == 0:
					NWarning  = 1
				elif aPos == len(self) - 1:
					if NWarning < 1:
						NWarning  = 1
				else:
					NWarning  = 2
					cpt=1
					residuNameManquant.append(aPos)
					
			if aRes.CApos() == None:
				if aPos == 0:
					CAWarning  = 1
				elif aPos == len(self) - 1:
					if CAWarning < 1:
						CAWarning  = 1
				else:
					CAWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			if aRes.Cpos() == None:
				if aPos == 0:
					CWarning  = 1
				elif aPos == len(self) - 1:
					if CWarning < 1:
						CWarning  = 1
				else:
					CWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			if aRes.Opos() == None:
				if aPos == 0:
					OWarning  = 1
				elif aPos == len(self) - 1:
					if OWarning < 1:
						OWarning  = 1
				else:
					OWarning  = 2
					if cp==0:
						cp=1
						residuNameManquant.append(aPos)
			
			cp=0

			
		if OWarning == 2 or NWarning == 2 or CAWarning == 2 or CWarning == 2:
			BBAtmMiss = "Yes"
		elif OWarning == 1 or NWarning == 1 or CAWarning == 1 or CWarning == 1:
			BBAtmMiss = "Ext"
		else:
			BBAtmMiss = "No"

		return BBAtmMiss,residuNameManquant


	# Check if BB peptidic geometry is correct (distance)
	def geomCheck(self,verbose):

		aN = None
		aC = None
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		
		for aPos in range(0,len(self)):
			# aRes = atmList(self[aPos])
			aRes = self[aPos]
			aN = aRes.Npos()
			if aN != None:
				# Nx, Ny, Nz = atmLine.atmCrds(aRes[aN])
				Nx, Ny, Nz = aRes[aN].xyz()
                        if aC != None:
                                aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
                                if aDist > 1.50 and aDist < 3.:
                                        if verbose:
                                                print "Poor peptidic bond of ",aDist," for ", resName(theChain[aC]), resNum(theChain[aC]), resName(theChain[aN]), resNum(theChain[aN])
                                        if BBGeoOK == "Ok":
                                                BBGeoOK = "Poor"
                                elif aDist > 3.:
                                        if verbose:
                                                print "Bad peptidic bond  of ",aDist," for :", resName(theChain[aC]), resNum(theChain[aC]), resName(theChain[aN]), resNum(theChain[aN])
                                        BBGeoOK = "Bad"
			aC  = aRes.Cpos()
			if aC != None:
				# Cx, Cy, Cz =atmLine.atmCrds(aRes[aC])
				Cx, Cy, Cz = aRes[aC].xyz()

		return BBGeoOK

	# Check if BB peptidic geometry is correct (distance)
	def traceCheck(self,hetSkip = 0, verbose = 0):
		theTrace = self.trace("","",hetSkip, verbose)
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Yes"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""

		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes - 1]
                                return CisWarning,"No"

			try:
				x2, y2, z2 = theTrace[aRes].xyz()
			except ValueError:
				if verbose:
					print fname," Sorry: fname incorrect ATOM format for:", theTrace[aRes]
                                return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			
			if aDist < 3.60: # CIS peptide
				res    = atmLine(self[aRes].atms[0]).resName()
				resNum = atmLine(self[aRes].atms[0]).resNum()
				#print i, resNum
				icode  = atmLine(self[aRes].atms[0]).icode()
				lbl  = atmLine(self[aRes].atms[0]).chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if theTrace[aRes][17:20] != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + Res
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + Res

			if aDist > 4.10: # mauvaise geometrie
				traceOK = "No"
				if verbose:
					print "Bad Trace for ",theTrace[aRes-1]

		return CisWarning, hasCisPRO, hasCisPEP, traceOK, nCISPRO, CISPRO, nCISPep, CISPep


	def BBAngles(self,aRes = -1000):
		res = []
		if aRes == -1000:
			rFrom = 0
			rTo = len(self)
		else:
			rFrom = aRes
			rTo = aRes+1
		for aPos in range(rFrom,rTo):
			phi = -1000.
			psi = -1000.
			ome = -1000.

			if aPos > 0:
				OK = 1
				aAtm = self[aPos-1].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos-1].theAtm("C").xyz()
					b = self[aPos].theAtm("N").xyz()
					c = self[aPos].theAtm("CA").xyz()
					d = self[aPos].theAtm("C").xyz()
					phi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("N").xyz()
					b = self[aPos].theAtm("CA").xyz()
					c = self[aPos].theAtm("C").xyz()
					d = self[aPos+1].theAtm("N").xyz()
					psi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("CA")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("CA").xyz()
					b = self[aPos].theAtm("C").xyz()
					c = self[aPos+1].theAtm("N").xyz()
					d = self[aPos+1].theAtm("CA").xyz()
					ome = apply(dihedral,a+b+c+d)
			## print phi,psi,ome
			res.append([phi,psi,ome])
 		return res


	def SGList(self):
		SGList = []
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			if aRes[0].resName() == "CYS":
				lSGList = []
				for aAtm in range(0,len(aRes.atms)):
					if atmLine(aRes.atms[aAtm]).atmName() == "SG":
						# print str(aRes[aAtm])
						lSGList.append(atmLine(aRes.atms[aAtm]).xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList
		
	def nSSIntra(self):
		
		nSSBond = 0
		aSGList = self.SGList()
		# print aSGList, len(aSGList)
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				# print aSGList[aRes1][aSG1]
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < 2.35:
							nSSBond = nSSBond + 1
							break


		return nSSBond


	def BB(self):
		# print "PDB5.BB: chn len ",len(self)
		res = []
		for aLine in self.atms:
			theName = aLine.atmName()
			if BBATMS.count(theName) > 0:
				res.append(aLine)
		#print res.__class__
		return  protein(res)
			
	def SC(self):
		res = []
		for aLine in self.atms:
			theName = atmLine(aLine).atmName()
			if BBATMS.count(theName) == 0:
				res.append(aLine)
		return  PDB(res)			


def PDBList(input = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	This is to organize the iterative treatment of PDB instances.
	As input, we take a list (of lines), or a file (list of lines).
	Each line can specify:
	- a local file
	  the local file may contain a multi PDB separated with HEADER / END lines
	- aPDB Id compatible with the PDB class
	- an url
	We return a list of PDB instances. The list length is the number of
	lines of the input
	"""
	#print "PDBList","<br>"
	#sys.stdout.flush()
	# return
	from urllib import urlretrieve

	rs = []
	if input == None:
		return rs

	# print input
	# return	
	# Try as if input already specifies the data itself
	# try reading PDB(s) as local PDB file

	#
	# String can be:
	#    file
	#    PDBid
	#    url
	#
	# file can be:
	#    PDBList
	#    id/url list
	#
	# url can be:
	#    PDBList
	#    id/url list
	#
	if isinstance(input,types.StringType):
		if verbose:
			print "PDBList: Input is string: %s" % input
		try:
			if verbose:
				print "PDBList: Attempting local multiPDB file"
			open(input).close()
			rs = fileInput(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if rs != []:
				try:
					id = rs[0].rt
					return rs
				except:
					pass
		except:
			pass
		# direct PDBid
		if verbose:
			print "PDBList: Not a local multiPDB file"
		try:
			if verbose:
				print "PDBList: Attempting PDBid", input
			x = PDB(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (x != None) and (len(x)):
				return [x]
		except:
			if verbose:
				print "PDBList: Not a PDB"
			try:
				if verbose:
					print "PDBList: Attempting PDBid at pdb.org"
				file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(input[:4],input[:4]))
				x = PDB(file, chId = string.replace(input[4:]," ",""), hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
				if (x != None) and (len(x)):
					x.id = input
					return [x]
			except:
				pass
			try:
				if verbose:
					print "PDBList: Attempting URL to direct multiPDB or PDB"
				file, log = urlretrieve(input)
				if verbose:
					print file
				x = fileInput(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
				if len(x) == 0:
					if verbose:
						print "PDBList: Not a direct multiPDB URL"
				elif x != None:
					return x
			except:
				pass
			
## 			try: # here we try explicit url if pdb.org url did not returned exception
## 				if verbose:
## 					print "Attempting URL"
## 				file, log = urlretrieve(input)
## 				x = PDBList(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
## 				if x != None:
## 					return x
## 			except:
## 				pass

	# a list of atoms
	elif isinstance(input,types.ListType):    
		if verbose:
			print "PDBList: Input is list"
		try:
			rs = parseInput(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if rs != []:
				return rs
		except:
			pass

	if verbose:
		print "PDBList: Input might be complex input"
	# print "PDBList: for not raw data "

	# From here, we consider we have indirections to the data
	# We expect 1 dataset each valid line
	rs = []
	# try file of Ids or URLs
	if isinstance(input,types.StringType):
		try:
			input = open(input).readlines()
		except:
			try:
				if verbose:
					print "PDBList: Attempting URL to list of proteins"
				file, log = urlretrieve(input)
				input = open(file).readlines()
				if verbose:
					print input
			except:
				pass

	# try list of Ids or URLs
	if isinstance(input,types.ListType):    
		for aInput in input:
			x = None
			# remove \n, \r if any
			try:
				aInput = string.replace(aInput,"\n","")
				aInput = string.replace(aInput,"\r","")
				if aInput == "":
					continue
			except:
				pass
			if aInput[0] == "#":
				continue
			# print "\"",aInput,"\""
			if isinstance(aInput,PDB):                 # already a PDB instance
				rs.append(aInput)
				continue
			elif isinstance(aInput,types.StringType):  # read file from disk
				if verbose:
					print "PDBList: Considering: ",aInput
				try:
					if verbose:
						print "PDBList: Trying local file"
					lines = open(aInput).readlines()
					x = parseInput(lines, Id = aInput, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
					# print "Got local file",len(x),len(x[0])
				except:
					try:
						if verbose:
							print "PDBList: Trying PDBId"
						x = PDB(aInput, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
						x.id = aInput
						if (x != None) and len(x):
							x = [x]
						else:
							x = None
					except:
						pass
					# if len(aInput) == 4:
					# if 1:
					if x == None:
						if verbose:
							print "PDBList: Trying PDB entry (at PDB), ChnIds: \"%s\"" % aInput[4:]
						try:
							file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(aInput[:4],aInput[:4]))
							x = PDB(file, chId = string.replace(aInput[4:]," ",""), hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
							if verbose:
								print "PDBList: Completed PDB at pdb.org"
							if (x != None) and len(x):
								x.id = aInput
								x = [x]
							else:
								x = None
								# print "URL: ",file,len(x),len(x[0])
						except:
							pass
					if x == None and (string.lower(aInput[:4]) == "http"):
						try:
							if verbose:
								print "PDBList: Trying URL file"
							file, log = urlretrieve(aInput)
							x = PDBList(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
							if (x == []):
								x = None
						except:
							pass
						
			if x != None:
				rs += x
	del urlretrieve
	return rs


def fileInput(input = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	fileInput: to read a multiPDB file from disk.
	"""
	
	rs = []
	try: 
		inputs = open(input).readlines()
	except:
		return rs
	return parseInput(inputs, Id = input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
	
def parseInput(inputs, Id = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	parseInputs:
	We have a list of PDBlines. We parse them and return a list of PDBs.
	"""

	rs = []
	count = 0
	hcount = 0
	ecount = 0
	for i in inputs:
		if string.count(i[:6],"HEADER"):
			hcount +=1
		if string.count(i[:3],"END"):
			ecount +=1
	if verbose:
		print "parseInput: detected %d HEADER lines" % hcount
	if hcount < 2:
		# print "Got ",hcount,"HEADER"
		try:
			x = PDB(inputs, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (Id != None) and (x.id == "unkwn"):
				x.id = Id
		except:
			return []
		if (x != None) and len(x):
			return [x]
		return []
	else:
		# print "MultiPDB",hcount,len(inputs)
		hcount = 0
		for i in range(0,len(inputs)):
			if string.count(inputs[i][:6],"HEADER"):
				# print "PDB",hcount,i
				if hcount != 0:
					try:
						x = PDB(inputs[hcount-1:i], hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
						if (x != None) and len(x):
							count += 1
							# print "Id:",x.id
							if x.id == "unkwn":
								if (Id != None):
									x.id = "%s_%d" % (Id,count)
								else:
									x.id = "%s_%d" % ("unkwn",count)
							rs.append(x)
					except:
						pass
				hcount = i+1 # +1 to avoid 0 again
		try:
			x = PDB(inputs[hcount-1:i], hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (x != None) and len(x):
				count += 1
				# print "Id:",x.id, Id
				if x.id == "unkwn":
					if (Id != None):
						x.id = "%s_%d" % (Id,count)
					else:
						x.id = "%s_%d" % ("unkwn",count)
				rs.append(x)
		except:
			pass

	return rs

def outPDBList(pdbList, outName = "", initMode =  "w", altCare = 0, altLbl = "", OXTCare = 0, hetSkip = 0, fmode = "w", header = 1, ter = 1, end = 0, info = 0, verbose = 0):
	for i in range(0,len(pdbList)):
		if i == 0:
			pdbList[i].out(outName = outName, fmode=initMode, header = 1, end = 1, ter = 0)
		else:
			pdbList[i].out(outName = outName, fmode="a", header = 1, end = 1, ter = 0)


def PDBSumHeaders(what):
    url = 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/SearchHeaders.pl'
    values = {'string' : 'hydrolase',
              'toolbar' : 'biobar'}
    values['all'] = "TRUE"
    values['string'] = what
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    the_page = response.read()
    # return the_page
    the_lines =  string.split(html2text(the_page),"\n")
    # print len(the_lines)
    # sys.exit(0)
    rs = []
    for i in the_lines:
        # print i
        if string.count(i,"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode="):
            rs.append(string.split(i,"=")[1].encode())
    return rs

def PDBListFromPDBSum(what, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
    rs = PDBSumHeaders(what)
    if verbose:
	    print what
	    print rs
    pdbrs = PDBList(rs, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
    return pdbrs

def subPDB(pdb, seedSeq):
	"""
	Given a PDB instance, return a part corresponding to the seed sequence 
	@author: P. Tuffery
	@param: pdb : a PPDB instance
	@param: seedSeq: the sequence to fetch (a string)
	@return: a PDB instance corresponding to the seeSeq
	         or None if the seedSeq is not found
	"""
	aas = pdb.aaseq()
	if aas.count(seedSeq):
		pos = aas.index(seedSeq)
		return pdb[pos: pos+len(seedSeq)]
	return None

def distance(coordA,coordB):
	return math.sqrt((float(coordA[0])-float(coordB[0]))**2+(float(coordA[1])-float(coordB[1]))**2+(float(coordA[2])-float(coordB[2]))**2)

		
def RMSD(path_struct1,path_struct2):
	n=0
	resn_s1=PDB(path_struct1)
	resn_s2=PDB(path_struct2)
	somme=0
	for i in xrange(len(resn_s1)-1):
		if resn_s1[i][0][0]=="A" and resn_s2[i][0][0]=="A" :
			n+=1
			coordA=[float(resn_s1[i][1][30:38]),float(resn_s1[i][1][38:46]),float(resn_s1[i][1][46:54])]
			coordB=[float(resn_s2[i][1][30:38]),float(resn_s2[i][1][38:46]),float(resn_s2[i][1][46:54])]			
			somme+=distance(coordA,coordB)
	return somme/n
