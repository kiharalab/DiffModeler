#!/usr/bin/env python
# FileBasics.py

import string
import sys
import os
import copy
import math
import gzip
## import bz2

def simpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = open(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	return lines

def load(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = open(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def gzsimpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = gzip.GzipFile(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	return lines

def gzload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = gzip.GzipFile(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def bz2simpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = bz2.BZ2File(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	return lines

def bz2load(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	lines = bz2.BZ2File(fname,'r').readlines()
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def zsimpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	cmd = "gunzip -c "+fname
	try:
		lines = os.popen(cmd).readlines()
	except:
		raise IOError
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	return lines

def zload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError
	cmd = "gunzip -c "+fname
	try:
		lines = os.popen(cmd).readlines()
	except:
		raise IOError
	if verbose:
		sys.stderr.write("%s : Read %d lines\n" % (fname, len(lines)))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def gsimpleload(fname, verbose = 0):
	if fname[-3:] == ".gz":
		lines = gzsimpleload(fname, 0)
##	elif fname[-2:] == ".bz2":
##		lines = bz2simpleload(fname, 0)
	elif fname[-2:] == ".Z":
		lines = zsimpleload(fname, 0)
	else:
		lines = simpleload(fname, 0)
	return lines

def gload(fname, verbose = 0):
	if fname[-3:] == ".gz":
		lines = gzload(fname, 0)
##	elif fname[-2:] == ".bz2":
##		lines = bz2load(fname, 0)
	elif fname[-2:] == ".Z":
		lines = zload(fname, 0)
	else:
		lines = load(fname, 0)
	return lines

