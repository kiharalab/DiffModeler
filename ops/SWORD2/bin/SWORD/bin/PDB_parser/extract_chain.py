#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# DSIMB
# UNIVERSITE PARIS DIDEROT



import sys
import copy
import os
import re
import argparse
sys.path.append("/home/www/www-tools/sword/bin/PDB_parser/")
import PDB


def extract_chain(filename, outfile, chain_name):
        pdb_name = filename.split("/")[-1].split('.')[0]
	
        try :
                pdb = PDB.PDB(filename)
        except :
                PDB.clean_PDB(filename)
                pdb = PDB.PDB(filename)
	
	fichier = open(outfile,"w")
	fichier.write(str(pdb[chain_name]))
	fichier.close()


def get_args():
        usage = ("extract_chain.py -i pdbfile -o outfile -c chain")
        parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('-i', dest = "infile", type = str,
        help = "Input pdbfile:")
	parser.add_argument('-o', dest = "outfile", type = str,
        help = "Output pdb")
	parser.add_argument('-c', dest = "chain_name", type = str,
        help = "chain name")
	
	args = parser.parse_args()
	
        return (args.infile, args.outfile, args.chain_name)


def main():
        infile, outfile, chain_name = get_args()
	extract_chain(infile, outfile, chain_name)


if __name__=="__main__":
        main()

