
                           DSSP
            by Wolfgang Kabsch and Chris Sander

Defines secondary structure and solvent exposure of proteins from
atomic coordinates as given in Protein Data Bank format.
_________________________________________________________________

Dear DSSP user!

To compile DSSP, type:

DsspCompileCC (uses standard C compiler, has been tested under SGI IRIX)

or

DsspCompileGCC (uses GNU C, has been tested under Linux and SGI IRIX)


To run DSSP, type:

dsspcmbi PDBSourcefile DSSPDestinationfile

The CMBI version of DSSP contains some additional code to make sure you get a
correct output even if the PDB file is a bit messy. See the DsspCMBI.c source
code for details. It has been tested with all the official PDB files available.

In the (very unlikely ;-) case that you find any bugs (sorry, unexpected features)
in the program, please report them to Elmar.Krieger@cmbi.ru.nl , as long as they
affect the output given for true PDB files obtained from the RCSB. A bug-fixed
source file will then be mailed back to you.

CMBI
Centre for Molecular and Biomolecular Informatics
KU Nijmegen, Netherlands
http://www.cmbi.ru.nl

__________________________________________________________________________________

The executable may be used freely for academic purposes.
Do not redistribute.

Original algorithm in 1982 with minor
changes in 1988. Fast accessibility calculation in 1994 by Michael Scharf. 
Refer to W.Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637.

Copyright by Wolfgang Kabsch and Chris Sander, 1983, 1985, 1988, 1994, 1995. 
Max Planck Institut fuer Medizinische Forschung
and EMBL, D-69012 Heidelberg, Germany.
