
  _____  _    _  _________________  
 /  ___|| |  | ||  _  | ___ \  _  \ 
 \ `--. | |  | || | | | |_/ / | | | 
  `--. \| |/\| || | | |    /| | | | 
 /\__/ /\  /\  /\ \_/ / |\ \| |/ /  
 \____/  \/  \/  \___/\_| \_|___/   
                            v1.0    


============
Description:
============
SWORD (SWift and Optimized Recognition of structural Domains) is an automated
method that identifies protein domains through the hierarchical clustering of
Protein Units (PUs), which are substructures describing the protein architecture
at an intermediate level, between secondary structures and domains. For a given
protein structure, SWORD can provide multiple alternative decompositions into
domains.

Licence: CeCILL v2


====================
System requirements:
====================
- Operating system: Linux or Mac OS X
- Perl 5.10 or higher
- ‘Getopt::Long' Perl module
- ParsePDB.pm Perl module (B. Bulheller)
- Users must download the hydrogen bond estimation algorithm DSSP (W. Kabsch & C. Sander)
  (http://swift.cmbi.ru.nl/gv/dssp/)


=================
How to run SWORD?
=================
Usage:
./SWORD -i structure [-d] [-p I/O_dir] [-c chain] [-m max assignments] [-v]

Options:
-i, --input          input: either a PDB entry (4 characters) or a structure file; defines I/O_dir
-d, --download       required if query PDB file not in I/O_dir (will download from ftp.wwpdb.org)
-p, --path           if --download: input/output directory (default: ./)
-c, --chain          chain name (1 character, facultative); if not given, all chains will be processed
-m, --max            3, 9 or 15: maximum number of alternative assignments (default = 3)
-v, --verbose        RECOMMANDED: all steps printed on STDOUT

Example #1:
./SWORD -i 1A8Y -d -p myResults/

Example #2:
./SWORD -i myModels/model.B99990001.pdb -m 15 -c A -v


=================
Versions history:
=================
v1.0 (10/02/2015) G. Postic
v0.1 (07/04/2011) R. Chebrek




          |__________________________________
(/////////|_________________________________/
          |

	guillaume.postic@univ-paris-diderot.fr

