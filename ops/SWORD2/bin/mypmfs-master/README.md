# MyPMFs

Postic G., Hamelryck T., Chomilier J., Stratmann D.

## Generate statistical potentials from a user-defined list of protein structures


### INSTALL

Type 'make' in the terminal.
This will create executable binaries named 'scoring' and 'training'.




### GET HELP

Run each program without any argument (or with -h option).




### EXAMPLES


#### Case#1
```
$ ./training -l example/list1.txt -d example/dataset/ -o myPotentials
```

This will create a statistical potential for each residue pair represented by the carbons alpha (n=210; *.nrg files).

The 'myPotentials/' output directory will also contain 3 Tab-Separated Values (.tsv) files with some statistics about the training dataset:
- the atomic pairs ranked by their lowest energy peaks (top_energies.tsv);
- the atomic pairs ranked by their frequencies (top_occurrences.tsv);
- the 100 shortest distances (top_distances.tsv).

Note: The same results can be obtained with the following command:
```
$ ./training -L example/list3.txt -d example/dataset/ -o myPotentials
```
Unlike the -l argument, -L does not require using a list of native protein structures (i.e. a list of PDB codes).
This allows using a set of decoys as an input (each having any type of filename).


```
$ ./scoring -i example/dataset/1BKR.pdb -d myPotentials/
```
This will calculate the pseudo-energy of the structure 1BKR by using the previously computed potentials.




#### Case#2
```
$ ./training -l example/list1.txt -d example/dataset/ -o myPotentials -r CB -p -g
```
This will create statistical potentials, with residues represented by their carbons beta (-r CB)
Each potential will be plotted as a SVG file (-p).
This interatomic squared distances used for the calculations are written into *.dat files (-g).

Note: Any previously created 'myPotentials/parameters.log' file will be overwritten.


```
$ ./scoring -i example/dataset/1BKR.pdb -d myPotentials/ -c -p -w -o myResults
```
The pseudo-energy of 1BKR will be calculated with cubic-interpolated potentials (-c).
These interpolated potentials will be plotted as SVG files (-p).
Two TSV files will be written (-w):
- the pseudo-energy and distance for each atomic pair (data.tsv);
- the pseudo-energy for each residue of the protein sequence (energy_[WINDOW_SIZE].tsv).
All these data are written into 'myResults' directory (-o myResults).

Notes:
- the default representation is now CB (carbons beta), as defined in 'myPotentials/parameters.log';
- the WINDOW_SIZE is defined as in https://prosa.services.came.sbg.ac.at/prosa_help.html




#### Case#3
```
$ ./training -l example/list1.txt -d example/dataset/ -o myPotentials -k e -b SJ-dpi -p
```
Same training as case#1 but with Kernel Density Estimations (KDE)
Here, we use an Epanechnikov kernel (-k e), and the kernel bandwidth is selected with the Sheather-Jones direct plug-in (-b SJ-dpi) method.
Each potential will be plotted as a SVG file (-p).


```
$ ./scoring -i example/dataset/1BKR.pdb -d myPotentials/ -q 10A,11A,12A,13A,14A,15A,16A,17A,18A,19A,20A -z -s 2000
```
Only the residues 10A to 20A of 1BKR will be processed (-q).
A Z-score will be computed to evaluate the absolute structural quality (-z); the more negative, the better the model.
This Z-score will be computed on 2000 random sequence decoys (-s 2000).




#### Case#4

After any training:
```
$ ./scoring -l example/list2.txt -d myPotentials/
```
Multiple inputs: a pseudo-energy will be calculated for each of the 25 structures of the 'example/list2.txt'.
The chain name is provided for 2 structures in this list. By default, all chains found will be processed.




#### Case#5
```
$ ./training -l example/list1.txt -d example/dataset/ -o myRefState -r allatom -W
```
This trains the reference state separately (-W) on all atoms (-r allatom).
A 'frequencies.ref' file is created, which can then be used (-R) to train a statistical potential.


```
$ ./training -l example/list1.txt -d example/dataset/ -o myPotentials -R myRefState/frequencies.ref -r backbone
```
Thus, the observed frequencies are trained on backbones, while the reference state is trained on all atoms.




Contact: guillaume.postic@u-paris.fr
