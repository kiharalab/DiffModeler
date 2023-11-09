
# DiffModeler
<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/DiffModeler-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>  

DiffModeler is a computational tool using a diffusion model to automatically build full protein complex structure from cryo-EM maps at 0-20A resolution.  

Copyright (C) 2023 Xiao Wang, Han Zhu, Genki Terashi, Daisuke Kihara, and Purdue University. 

License: GPL v3. (If you are interested in a different license, for example, for commercial use, please contact us.) 

Contact: Daisuke Kihara (dkihara@purdue.edu)

For technical problems or questions, please reach to Xiao Wang (wang3702@purdue.edu).

## Citation:

Xiao Wang, Han Zhu, Genki Terashi & Daisuke Kihara. Protein Complex Structure Modeling with Diffusion Model and AlphaFold in cryo-EM maps.bioArxiv, 2023.
```
@article{wang2023DiffModeler,   
  title={Protein Complex Structure Modeling with Diffusion Model and AlphaFold in cryo-EM maps},   
  author={Xiao Wang, Han Zhu, Genki Terashi, and Daisuke Kihara},    
  journal={bioArxiv},    
  year={2023}    
}   
```

## Free Online Server: 
### Input map+single-chain structures: https://em.kiharalab.org/algorithm/DiffModeler
### Input map+sequence: https://em.kiharalab.org/algorithm/DiffModeler(seq)

## Introduction

<details>
   <summary>DiffModeler is a computational tool using a diffusion model to automatically build full protein complex structure from cryo-EM maps at 0-20A resolution. </summary>

Cryogenic electron microscopy (cryo-EM) has been widely employed in experimental settings to determine 
multi-chain protein complexes, but modeling accuracy greatly diminishes when resolution 
decreases. At intermediate resolutions of 5-10 Å, even template-based structure fitting presents
significant challenges. To tackle this issue, we introduce DiffModeler, a fully automated protein complex structure modeling
method that leverages a diffusion model for backbone tracing and structure fitting with AlphaFold predicted single-chain structure.
In extensive testing on cryo-EM maps at intermediate resolution, DiffModeler showcased remarkably accurate
structure modeling, surpassing existing methods significantly. 
Notably, we successfully modeled a protein complex consisting of 47 chains,
comprising 13,462 residues, with an impressive TM-Score of 0.9. 
We also further benchmarked DiffModeler for maps at low resolution of 10-20 Å and 
validated its generalizability with plausible performances. 
</details>

## Overall Protocol 

<details>

1) Backbone tracing from cryo-EM maps at intermediate resolution via diffusion model.  
2) Single-chain structure prediction by AlphaFold.  
3) Single-chain structure fitting using VESPER.  
4) Protein complex modeling by assembling algorithms.  

<p align="center">
  <img src="framework.png" alt="DiffModeler framework" width="70%">
</p>
</details>

## Installation

<details>



### System Requirements
CPU: >=4 cores <br>
Memory (RAM): >=12Gb. <br>
GPU: any GPU supports CUDA with at least 12GB memory. <br>
GPU is required for DiffModeler since most computations are done on GPU.

## Installation  
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone git@github.com:kiharalab/DiffModeler.git && cd DiffModeler
```

### 3. Configure environment for DiffModeler.
#### 3.1.1 Install anaconda
Install anaconda from https://www.anaconda.com/download#downloads.
#### 3.1.2 Install environment via yml file
Then create the environment via
```commandline
conda env create -f environment.yml
```
#### 3.1.3 Activate environment for running
Each time when you want to run this software, simply activate the environment by
```
conda activate DiffModeler
conda deactivate(If you want to exit) 
```

### 4. Download the pre-trained diffusion model
make a directory ``best_model`` and then download our pretrained model in this directory. <br>
Diffusion model weights (trained on 5-10A, can be used for 2-5A (very good backbone tracing) and 10-20A): [diffusion_model](https://huggingface.co/zhtronics/DiffModelerWeight/resolve/main/diffusion_best.pth.tar) <br>

You can also use command line to do this
```commandline
mkdir best_model
cd best_model
wget https://huggingface.co/zhtronics/DiffModelerWeight/resolve/main/diffusion_best.pth.tar
cd ..
```

If the link failed, you can also download our model files via our [lab server](https://kiharalab.org/emsuites/diffmodeler_model/) to ``best_model`` directory. 

### 5. (Optional) Visualization software
Pymol (for structure visualization): https://pymol.org/2/    
Chimera (for map visualization): https://www.cgl.ucsf.edu/chimera/download.html  

</details>

# Usage

<details>
<summary>Command Parameters</summary>

```commandline
usage: main.py [-h] --mode MODE [-F F] [-M M] [--config CONFIG] [--gpu GPU] [--output OUTPUT]
               [--contour CONTOUR]

options:
  -h, --help            show this help message and exit
  --mode MODE           control mode, mode 0: template mode; mode 1: sequence mode (online search); mode 2: sequence mode (local db search)
  -F F                  input map path
  -P P                  directory or zipped file of Single-Chain PDB files
  -M M                  txt file path which records protein information
  --resolution RESOLUTION
                        specify the resolution to skip diffusion for super high resolution maps (better than 2A)
  --config CONFIG       specifying the config path
  --gpu GPU             specify the gpu we will use
  --output OUTPUT       Output directory
  --contour CONTOUR     Contour level for input map, suggested 0.5*[author_contour]. (Float), Default value: 0.0
```
</details>

<details>
<summary>Protein Structure Complex Modeling with provided template</summary>

### Protein Structure Complex Modeling with provided template 
This is for DiffModeler running if you have map and available template candiates (either from experimental or AlphaFold predicted structure).
It is fine to run if you only know some of the template structures..
```commandline
python3 main.py --mode=0 -F=[Map_Path] -P=[Single_chain_structure_dir] -M=[protin_config_path] --config=[pipeline_config_file] --contour=[Contour_Level] --gpu=[GPU_ID] --resolution=[resolution]
```
[Map_Path] is the path of the input experimental cryo-EM map, [Single_chain_structure_dir] specifis the directory of all single-chain PDB files. You can also zip them into a zip/.tar/.tar.gz file here to pass the zip file path here. [protin_config_path] is the text file path that records the protein single chain PDB name and corressponding chains, [pipeline_config_file] is the pipeline's parameter configuration file, saved in ``config`` directory; [Contour_Level] is the map density threshold to remove outside regions to save processing time (suggested to use half author recommended contour level), [GPU_ID] specifies the gpu used for inference. [resolution] specified the map resolution, where 0-2A will skip the diffusion model. Therefore, you can use an approximate resolution value here.

Example of PDB config file
```commandline
tmp1.pdb A B
tmp2.pdb C
tmp3.pdb D E
```
which indicates 5 single-chain structures are provided. tmp1.pdb is for identical chain A and B, tmp2.pdb is for chain C, tmp3.pdb is for identical chain D and chain E.
To obtain such template/alphafold predicted single-chain structure, please consider two options:

1. Please check <a href='https://alphafold.ebi.ac.uk/'>AlphaFold Database</a> for single-chain structure with UniProt ID.
2. Please simpy search <a href='https://www.ebi.ac.uk/Tools/sss/fasta/'>EBI Search Tool</a> aginst structure database to find most similar structures as templates for us to model protein complex. Here you can get similar (or even identical) experimental single-chain structures or AlphaFold predicted structures.

### Example Command
```commandline
python3 main.py --mode=0 -F=example/6824.mrc -P=example -M=example/input_info.txt --config=config/diffmodeler.json --contour=2 --gpu=0 --resolution=5.8
```
This is the example command with ``-P`` specify the directory of single-chain pdbs.

You can also use this command line to specify the zip file including all single-chain PDB files for ``-P``
```commandline
python3 main.py --mode=0 -F=example/6824.mrc -P=example/6824.zip -M=example/input_info.txt --config=config/diffmodeler.json --contour=2 --gpu=0 --resolution=5.8
```
</details>

<details>
<summary>Protein Structure Complex Modeling with sequence (EBI-search)</summary>

### Protein Structure Complex Modeling with sequence (EBI-search)
This is for DiffModeler running if you have map and corresponding sequence. It is fine to run if you only know some of the sequences.
<br><b>Please use our [server](https://em.kiharalab.org/algorithm/DiffModeler(seq)) if with more than 4 non-identical chains. </b> EBI's API is too slow to respond when you have many non-identical sequences. <br>
This mode can only support structure modeling with less than 4 non-identical chains because of EBI's search tool limit.
```commandline
python3 main.py --mode=1 -F=[Map_Path] -P=[fasta_path] --config=[pipeline_config_file] --contour=[Contour_Level] --gpu=[GPU_ID] --resolution=[resolution]
```
[Map_Path] is the path of the input experimental cryo-EM map, [fasta_path] specifis the path of sequence file with .fasta format. [pipeline_config_file] is the pipeline's parameter configuration file, saved in ``config`` directory; [Contour_Level] is the map density threshold to remove outside regions to save processing time (suggested to use half author recommended contour level), [GPU_ID] specifies the gpu used for inference. [resolution] specified the map resolution, where 0-2A will skip the diffusion model. Therefore, you can use an approximate resolution value here.
<b>Please update ``email`` field in config/diffmodeler.json to your email.</b> 

This is based on <a href='https://www.ebi.ac.uk/Tools/sss/fasta/'>EBI Search Tool</a> aginst structure database to find most similar structures as templates for us to model protein complex. 
<br><b>Please use our [server](https://em.kiharalab.org/algorithm/DiffModeler(seq)) if with more than 4 non-identical chains. </b> EBI's API is too slow to respond when you have many non-identical sequences.

Example of fasta file
```
>A,B,C,D
MATPAGRRASETERLLTPNPGYGTQVGTSPAPTTPTEEEDLRR
>E,F
VVTFREENTIAFRHLFLLGYSDGSDDTFAAYTQEQLYQ
```
For ID line, please only include the chain id without any other information. If multiple chains include the identical sequences, please use comma "," to split different chains.
<br> In this example, we have 6 chains in total, with A,B,C,D share the identical sequences and E,F share another identical sequences.

### Example Command
```commandline
python3 main.py --mode=1 -F=example/6824.mrc -P=example/6824.fasta --config=config/diffmodeler.json --contour=2 --gpu=0 --resolution=5.8
```
</details>

<details>
<summary>Protein Structure Complex Modeling with sequence (Local Sequence Database)</summary>

### Protein Structure Complex Modeling with sequence (Local Sequence Database)
This is for DiffModeler running if you have map and corresponding sequence. It is fine to run if you only know some of the sequences. If you have more 4 non-identical chains, you can set up local sequence database to run DiffModeler by yourself.
<br><b>This is the mode used in server. We highly recommend to directly user [server](https://em.kiharalab.org/algorithm/DiffModeler(seq)) </b>
### 1 Install Blast
Please follow the instructions in [NCBI website](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) to install Blast locally.
### 2 Download and install database
Download the processed database from https://huggingface.co/datasets/zhtronics/BLAST_RCSB_AFDB/tree/main.
You can also use command line
```commandline
wget https://huggingface.co/datasets/zhtronics/BLAST_RCSB_AFDB/resolve/main/data.tar.gz.aa
wget https://huggingface.co/datasets/zhtronics/BLAST_RCSB_AFDB/resolve/main/data.tar.gz.ab
wget https://huggingface.co/datasets/zhtronics/BLAST_RCSB_AFDB/resolve/main/data.tar.gz.ac
```
After downloading finished, unzip the database to ``data`` subdirectory under this project with following command
```commandline
cat data.tar.gz.aa data.tar.gz.ab data.tar.gz.ac >data.tar.gz
tar -xzvf data.tar.gz
```
### 3 Run DiffModeler
After configuring the environment, please run 
```commandline
python3 main.py --mode=2 -F=[Map_Path] -P=[fasta_path] --config=[pipeline_config_file] --contour=[Contour_Level] --gpu=[GPU_ID] --resolution=[resolution]
```
[Map_Path] is the path of the input experimental cryo-EM map, [fasta_path] specifis the path of sequence file with .fasta format. [pipeline_config_file] is the pipeline's parameter configuration file, saved in ``config`` directory; [Contour_Level] is the map density threshold to remove outside regions to save processing time (suggested to use half author recommended contour level), [GPU_ID] specifies the gpu used for inference. [resolution] specified the map resolution, where 0-2A will skip the diffusion model. Therefore, you can use an approximate resolution value here.

Example of fasta file
```
>A,B,C,D
MATPAGRRASETERLLTPNPGYGTQVGTSPAPTTPTEEEDLRR
>E,F
VVTFREENTIAFRHLFLLGYSDGSDDTFAAYTQEQLYQ
```
For ID line, please only include the chain id without any other information. If multiple chains include the identical sequences, please use comma "," to split different chains.
<br> In this example, we have 6 chains in total, with A,B,C,D share the identical sequences and E,F share another identical sequences.


### Example Command
```commandline
python3 main.py --mode=2 -F=example/6824.mrc -P=example/6824.fasta --config=config/diffmodeler.json --contour=2 --gpu=0 --resolution=5.8
```

</details>

## Example

<details>

### Input File
Cryo-EM map with mrc format. 
AlphaFold/Template single-chain structure and information file to indicate the path.
Our example input can be found [here](https://github.com/kiharalab/DiffModeler/tree/master/example)

### Output File 
DiffModeler.cif: a CIF file that records the final modeled protein complex structure.
Our example output can be found [here](https://kiharalab.org/emsuites/diffmodelder_example/output). All the intermediate results are also kept here. 

</details>
