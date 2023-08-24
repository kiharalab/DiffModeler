
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

## Free Online Server: https://em.kiharalab.org/algorithm/DiffModeler

## Introduction
Cryogenic electron microscopy (cryo-EM) has been widely employed in experimental settings to determine multi-chain protein complexes, but modeling accuracy greatly diminishes when resolution decreases. At intermediate resolutions of 5-10 Å, even template-based structure fitting presents significant challenges. To tackle this issue, we introduce DiffModeler, a fully automated protein complex structure modeling method that leverages a diffusion model for backbone tracing and structure fitting with AlphaFold predicted single-chain structure. In extensive testing on cryo-EM maps at intermediate resolution, DiffModeler showcased remarkably accurate structure modeling, surpassing existing methods significantly. Notably, we successfully modeled a protein complex consisting of 47 chains, comprising 13,462 residues, with an impressive TM-Score of 0.9. We also further benchmarked DiffModeler for maps at low resolution of 10-20 Å and validated its generalizability with plausible performances. 

## Overall Protocol 
1) Backbone tracing from cryo-EM maps at intermediate resolution via diffusion model. 
2) Single-chain structure prediction by AlphaFold. 
3) Single-chain structure fitting using VESPER. 
4) Protein complex modeling by assembling algorithms. 

## Pre-required software
### Required 
Python 3 : https://www.python.org/downloads/   
### Optional
Pymol (for structure visualization): https://pymol.org/2/    
Chimera (for map visualization): https://www.cgl.ucsf.edu/chimera/download.html  

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

## Usage
```commandline
usage: main.py [-h] --mode MODE [-F F] [-M M] [--config CONFIG] [--gpu GPU] [--output OUTPUT]
               [--contour CONTOUR]

options:
  -h, --help            show this help message and exit
  --mode MODE           control mode
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
### 1. Protein Structure Complex Modeling
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
## Example
### Input File
Cryo-EM map with mrc format. 
AlphaFold/Template single-chain structure and information file to indicate the path.
Our example input can be found [here](https://github.com/kiharalab/DiffModeler/tree/master/example)

### Output File 
DiffModeler.cif: a CIF file that records the final modeled protein complex structure.
Our example output can be found [here](https://kiharalab.org/emsuites/diffmodelder_example/output). All the intermediate results are also kept here. 

## Benchmark Dataset
All input and output of the benchmarked datasets are maintained [here](https://kiharalab.org/emsuites/diffmodelder_benchmark)
