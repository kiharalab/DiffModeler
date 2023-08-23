
# DiffModeler
<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/DiffModeler-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>  

DiffModeler is a computational tool using a diffusion model to automatically build full protein complex structure from cryo-EM maps at intermediate and low resolution.  

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

### 3. Build dependencies.   
### 3.1 Install with anaconda
Install anaconda from https://www.anaconda.com/download#downloads.
Then create the environment via
```commandline
conda env create -f environment.yml
```

Each time when you want to run this software, simply activate the environment by
```
conda activate DiffModeler
conda deactivate(If you want to exit) 
```
#### 3.1 Install with pip and python.
##### 3.1.1[`install pip`](https://pip.pypa.io/en/stable/installing/).
##### 3.1.2  Install dependency in command line.
```
pip3 install -r requirements.txt --user
```
If you encounter any errors, you can install each library one by one:
```
pip3 install biopython
pip3 install numpy
pip3 install numba
pip3 install scipy
pip3 install mrcfile
pip3 install torch==1.9.0
pip3 install tqdm
pip3 install progress
pip3 install scikit-learn
pip3 install pyFFTW
pip3 install BioTEMPy>=2.0.0
```

## Usage
```commandline
usage: main.py [-h] --mode MODE [-F F] [-M M] [--config CONFIG] [--gpu GPU] [--output OUTPUT]
               [--contour CONTOUR]

options:
  -h, --help         show this help message and exit
  --mode MODE        control mode
  -F F               input map path
  -M M               txt file path which records protein information
  --config CONFIG    specifying the config path
  --gpu GPU          specify the gpu we will use
  --output OUTPUT    Output directory
  --contour CONTOUR  Contour level for input map, suggested 0.5*[author_contour]. (Float),
                     Default value: 0.0
```
### 1. Protein Structure Complex Modeling
```commandline
python3 main.py --mode=0 -F=[Map_Path] -M=[protin_config_path] --config=[pipeline_config_file] --contour=[Contour_Level] --gpu=[GPU_ID]
```
[Map_Path] is the path of the input experimental cryo-EM map, [protin_config_path] is the text file path that records the protein single chain path, [pipeline_config_file] is the pipeline's parameter configuration file, saved in config directory; [Contour_Level] is the map density threshold to remove outside regions to save processing time (suggested to use half author recommended contour level), [GPU_ID] specifies the gpu used for inference.

### Example Command
```commandline
python3 main.py --mode=0 -F=example/6824.mrc -M=example/input_info.txt --config=config/diffmodeler.json --contour=2 --gpu=0
```
## Example
### Input File
Cryo-EM map with mrc format. 
AlphaFold/Template single-chain structure and information file to indicate the path.
Our example input can be found [here](https://github.com/kiharalab/DiffModeler/tree/master/example)

### Output File 
DiffModeler.cif: a CIF file that records the final modeled protein complex structure.
Our example output can be found [here](). All the intermediate results are also kept here. 
