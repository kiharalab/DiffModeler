


def filter_backbone(input_pdb_path,backbone_pdb_path):
    backbone_list=["CA","C","N"]
    with open(input_pdb_path,'r') as file:
        with open(backbone_pdb_path,'w') as wfile:
            for line in file:
                if line.startswith("ATOM"):
                    chain_name = line[21]
                    atom_name = line[12:16].replace(" ","")
                    if atom_name in backbone_list:
                        wfile.write(line)


from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from scipy.spatial.transform import Rotation as R
import os
import numpy as np
def generate_move_structure(input_pdb_path,output_dir,vesper_path,count_model):
    parser = PDBParser(QUIET=True)
    with open(vesper_path) as f:

        model_line_start = ["#0", "#1", "#2", "#3", "#4", "#5", "#6", "#7", "#8", "#9"]
        pdbio = PDBIO()
        for line in f:
            if line[:2] in model_line_start:
                model_num = str(int(line[1]) + 1)
                r_info = line.split()[2:5]
                t_info = line.split()[6:9]

                rot_vec = np.array([float(r_info[0].replace("(","")), float(r_info[1]), float(r_info[2].replace(")",""))])
                trans_vec = np.array([float(t_info[0].replace("(","")), float(t_info[1]), float(t_info[2].replace(")",""))])
                rotation = R.from_euler('xyz', rot_vec, degrees=True)
                rotation = rotation.inv().as_matrix()

                search_pdb = parser.get_structure("search", input_pdb_path)
                search_pdb.transform(rotation, trans_vec)
                pdbio.set_structure(search_pdb)
                pdbio.save(os.path.join(output_dir,"vesper_%d.pdb"%count_model))
                count_model+=1
    return count_model
