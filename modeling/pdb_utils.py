
from ops.pdb_utils import reindex_cif,rename_chains_cif

def remove_hetatm_lines(input_file, output_file):
    """
    Remove HETATM lines from a PDB file and save the result to a new file.

    Args:
        input_file (str): Path to the input PDB file.
        output_file (str): Path to the output PDB file.
    """
    with open(input_file, 'r') as pdb_file:
        lines = pdb_file.readlines()

    filtered_lines = [line for line in lines if not line.startswith('HETATM')]

    with open(output_file, 'w') as output_file:
        output_file.writelines(filtered_lines)

def filter_backbone(input_pdb_path,backbone_pdb_path):
    backbone_list=["CA","C","N"]
    #support DNA/RNA template fitting
    backbone_list_drna=["OP3","P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
    backbone_list+=backbone_list_drna
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
def extract_residue_locations_slow(pdb_file):
    # Load the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Extract the locations of all residues
    residue_locations = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                #only extract ca positions
                atom_name = atom.get_name().replace(" ","")
                if atom_name=="CA" or atom_name=="P":
                    residue_locations.append(atom.get_coord())

    return np.array(residue_locations)
def extract_residue_locations(pdb_file):
    residue_locations = []
    with open(pdb_file) as f:
        for line in f:
            # CA for protein, P for dna/rna
            if line.startswith("ATOM") and (line[12:16].strip() == "CA" or line[12:16].strip() == "P"):
                residue_locations.append(
                    np.array(
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    )
                )

    return np.array(residue_locations)

def find_identical_chain(fit_dict,query_chain,visit_dict):
    for key in fit_dict:
        current_chain_list=fit_dict[key]
        if query_chain in current_chain_list:
            new_list = [chain_name for chain_name in current_chain_list if visit_dict[chain_name]==0]
            return new_list
    return None

def clean_fit_dict(fit_dict,query_chain):
    new_fit_dict = {}
    for key in fit_dict:
        current_chain_list=fit_dict[key]
        if query_chain in current_chain_list:
            current_chain_list.remove(query_chain)
            if len(current_chain_list)>0:
                new_fit_dict[key]= current_chain_list
        else:
            new_fit_dict[key]=current_chain_list
    return new_fit_dict
def find_chain_pdb_path(fit_dict,query_chain):
    """
    find the pdb path of specified chain
    :param fit_dict:
    :param query_chain:
    :return:
    """
    for key in fit_dict:
        current_chain_list=fit_dict[key]
        if query_chain in current_chain_list:
            return key
    return None
import torch

def remove_overlap_pdb(score_dict,fitted_pdb,clash_distance=3,ratio_cutoff=0.05,split_key=True):
    new_score_dict= {}
    current_atom_locations= extract_residue_locations(fitted_pdb)
    current_atom_locations = torch.from_numpy(current_atom_locations).cuda()
    for key in score_dict:
        if split_key==False:
            current_pdb_path = key
        else:
            current_pdb_path = key.split(",")[1]
        query_atom_locations = extract_residue_locations(current_pdb_path)
        query_atom_locations = torch.from_numpy(query_atom_locations).cuda()
        distance_array = torch.cdist(query_atom_locations,current_atom_locations, p=2)
        distance_array = torch.amin(distance_array,dim=1)
        ratio_close = (distance_array < clash_distance).sum().item() / len(query_atom_locations)#len(distance_array[distance_array<=clash_distance])/len(query_atom_locations)
        if ratio_close<ratio_cutoff:
            new_score_dict[key]=score_dict[key]
        else:
            print("%s remove because of the overlap too big"%(key),ratio_close)
    return new_score_dict
def calculate_overlap_score(score_dict,clash_dict,fitted_pdb,clash_distance=3,split_key=True):
    current_atom_locations= extract_residue_locations(fitted_pdb)
    current_atom_locations = torch.from_numpy(current_atom_locations).cuda()

    for key in score_dict:
        if split_key==False:
            current_pdb_path = key
        else:
            current_pdb_path = key.split(",")[1]
        query_atom_locations = extract_residue_locations(current_pdb_path)
        query_atom_locations = torch.from_numpy(query_atom_locations).cuda()
        distance_array = torch.cdist(query_atom_locations,current_atom_locations, p=2)
        distance_array = torch.amin(distance_array,dim=1)
        ratio_close = (distance_array < clash_distance).sum().item() / len(query_atom_locations)#len(distance_array[distance_array<=clash_distance])/len(query_atom_locations)
        if key not in clash_dict:
            clash_dict[key]=ratio_close
        else:
            clash_dict[key]+=ratio_close
    return clash_dict

def rename_chains_pdb2cif(pdb_file,  new_chain_id,cif_file):
    # Read the PDB file
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()

    cif_lines = []
    current_model = 0
    atom_serial = 0

    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # x=float(line[30:38])
            # y=float(line[38:46])
            # z=float(line[46:54])
            # line = line[:21] + new_chain_id + line[22:30]+" %.3f %.3f %.3f "%(x,y,z)+line[54:]
            atom_serial += 1
            new_line=""
            new_line += line[:4]+"\t"
            new_line += line[6:11]+"\t"
            atom_type = line[12:16].replace(" ","")[0]
            new_line += atom_type+"\t"
            new_line += line[12:16]+"\t"
            # new_line += line[16]+"\t"
            # alt-id
            new_line += "." + "\t"
            new_line += line[17:20]+"\t"
            new_line += new_chain_id+"\t"
            new_line += line[22:26]+"\t"
            #new_line += line[26]+"\t"
            new_line += "1" + "\t"
            new_line += "?" + "\t"
            new_line += line[30:38]+"\t"
            new_line += line[38:46]+"\t"
            new_line += line[46:54]+"\t"
            new_line += line[54:60]+"\t"
            new_line += line[60:66]+"\t"
            new_line += "\n"
            cif_lines.append(new_line)

        # if line.startswith('ENDMDL'):
        #     current_model += 1
        #     atom_serial = 0

    # Write the modified structure to a new CIF file
    with open(cif_file, 'w') as f:
        f.write("data_" + new_chain_id + "\n#\n")
        # f.write("_audit_creation_method     'Diffusion Fitting'\n")
        # f.write("_audit_creation_date       \n")
        # f.write("_audit_author_name         'Xiao Wang'\n")
        # f.write("_entry.id                   %s\n" % new_chain_id)
        # f.write("\n")
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        # new column: alternative position
        f.write("_atom_site.label_alt_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.pdbx_PDB_model_num\n")
        f.write("_atom_site.pdbx_PDB_ins_code\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")

        for line in cif_lines:
            f.write(line)

def collect_final_pdb(modeling_dir,chain_visit_dict):
    final_pdb_path = os.path.join(modeling_dir,"Final_unformated.cif")
    first_structure=True
    with open(final_pdb_path,'w') as wfile:
        for chain in chain_visit_dict:
            current_dir=os.path.join(modeling_dir,"iterative_%s"%chain)
            tmp_cif_path = os.path.join(current_dir,"final_rechain.cif")
            with open(tmp_cif_path,'r') as rfile:
                if first_structure:
                    for line in rfile:
                        wfile.write(line)
                else:
                    for line in rfile:
                        if len(line)>=4 and line[:4]=="ATOM":
                            wfile.write(line)
            first_structure=False
    final_cif_path = os.path.join(modeling_dir,"Final.cif")
    reindex_cif(final_pdb_path,final_cif_path)
    return final_cif_path

def collect_domain_pdb(modeling_dir,chain_visit_dict):
    final_pdb_path = os.path.join(modeling_dir,"DiffModeler_domain.cif")
    first_structure=True
    with open(final_pdb_path,'w') as wfile:
        for chain in chain_visit_dict:
            current_dir=os.path.join(modeling_dir,"iterative_%s"%chain)
            old_cif_path = os.path.join(current_dir,"final_rechain.cif")
            tmp_cif_path = os.path.join(current_dir,"final_domain2chain.cif")
            real_chain_id = chain.split("_")[0]
            chain_map_dict = {chain:real_chain_id}#can only be used since old chain id includes it
            rename_chains_cif(old_cif_path, chain_map_dict, tmp_cif_path)
            
            with open(tmp_cif_path,'r') as rfile:
                if first_structure:
                    for line in rfile:
                        wfile.write(line)
                else:
                    for line in rfile:
                        if len(line)>=4 and line[:4]=="ATOM":
                            wfile.write(line)
            first_structure=False

    return final_pdb_path


def swap_cif_occupancy_bfactor(input_cif_path, output_cif_path):
    """
    input_cif_path: str
        Path to the input CIF file.
    output_cif_path: str
        Path to the output CIF file.

    """
    with open(input_cif_path,'r') as rfile:
        with open(output_cif_path,'w') as wfile:
            for line in rfile:
                if line.startswith("_atom_site.occupancy"):
                    #_atom_site.B_iso_or_equiv
                    line = line.replace("occupancy","B_iso_or_equiv")
                elif line.startswith("_atom_site.B_iso_or_equiv"):
                    line = line.replace("B_iso_or_equiv","occupancy")
                wfile.write(line)
