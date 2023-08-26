from Bio.PDB import PDBParser
import os
def read_proteinpdb_data(input_cif_path,run_type=0,atom_cutoff=1):
    """
    :param input_cif_path:
    :param run_type: 0: backbone mode, 1: all atom mode
    :return:
    """
    cif_name = os.path.split(input_cif_path)[1][:-4]
    parser = PDBParser()
    structure = parser.get_structure(cif_name, input_cif_path)
    residue_list = ["ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU", "ASP", "GLU", "LYS", "ARG", "SER", "THR", "TYR",
                    "HIS", "CYS", "ASN", "TRP", "GLN", "GLY"]
    dna_rna_set = {"A":0, "U":1, "T":1, "C":2, "G":3,
                   "DA":0, "DU":1, "DT":1, "DC":2, "DG":3}
    residue_dict = {}
    for k in range(len(residue_list)):
        res_name = residue_list[k]
        residue_dict[res_name]=k+1
    #atom class: N:1 CA:2 C:3 O:4 CB:5 others: 6
    if run_type==0:
        atom_map_set=['N',"CA","C"]
        atom_map_set= set(atom_map_set)
    elif run_type==2:
        atom_map_set=["CA"]
        atom_map_set= set(atom_map_set)

    Information_Dict = {}
    check_id = 0
    residue_id =0
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                res_name=residue.get_resname()
                if res_name not in residue_dict:
                    print("res name %s not included"%res_name)
                    continue
                else:
                    res_id = residue_dict[res_name]

                for atom in residue.get_list():
                    atom_name = atom.get_fullname().replace(" ","")
                    atom_coord = atom.get_coord()
                    format_coord = []
                    atom_coord = str(atom_coord)
                    atom_coord=atom_coord.replace("[","")
                    atom_coord=atom_coord.replace("]","")
                    atom_coord_split = atom_coord.split()
                    for k in range(3):
                        format_coord.append(float(atom_coord_split[k]))
                    #atom_radius = atom.get_radius()
                    # if "H" in atom_name:
                    #     try:
                    #         print("hydrogen %s is not considered in our case"
                    #           " with radius %.5f"%(atom_name,atom_radius))
                    #     except:
                    #         print("hydrogen %s is not considered in our case"% (atom_name))
                    #     continue
                    tmp_dict = {}
                    tmp_dict['coord'] = format_coord
                    if (run_type==0 or run_type==2) and atom_name not in atom_map_set:
                        continue



                    tmp_dict['atom'] = 1
                    #if atom_radius is None:
                        #print("resetting to 2 since atom radius is none")
                    atom_radius = atom_cutoff
                    tmp_dict['radius'] = atom_radius

                    if res_name in residue_dict:
                        res_id = residue_dict[res_name]
                        check_nuc = res_id # set as protein
                    else:
                        if "H" not in atom_name:
                            print("%s not in nuc dict" % atom_name)
                            check_nuc = 0
                    tmp_dict['nuc'] = check_nuc
                    tmp_dict['res_id'] = residue_id
                    Information_Dict[check_id]=tmp_dict

                    check_id +=1
                residue_id+=1
    return Information_Dict
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.mmcifio import MMCIFIO
def rename_chains(pdb_file, new_chain_name, output_file):
    """
    Reads a PDB file, renames all chains to the specified name, and writes the new PDB file.

    Parameters:
    pdb_file (str): The name of the input PDB file to read.
    new_chain_name (str): The new chain name to set for all chains in the structure.
    output_file (str): The name of the output PDB file to write.
    """
    # Parse the input PDB file
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)

    # Rename all chains to the specified name
    for chain in structure.get_chains():
        chain.id = new_chain_name

    # Write the new PDB file with the renamed chains
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

# def rename_chains_cif(pdb_file, new_chain_name, output_file):
#     """
#     Reads a PDB file, renames all chains to the specified name, and writes the new PDB file.
#
#     Parameters:
#     pdb_file (str): The name of the input PDB file to read.
#     new_chain_name (str): The new chain name to set for all chains in the structure.
#     output_file (str): The name of the output PDB file to write.
#     """
#     # Parse the input PDB file
#     parser = PDBParser()
#     structure = parser.get_structure('pdb', pdb_file)
#
#     # Rename all chains to the specified name
#     for chain in structure.get_chains():
#         chain.id = new_chain_name
#
#     # Write the new PDB file with the renamed chains
#     io = MMCIFIO() #PDBIO()
#     io.set_structure(structure)
#     io.save(output_file)
def rename_chains_cif(pdb_file,  new_chain_id,cif_file):
    # Read the PDB file
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()

    cif_lines = []
    current_model = 0
    atom_serial = 0

    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            x=float(line[30:38])
            y=float(line[38:46])
            z=float(line[46:54])
            line = line[:21] + new_chain_id + line[22:30]+" %.3f %.3f %.3f "%(x,y,z)+line[54:]
            atom_serial += 1
            cif_lines.append(line)

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
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")
        f.write("_atom_site.type_symbol\n")

        for line in cif_lines:
            f.write(line)
def rename_chains_naive(pdb_file_path, new_chain_name, output_file_path):
    with open(pdb_file_path, 'r') as f:
        lines = f.readlines()

    with open(output_file_path, 'w') as f:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
               # chain_id = line[21]

                line = line[:21] + new_chain_name + line[22:]
                f.write(line)
def replace_chain_ids(pdb_file_path, chain_dict, output_file_path):
    """
    Read a PDB file, replace each chain ID based on a dictionary mapping, and write the modified PDB file to disk.

    Args:
    pdb_file_path (str): The path to the input PDB file.
    chain_dict (dict): A dictionary mapping old chain IDs to new chain IDs.
    output_file_path (str): The path to the output PDB file.

    Returns:
    None
    """
    with open(pdb_file_path, 'r') as f:
        lines = f.readlines()

    with open(output_file_path, 'w') as f:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                if chain_id in chain_dict:
                    new_chain_id = chain_dict[chain_id]
                    line = line[:21] + new_chain_id + line[22:]
                f.write(line)
def count_atom_line(file_path):
    count=0
    with open(file_path,'r') as file:
        for line in file:
            if len(line)>5 and line[:4]=="ATOM":
                count+=1
    return count

def reindex_cif(final_pdb_path,final_cif_path):
    begin_check=False
    block_list=[]
    with open(final_pdb_path,'r') as rfile:
        for line in rfile:
            if "loop_" in line:
                begin_check=True
                continue

            if begin_check and "_atom_site" in line:
                block_list.append(line.strip("\n").replace(" ",""))
                continue
            if begin_check and "_atom_site" not in line:
                begin_check=False

    atom_ids = block_list.index('_atom_site.id')
    try:
        seq_ids = block_list.index('_atom_site.label_seq_id')
    except:
        seq_ids = block_list.index('_atom_site.auth_seq_id')
    atom_id=1
    seq_id=1
    prev_seq_id=None
    with open(final_pdb_path,'r') as rfile:
        with open(final_cif_path,'w') as wfile:
            for line in rfile:
                if len(line)>4 and line[:4]=="ATOM":
                    split_info=line.strip("\n").split()
                    current_seq_id = int(split_info[seq_ids])
                    if prev_seq_id is not None and current_seq_id!=prev_seq_id:
                        seq_id+=1
                    for j,item in enumerate(split_info):
                        if j==atom_ids:
                            wfile.write("%d  "%atom_id)

                        elif j!=seq_ids:
                            wfile.write("%s  "%item)
                        else:
                            wfile.write("%d  "%seq_id)
                    wfile.write("\n")
                    prev_seq_id=current_seq_id
                    atom_id+=1
                else:
                    wfile.write(line)

def filter_chain_cif(input_cif,select_chain_name,output_cif):
    begin_check=False
    block_list=[]
    with open(input_cif,'r') as rfile:
        for line in rfile:
            if "loop_" in line:
                begin_check=True
                continue

            if begin_check and "_atom_site" in line:
                block_list.append(line.strip("\n").replace(" ",""))
                continue
            if begin_check and "_atom_site" not in line:
                begin_check=False
    atom_ids = block_list.index('_atom_site.id')
    try:
        seq_ids = block_list.index('_atom_site.label_seq_id')
    except:
        seq_ids = block_list.index('_atom_site.auth_seq_id')
    chain_ids = block_list.index("_atom_site.label_asym_id")
    atom_id=1
    seq_id=1
    prev_seq_id=None
    with open(input_cif,'r') as rfile:
        with open(output_cif,'w') as wfile:
            #writee header
            # wfile.write("chain_%s\n"%select_chain_name)
            # wfile.write("#\nloop_\n")
            # for block_id in block_list:
            #     wfile.write("%s\n"%block_id)
            for line in rfile:
                if len(line)>4 and line[:4]=="ATOM":
                    split_info=line.strip("\n").split()
                    current_chain = split_info[chain_ids]
                    if current_chain!=select_chain_name:
                        continue
                    current_seq_id = int(split_info[seq_ids])
                    if prev_seq_id is not None and current_seq_id!=prev_seq_id:
                        seq_id+=1
                    for j,item in enumerate(split_info):
                        if j==atom_ids:
                            wfile.write("%d  "%atom_id)

                        elif j!=seq_ids:
                            wfile.write("%s  "%item)
                        else:
                            wfile.write("%d  "%seq_id)
                    wfile.write("\n")
                    prev_seq_id=current_seq_id
                    atom_id+=1
                else:
                    wfile.write(line)

def cif2pdb(input_cif_path,final_pdb_path):
    begin_check=False
    block_list=[]
    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if "loop_" in line:
                begin_check=True
                continue

            if begin_check and "_atom_site" in line:
                block_list.append(line.strip("\n").replace(" ",""))
                continue
            if begin_check and "_atom_site" not in line:
                begin_check=False
    atom_ids = block_list.index('_atom_site.id')
    try:
        seq_ids = block_list.index('_atom_site.label_seq_id')
    except:
        seq_ids = block_list.index('_atom_site.auth_seq_id')
    chain_ids = block_list.index("_atom_site.label_asym_id")
    atom_type_ids = block_list.index("_atom_site.label_atom_id")
    res_name_ids = block_list.index("_atom_site.label_comp_id")
    x_ids = block_list.index("_atom_site.Cartn_x")
    y_ids = block_list.index("_atom_site.Cartn_y")
    z_ids = block_list.index("_atom_site.Cartn_z")
    with open(input_cif_path,'r') as rfile:
        with open(final_pdb_path,'w') as wfile:

            for line in rfile:
                if len(line)>4 and line[:4]=="ATOM":
                    split_info=line.strip("\n").split()
                    current_chain = split_info[chain_ids]
                    current_atom_index = int(split_info[atom_ids])
                    current_atom_name = split_info[atom_type_ids]
                    current_res_index = int(split_info[seq_ids])
                    current_res_name = split_info[res_name_ids]
                    current_x = float(split_info[x_ids])
                    current_y = float(split_info[y_ids])
                    current_z = float(split_info[z_ids])
                    wline=""
                    wline += "ATOM%7d %-4s %3s%2s%4d    " % (current_atom_index, current_atom_name,
                                                             current_res_name, current_chain,current_res_index)
                    wline = wline + "%8.3f%8.3f%8.3f%6.2f\n" % (current_x,current_y,current_z, 1.0)
                    wfile.write(wline)
def filter_chain_pdb(pdb_file_path, chain_name, output_file_path):
    """
    Read a PDB file, replace each chain ID based on a dictionary mapping, and write the modified PDB file to disk.

    Args:
    pdb_file_path (str): The path to the input PDB file.
    chain_dict (dict): A dictionary mapping old chain IDs to new chain IDs.
    output_file_path (str): The path to the output PDB file.

    Returns:
    None
    """
    with open(pdb_file_path, 'r') as f:
        lines = f.readlines()

    with open(output_file_path, 'w') as f:
        for line in lines:
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id!=chain_name:
                    continue
                f.write(line)
