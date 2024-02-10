import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from TEMPy.maps.map_parser import MapParser
import TEMPy.math.vector as Vector

import os


def save_rotated_pdb(input_pdb, rot_mtx, real_trans, save_path, model_num):
    pdbio = PDBIO()
    # check input file format
    parser = None
    if input_pdb.split(".")[-1] == "pdb":
        parser = PDBParser(QUIET=True)
    elif input_pdb.split(".")[-1] == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise Exception("Input PDB/mmCIF file format not supported.")

    structure = parser.get_structure("target_pdb", input_pdb)
    structure.transform(rot_mtx, real_trans)

    pdbio.set_structure(structure)
    pdbio.save(save_path)

    # add model number
    with open(save_path, "r") as f:
        lines = f.readlines()

    with open(save_path, "w") as f:
        f.write(f"MODEL        {model_num}\n")
        for line in lines:
            f.write(line)
        f.write("ENDMDL\n")


def save_rotated_mrc(mrc_path, rot_vec, real_trans, save_path):

    emap = MapParser.readMRC(mrc_path)
    center = Vector.Vector(0, 0, 0)

    # get max distance of diagonal
    max_dist = int(np.ceil(np.sqrt(emap.x_size() ** 2 + emap.y_size() ** 2 + emap.z_size() ** 2)))

    # pad map to make sure it is big enough
    emap.pad_map(max_dist, max_dist, max_dist)

    emap = emap.map_rotate_by_axis_angle(1, 0, 0, rot_vec[0], center)
    emap = emap.map_rotate_by_axis_angle(0, 1, 0, rot_vec[1], center)
    emap = emap.map_rotate_by_axis_angle(0, 0, 1, rot_vec[2], center)

    # add translation to origin

    new_orig = emap.origin + real_trans
    emap.change_origin(new_orig[0], new_orig[1], new_orig[2])

    # map = map.translate(real_trans[0], real_trans[1], real_trans[2])
    emap.write_to_MRC_file(save_path)


def save_vec_as_pdb(
    origin,
    sampled_mrc_vec,
    sampled_mrc_data,
    score_arr,
    score,
    sample_width,
    trans,
    save_path,
    rank,
):
    dim = sampled_mrc_data.shape[0]

    origin = np.array([origin[0], origin[1], origin[2]])
    trans = np.array(trans)

    trans = np.where(trans > 0.5 * dim, trans - dim, trans)

    add = origin - trans * sample_width

    natm = 1
    nres = 1

    with open(save_path, "w") as pdb_file:
        non_zero_dens_index = np.transpose(np.nonzero(sampled_mrc_data))
        for position in non_zero_dens_index:
            real_coord = position * sample_width + add  # center of voxel
            atom_header = f"ATOM{natm:>7d}  CA  ALA{nres:>6d}    "
            atom_content = f"{real_coord[0]:8.3f}{real_coord[1]:8.3f}{real_coord[2]:8.3f}{1.0:6.2f}{score_arr[position[0], position[1], position[2]]:6.2f}"
            pdb_file.write(atom_header + atom_content + "\n")
            natm += 1

            real_coord = (
                position + sampled_mrc_vec[position[0]][position[1]][position[2]]
            ) * sample_width + add  # center of voxel plus unit vector
            atom_header = f"ATOM{natm:>7d}  CB  ALA{nres:>6d}    "
            atom_content = f"{real_coord[0]:8.3f}{real_coord[1]:8.3f}{real_coord[2]:8.3f}{1.0:6.2f}{score_arr[position[0], position[1], position[2]]:6.2f}"
            pdb_file.write(atom_header + atom_content + "\n")
            natm += 1
            nres += 1


def save_map_as_pdb(origin, data, sample_width, trans, save_path):
    dim = data.shape[0] # dimension of map
    origin = np.array([origin[0], origin[1], origin[2]])
    trans = np.array(trans)
    trans = np.where(trans > 0.5 * dim, trans - dim, trans)
    offset = origin - trans * sample_width

    natm = 1    # atom number
    nres = 1    # residue number

    with open(save_path, "w") as pdb_file:
        non_zero_dens_index = np.transpose(np.nonzero(data))
        for position in non_zero_dens_index:
            real_coord = position * sample_width + offset
            atom_header = f"ATOM{natm:>7d}  CA  ALA{nres:>6d}    "
            atom_content = f"{real_coord[0]:8.3f}{real_coord[1]:8.3f}{real_coord[2]:8.3f}{1.0:6.2f}{data[position[0], position[1], position[2]]:6.2f}"
            pdb_file.write(atom_header + atom_content + "\n")
            natm += 1
            nres += 1