import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from TEMPy.maps.map_parser import MapParser
import TEMPy.math.vector as Vector

import os


def save_rotated_pdb(input_pdb, rot_mtx, real_trans, save_path):
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
    file_path,
    rank,
):
    dim = sampled_mrc_data.shape[0]

    # filename = "R_{:02d}-S_{:.3f}.pdb".format(rank, score).replace(" ", "_")
    #
    # filepath = os.path.join(folder_path, filename)

    origin = np.array([origin[0], origin[1], origin[2]])
    trans = np.array(trans)

    if trans[0] > 0.5 * dim:
        trans[0] -= dim
    if trans[1] > 0.5 * dim:
        trans[1] -= dim
    if trans[2] > 0.5 * dim:
        trans[2] -= dim

    add = origin - trans * sample_width

    natm = 1
    nres = 1

    pdb_file = open(file_path, "w")
    non_zero_dens_index = np.transpose(np.nonzero(sampled_mrc_data))
    for idx in non_zero_dens_index:
        tmp = idx * sample_width + add
        atom_header = "ATOM{:>7d}  CA  ALA{:>6d}    ".format(natm, nres)
        atom_content = "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format(
            tmp[0], tmp[1], tmp[2], 1.0, score_arr[idx[0], idx[1], idx[2]]
        )
        pdb_file.write(atom_header + atom_content + "\n")
        natm += 1

        tmp = (idx + sampled_mrc_vec[idx[0]][idx[1]][idx[2]]) * sample_width + add
        atom_header = "ATOM{:>7d}  CB  ALA{:>6d}    ".format(natm, nres)
        atom_content = "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format(
            tmp[0], tmp[1], tmp[2], 1.0, score_arr[idx[0], idx[1], idx[2]]
        )
        pdb_file.write(atom_header + atom_content + "\n")
        natm += 1
        nres += 1
