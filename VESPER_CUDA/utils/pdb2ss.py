import os
import pathlib
import shutil
import tempfile

import mrcfile
# from TEMPy.maps.map_parser import MapParser
# from TEMPy.protein.structure_blurrer import StructureBlurrer
# from TEMPy.protein.structure_parser import PDBParser, mmCIFParser
from utils.pdb2vol import pdb2vol
import gemmi


# import biotite.structure.io as strucio
# import biotite.structure as struc


# def split_pdb_by_ss(pdb_path, output_dir):
#     array = strucio.load_structure(pdb_path)
#     residues = struc.get_residues(array)[0]
#     sse = struc.annotate_sse(array)
#
#     # get res ids for each ss class
#     a_res = residues[sse == "a"]
#     b_res = residues[sse == "b"]
#     c_res = residues[sse == "c"]
#
#     # create ss mask by residue id
#     a_mask = [(True if res_id in a_res else False) for res_id in array.res_id]
#     b_mask = [(True if res_id in b_res else False) for res_id in array.res_id]
#     c_mask = [(True if res_id in c_res else False) for res_id in array.res_id]
#
#     # apply mask to array
#     arr_a = array[a_mask]
#     arr_b = array[b_mask]
#     arr_c = array[c_mask]
#
#     os.makedirs(output_dir, exist_ok=True)
#     pdb_path = pathlib.Path(pdb_path)
#
#     strucio.save_structure(os.path.join(output_dir, f"{pdb_path.stem}_ssA.pdb"), arr_a)
#     strucio.save_structure(os.path.join(output_dir, f"{pdb_path.stem}_ssB.pdb"), arr_b)
#     strucio.save_structure(os.path.join(output_dir, f"{pdb_path.stem}_ssC.pdb"), arr_c)


def split_cif_by_ss(cif_path, output_dir):
    """
    Split the CIF file by secondary structure and save the resulting structures as separate CIF files.

    Args:
        cif_path (str): The path to the input CIF file.
        output_dir (str): The directory where the output CIF files will be saved.
    """
    cif_data_block = gemmi.cif.read(cif_path).sole_block()
    structure = gemmi.make_structure_from_block(cif_data_block)[0]
    residue_coords_list = []
    residue_list = []
    for chain in structure:
        for residue in chain:
            curr_res_coords = []
            curr_atom_types = set()
            for atom in residue:
                if atom.name in ["N", "CA", "C", "O"]:
                    curr_res_coords.append(np.array((atom.pos.x, atom.pos.y, atom.pos.z)))
                    curr_atom_types = curr_atom_types | {atom.name}
            if len(curr_res_coords) == 4 and curr_atom_types == {"N", "CA", "C", "O"}:
                residue_coords_list.append(np.array(curr_res_coords))
                residue_list.append(residue)
    residue_coords = np.array(residue_coords_list)
    dssp = pydssp.assign(residue_coords, out_type="c3")

    st_a = gemmi.Structure()
    st_b = gemmi.Structure()
    st_c = gemmi.Structure()
    model_a = gemmi.Model("0")
    model_b = gemmi.Model("0")
    model_c = gemmi.Model("0")
    chain_a = gemmi.Chain("ssA")
    chain_b = gemmi.Chain("ssB")
    chain_c = gemmi.Chain("ssC")

    # save cif files
    for res, ss in zip(residue_list, dssp):
        if ss == "H":
            chain_a.add_residue(res)
        elif ss == "E":
            chain_b.add_residue(res)
        elif ss == "-":
            chain_c.add_residue(res)

    model_a.add_chain(chain_a)
    model_b.add_chain(chain_b)
    model_c.add_chain(chain_c)
    st_a.add_model(model_a)
    st_b.add_model(model_b)
    st_c.add_model(model_c)

    groups = gemmi.MmcifOutputGroups(True)
    groups.group_pdb = True
    doc_A = st_a.make_mmcif_document(groups)
    doc_B = st_b.make_mmcif_document(groups)
    doc_C = st_c.make_mmcif_document(groups)

    pdb_path = pathlib.Path(cif_path)

    doc_A.write_file(os.path.join(output_dir, f"{pdb_path.stem}_ssA.cif"))
    doc_B.write_file(os.path.join(output_dir, f"{pdb_path.stem}_ssB.cif"))
    doc_C.write_file(os.path.join(output_dir, f"{pdb_path.stem}_ssC.cif"))


import pydssp
from Bio.PDB import PDBParser, Select, MMCIFIO
import numpy as np


def split_pdb_by_ss(pdb_path, output_dir):
    """
    Split a PDB file by secondary structure and save the resulting structures as CIF files.

    Args:
        pdb_path (str): Path to the PDB file.
        output_dir (str): Directory to save the resulting CIF files.
    """
    parser = PDBParser()
    structure = parser.get_structure("PDB1", pdb_path)
    residues = list(structure.get_residues())
    residue_coords_list = []

    for res in residues:
        curr_res_coords = []
        curr_atom_types = set()
        for atom in res:
            if atom.name in ["N", "CA", "C", "O"]:
                curr_res_coords.append(np.array(atom.coord))
                curr_atom_types = curr_atom_types | {atom.name}
        if len(curr_res_coords) == 4 and curr_atom_types == {"N", "CA", "C", "O"}:
            residue_coords_list.append(np.array(curr_res_coords))
    residue_coords = np.array(residue_coords_list)

    dssp = pydssp.assign(residue_coords, out_type="c3")

    res_a = [res for res, ss in zip(residues, dssp) if ss == "H"]
    res_b = [res for res, ss in zip(residues, dssp) if ss == "E"]
    res_c = [res for res, ss in zip(residues, dssp) if ss == "-"]

    # save pdb files
    io = MMCIFIO()

    class ResidueSelect(Select):
        def __init__(self, res_list):
            self.res_list = res_list

        def accept_residue(self, res):
            if res in self.res_list:
                return True
            else:
                return False

    io.set_structure(structure)
    pdb_path = pathlib.Path(pdb_path)
    io.save(os.path.join(output_dir, f"{pdb_path.stem}_ssA.cif"), ResidueSelect(res_a))
    io.save(os.path.join(output_dir, f"{pdb_path.stem}_ssB.cif"), ResidueSelect(res_b))
    io.save(os.path.join(output_dir, f"{pdb_path.stem}_ssC.cif"), ResidueSelect(res_c))


def gen_simu_map(file_path, res, output_path, ref_dens_map=None):
    """
    The gen_simu_map function takes a PDB file and generates a simulated map from it.

    :param file_path: Specify the path to the pdb file
    :param res: Set the resolution of the simulated map
    :param output_path: Specify the path to where the simulated map will be saved
    :param ref_dens_map: Specify a density map to use as a reference for output dimensions
    :return: A simulated map based on a pdb file
    """

    is_cif = file_path.split(".")[-1] == "cif"

    # check number of atoms in pdb file, if none, return new map with dimensions of densMap\

    if is_cif:
        from Bio.PDB import MMCIF2Dict

        cif_dict = MMCIF2Dict.MMCIF2Dict(file_path)

        if "_atom_site.label_atom_id" in cif_dict:
            atom_site_data = cif_dict["_atom_site.label_atom_id"]
            atom_count = len(atom_site_data)
        elif "_atom_site.auth_atom_id" in cif_dict:
            atom_site_data = cif_dict["_atom_site.auth_atom_id"]
            atom_count = len(atom_site_data)
        else:
            atom_count = 0
    else:
        with open(file_path) as f:
            atom_count = 0
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1

    if atom_count == 0:
        # handle no atoms in ss category
        if ref_dens_map:
            # create new map with dimensions of ref_dens_map but all zero density values
            with mrcfile.open(ref_dens_map, permissive=True) as mrc:
                with mrcfile.new(output_path, overwrite=True) as mrc_new:
                    mrc_new.set_data(np.zeros(mrc.data.shape, dtype=np.float32))
                    mrc_new.voxel_size = mrc.voxel_size
                    mrc_new.update_header_from_data()
                    mrc_new.header.nxstart = mrc.header.nxstart
                    mrc_new.header.nystart = mrc.header.nystart
                    mrc_new.header.nzstart = mrc.header.nzstart
                    mrc_new.header.origin = mrc.header.origin
                    mrc_new.header.mapc = mrc.header.mapc
                    mrc_new.header.mapr = mrc.header.mapr
                    mrc_new.header.maps = mrc.header.maps
                    mrc_new.update_header_stats()
                    mrc_new.flush()  # write to disk
        else:
            raise Exception("No atoms in PDB file and no density map specified.")
    else:
        # ref_dens_map = MapParser.readMRC(ref_dens_map) if ref_dens_map else None
        # sb = StructureBlurrer()
        pdb_path = os.path.abspath(file_path)
        # output_path = os.path.abspath(output_path)
        # if file_path.split(".")[-1] == "cif":
        #     st = mmCIFParser.read_mmCIF_file(pdb_path, hetatm=True)
        # elif file_path.split(".")[-1] == "pdb":
        #     st = PDBParser.read_PDB_file("pdb1", pdb_path)
        # else:
        #     raise Exception("Make sure the input file is a PDB or mmCIF file.")
        # simu_map = sb.gaussian_blur_real_space(st, res, densMap=ref_dens_map)
        # simu_map.write_to_MRC_file(output_path)
        pdb2vol(pdb_path, res, output_path, ref_map=ref_dens_map)
        assert os.path.exists(output_path), "Simulated map not generated."


def gen_npy(pdb_path, sample_res, npy_path=None, verbose=False):
    if verbose:
        print("Combining MRC files into a Numpy array...")

    # get stem of pdb file
    pdb_stem = str(pathlib.Path(pdb_path).stem)
    tmp_dir = f"./{tempfile.gettempdir()}/{pdb_stem}/"
    try:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
    except:
        pass
    os.makedirs(tmp_dir, exist_ok=True)

    pdb_dir = os.path.join(tmp_dir, "pdb")
    simu_mrc_dir = os.path.join(tmp_dir, "simu_mrc")

    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs(simu_mrc_dir, exist_ok=True)

    if pdb_path.split(".")[-1] == "cif":
        split_cif_by_ss(pdb_path, pdb_dir)
    elif pdb_path.split(".")[-1] == "pdb":
        split_pdb_by_ss(pdb_path, pdb_dir)
    else:
        raise Exception("Make sure the input file is a PDB or mmCIF file.")

    pdb_simu_map_path = os.path.join(tmp_dir, f"{pdb_stem}_simu_map.mrc")

    gen_simu_map(pdb_path, sample_res, pdb_simu_map_path, ref_dens_map=None)

    for file in os.listdir(pdb_dir):
        filename = str(pathlib.Path(file).stem)
        if "ssA" in file or "ssB" in file or "ssC" in file:
            gen_simu_map(
                os.path.join(pdb_dir, file),
                sample_res,
                os.path.join(simu_mrc_dir, filename + ".mrc"),
                ref_dens_map=pdb_simu_map_path,
            )

    with mrcfile.open(pdb_simu_map_path) as mrc:
        dims = mrc.data.shape

        arr = np.zeros((4, dims[0], dims[1], dims[2]))

        for file in os.listdir(simu_mrc_dir):
            if "ssC" in file:
                arr[0] = mrcfile.open(os.path.join(simu_mrc_dir, file)).data.copy()
            elif "ssB" in file:
                arr[1] = mrcfile.open(os.path.join(simu_mrc_dir, file)).data.copy()
            elif "ssA" in file:
                arr[2] = mrcfile.open(os.path.join(simu_mrc_dir, file)).data.copy()

        arr = np.transpose(arr, (1, 2, 3, 0))

    # get stem of pdb file
    pdb_path = pathlib.Path(pdb_path)
    pdb_stem = str(pdb_path.stem)

    if npy_path:
        os.makedirs(npy_path, exist_ok=True)
        save_pth = os.path.join(npy_path, pdb_stem + "_prob.npy")
        np.save(save_pth, arr)
        print("Numpy array saved to: " + save_pth)

    if verbose:
        print("size of array: ", arr.shape)
        # print min, max, mean, std
        for idx, ss_class in enumerate(["Coil", "Beta", "Alpha"]):
            curr_arr = arr[..., idx]
            print(
                f"{ss_class}: ",
                "Count: ", np.count_nonzero(curr_arr),
                "Min: ", np.min(curr_arr),
                "Max: ", np.max(curr_arr),
                "Mean: ", np.mean(curr_arr),
                "Std: ", np.std(curr_arr),
            )

    return arr
