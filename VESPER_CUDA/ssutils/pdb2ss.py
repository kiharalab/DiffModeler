import os
import pathlib
import shutil
import subprocess
import urllib.request

import mrcfile
import numpy as np
import pandas as pd
from TEMPy.maps.map_parser import MapParser
from TEMPy.protein.structure_blurrer import StructureBlurrer
from TEMPy.protein.structure_parser import PDBParser


def donwload_and_compile_stride():
    """
    The donwload_and_compile_stride function downloads and compiles the stride program if it is not already present.

    :return: None
    """
    curr_dir = os.getcwd()
    ssutils_path = pathlib.Path(__file__).parent.resolve()
    if not os.path.exists(ssutils_path / "stride"):
        print("No Stride binary found, downloading Stride...")
        os.chdir(ssutils_path)
        urllib.request.urlretrieve(
            "http://webclu.bio.wzw.tum.de/stride/stride.tar.gz",
            ssutils_path / "stride.tar.gz",
        )
        print("Extracting Stride...")
        os.makedirs(ssutils_path / "stride_src", exist_ok=True)
        os.system(f"tar -zxf stride.tar.gz -C {ssutils_path / 'stride_src'}")
        os.system("rm stride.tar.gz")
        os.chdir("stride_src")
        print("Compiling Stride...")
        os.system("make")
        os.chdir(ssutils_path)
        shutil.copyfile(ssutils_path / "stride_src" / "stride", ssutils_path / "stride")
        os.chmod(ssutils_path / "stride", 0o755)
        shutil.rmtree(ssutils_path / "stride_src")
        if os.path.exists("./stride"):
            print("Stride downloaded and compiled!")
        else:
            print("Stride download and compile failed!")
    os.chdir(curr_dir)


def contains_number(s):
    return any(i.isdigit() for i in s)


def assign_ss(pdb_path, output_dir, verbose=False):
    """
    The assign_ss function takes a PDB file and assigns secondary structure to each residue.
    It does this by running the program Stride, which is included in the repository.
    The output of assign_ss is a text file containing one line for each residue in the PDB file, with information about that residue's secondary structure assignment.

    :param pdb_path: Specify the path to the pdb file
    :param output_dir: Specify the directory where the output file will be saved
    :param verbose: Print out the progress of the function
    :return: The path to the output file
    """
    if verbose:
        print("Assigning secondary structure using Stride...")

    os.makedirs(output_dir, exist_ok=True)

    pdb_path = os.path.abspath(pdb_path)
    filename = pathlib.Path(pdb_path).stem
    stride_path = pathlib.Path(__file__).parent.resolve() / "stride"
    # os.system('./stride "' + pdb_path + '" > ' + '"' + output_dir + filename + ".ss" + '"')
    subprocess.run(
        [str(stride_path), pdb_path],
        stdout=open(output_dir + filename + ".ss", "w"),
        stderr=subprocess.PIPE,
    )

    if verbose:
        print("SS assignment file save to: " + output_dir + filename + ".ss")

    return output_dir + filename + ".ss"


def split_pdb_by_ss(pdb_path, ss_path, output_dir, verbose=False):
    """
    The split_pdb_by_ss function takes in a PDB file and a secondary structure prediction file,
    and splits the PDB into three files based on the predicted secondary structure. The output is
    three new PDB files with names ending in _ssA.pdb, _ssB.pdb, and _ssC.pdb.

    :param pdb_path: Specify the path to the pdb file
    :param ss_path: Specify the path to the secondary structure file
    :param output_dir: Specify the directory where the pdb files will be saved
    :param verbose: More verbose output
    :return: None
    """
    if verbose:
        print("Splitting PDB file by secondary structure...")

    os.makedirs(output_dir, exist_ok=True)

    file_ss = open(ss_path, mode="r")
    file_pdb = open(pdb_path, mode="r")

    ss_lines = file_ss.readlines()
    pdb_lines = file_pdb.readlines()

    ss_lines_pred_a = []
    ss_lines_pred_b = []
    ss_lines_pred_c = []

    for line in ss_lines:
        if line.startswith("ASG"):
            entries = line.split()
            if entries[5] == "H" or entries[5] == "G" or entries[5] == "I":
                ss_lines_pred_a.append(entries[3])
            elif entries[5] == "B" or entries[5] == "E":
                ss_lines_pred_b.append(entries[3])
            else:
                ss_lines_pred_c.append(entries[3])

    pdb_lines_atoms = []
    residual_nums = []
    for line in pdb_lines:
        if line.startswith("ATOM"):
            pdb_lines_atoms.append(line)
            if not contains_number(line.split()[4]):
                residual_nums.append(line.split()[5])
            else:
                residual_nums.append("".join(filter(str.isdigit, line.split()[4])))

    data_list = pd.Series(pdb_lines_atoms)
    data_res = pd.Series(residual_nums, dtype=str)
    df = pd.concat((data_list, data_res), axis=1)
    df.columns = ["str", "res_num"]

    # Create a new directory is not exists
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    a_strs = df[df["res_num"].isin(ss_lines_pred_a)]
    b_strs = df[df["res_num"].isin(ss_lines_pred_b)]
    c_strs = df[df["res_num"].isin(ss_lines_pred_c)]

    pdb_path = pathlib.Path(pdb_path)

    with open(output_dir + pdb_path.stem + "_ssA.pdb", "w") as fp:
        for item in a_strs["str"].to_list():
            fp.write("%s" % item)

    with open(output_dir + pdb_path.stem + "_ssB.pdb", "w") as fp:
        for item in b_strs["str"].to_list():
            fp.write("%s" % item)

    with open(output_dir + pdb_path.stem + "_ssC.pdb", "w") as fp:
        for item in c_strs["str"].to_list():
            fp.write("%s" % item)

    if verbose:
        print("PDB files saved to: " + output_dir)


def gen_simu_map(pdb_path, res, output_path, densMap=False):
    """
    The gen_simu_map function takes a PDB file and generates a simulated map from it.

    :param pdb_path: Specify the path to the pdb file
    :param res: Set the resolution of the simulated map
    :param output_path: Specify the path to where the simulated map will be saved
    :param densMap: Specify a density map to use as a reference for output dimensions
    :return: A simulated map based on a pdb file
    """
    sb = StructureBlurrer()
    pdb_path = os.path.abspath(pdb_path)
    output_path = os.path.abspath(output_path)
    pdb = PDBParser.read_PDB_file("pdb1", pdb_path)
    simu_map = sb.gaussian_blur_real_space(pdb, res, densMap=densMap)
    simu_map.write_to_MRC_file(output_path)


def gen_npy(pdb_path, sample_res, npy_path=None, verbose=False):
    if verbose:
        print("Combining MRC files into a Numpy array...")

    # get stem of pdb file
    pdb_stem = str(pathlib.Path(pdb_path).stem)
    tmp_dir = f"./tmp_data/{pdb_stem}/"
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

    ss_dir = os.path.join(tmp_dir, "ss")
    pdb_dir = os.path.join(tmp_dir, "pdb")
    simu_mrc_dir = os.path.join(tmp_dir, "simu_mrc")

    os.makedirs(ss_dir, exist_ok=True)
    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs(simu_mrc_dir, exist_ok=True)

    out_ss = assign_ss(pdb_path, ss_dir, verbose)
    split_pdb_by_ss(pdb_path, out_ss, pdb_dir, verbose)

    pdb_simu_map_path = os.path.join(tmp_dir, f"{pdb_stem}_simu_map.mrc")

    gen_simu_map(pdb_path, sample_res, pdb_simu_map_path)

    pdb_simu_map = MapParser.readMRC(pdb_simu_map_path)

    for file in os.listdir(pdb_dir):
        filename = str(pathlib.Path(file).stem)
        if (
            file.endswith("ssA.pdb")
            or file.endswith("ssB.pdb")
            or file.endswith("ssC.pdb")
        ):
            gen_simu_map(
                os.path.join(pdb_dir, file),
                sample_res,
                os.path.join(simu_mrc_dir, filename + ".mrc"),
                densMap=pdb_simu_map,
            )

    with mrcfile.open(pdb_simu_map_path) as mrc:
        dims = mrc.data.shape

        arr = np.zeros((4, dims[0], dims[1], dims[2]))

        for file in os.listdir(simu_mrc_dir):
            if file.endswith("ssC.pdb.mrc"):
                arr[0] = mrcfile.open(os.path.join(simu_mrc_dir, file)).data.copy()
            elif file.endswith("ssB.pdb.mrc"):
                arr[1] = mrcfile.open(os.path.join(simu_mrc_dir, file)).data.copy()
            elif file.endswith("ssA.pdb.mrc"):
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
        print(arr.shape)
        # print stats
        print("Number in SS class Coil: ", np.count_nonzero(arr[..., 0]))
        print("Number in SS class Beta: ", np.count_nonzero(arr[..., 1]))
        print("Number in SS class Alpha: ", np.count_nonzero(arr[..., 2]))

        # print min, max, mean, std
        print(
            "Coil: ",
            np.min(arr[..., 0]),
            np.max(arr[..., 0]),
            np.mean(arr[..., 0]),
            np.std(arr[..., 0]),
        )
        print(
            "Beta: ",
            np.min(arr[..., 1]),
            np.max(arr[..., 1]),
            np.mean(arr[..., 1]),
            np.std(arr[..., 1]),
        )
        print(
            "Alpha: ",
            np.min(arr[..., 2]),
            np.max(arr[..., 2]),
            np.mean(arr[..., 2]),
            np.std(arr[..., 2]),
        )

    return arr
