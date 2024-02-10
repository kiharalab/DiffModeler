# pdb2vol.py adopted from BioTEMPy

from Bio.PDB import PDBParser, MMCIFParser
import numpy as np

import mrcfile
from scipy.ndimage import fourier_gaussian, gaussian_filter, zoom
from scipy.fftpack import fftn, ifftn

# Dictionary of atom types and their corresponding masses.
atom_mass_dict = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,  # for DNA/RNA
    "S": 32.066,
}


def get_atom_list(pdb_file):
    """
    Retrieve the coordinates and atom types from a PDB or CIF file.

    Parameters:
    pdb_file (str): The path to the PDB or CIF file.

    Returns:
    tuple: A tuple containing two elements:
        - np.array: An array of atom coordinates.
        - list: A list of atom types.
    """
    if pdb_file.endswith(".pdb"):
        st_parser = PDBParser(QUIET=True)
    elif pdb_file.endswith(".cif"):
        st_parser = MMCIFParser(QUIET=True)
    structure = st_parser.get_structure("protein", pdb_file)
    atom_list = []
    atom_type_list = []
    for atom in structure.get_atoms():
        atom_list.append(atom.get_coord())
        atom_type_list.append(atom.element)
    return np.array(atom_list), atom_type_list


def calculate_centre_of_mass(atom_list, atom_type_list):
    """
    Calculates the centre of mass for a given list of atoms and their types.

    Args:
        atom_list (np.ndarray): List of atom coordinates in the form (x, y, z).
        atom_type_list (list): List of atom types.

    Returns:
        tuple: The coordinates of the centre of mass in the form (x_co_m, y_co_m, z_co_m).
    """
    atom_list = np.array(atom_list)
    atom_type_list = np.array(atom_type_list)
    x = atom_list[:, 0]
    y = atom_list[:, 1]
    z = atom_list[:, 2]
    m = np.array([atom_mass_dict.get(atom_type, 0.0) for atom_type in atom_type_list])
    mass_total = np.sum(m)
    x_co_m = np.sum(x * m) / mass_total
    y_co_m = np.sum(y * m) / mass_total
    z_co_m = np.sum(z * m) / mass_total
    return x_co_m, y_co_m, z_co_m


def prot2map(
    atom_list,
    atom_type_list,
    voxel_size,
    resolution=None,
):
    """
    Calculate the size and origin of a protein map based on the given atom list,
    atom type list, voxel size, and optional resolution.

    Args:
        atom_list (numpy.ndarray): Array of atom coordinates.
        atom_type_list (list): Array of atom types.
        voxel_size (numpy.ndarray): The size of each voxel in x, y, and z dimensions.
        resolution (float, optional): Resolution of the map. Defaults to None.

    Returns:
        tuple: A tuple containing the size of the map in each dimension (z, y, x)
               and the origin of the map (x_origin, y_origin, z_origin).
    """
    max_x, max_y, max_z = atom_list.max(axis=0)
    min_x, min_y, min_z = atom_list.min(axis=0)

    if resolution is not None:
        edge = np.array(2 * resolution / voxel_size, dtype=int) + 4
    else:
        edge = np.array([10, 10, 10], dtype=int)

    x_size = int((max_x - min_x) / voxel_size[0]) + edge[0]
    y_size = int((max_y - min_y) / voxel_size[1]) + edge[1]
    z_size = int((max_z - min_z) / voxel_size[2]) + edge[2]

    CoM = calculate_centre_of_mass(atom_list, atom_type_list)

    # Origin calculated such that the centre of the map is the centre of
    # mass of the protein.
    half_x = max(CoM[0] - min_x, max_x - CoM[0])

    if half_x < (voxel_size[0] * x_size / 2.0):
        half_x = voxel_size[0] * x_size / 2.0
    x_origin = CoM[0] - half_x - edge[0] * voxel_size[0]
    x_size = int(half_x * 2.0 / voxel_size[0] + 2 * edge[0])
    half_y = max(CoM[1] - min_y, max_y - CoM[1])

    if half_y < (voxel_size[1] * y_size / 2.0):
        half_y = voxel_size[1] * y_size / 2.0
    y_origin = CoM[1] - half_y - edge[1] * voxel_size[1]
    y_size = int(half_y * 2.0 / voxel_size[1] + 2 * edge[1])
    half_z = max(CoM[2] - min_z, max_z - CoM[2])

    if half_z < (voxel_size[2] * z_size / 2.0):
        half_z = voxel_size[2] * z_size / 2.0
    z_origin = CoM[2] - half_z - edge[2] * voxel_size[2]
    z_size = int(half_z * 2.0 / voxel_size[2] + 2 * edge[2])

    return (z_size, y_size, x_size), (x_origin, y_origin, z_origin)


def mapGridPosition(origin, voxel_size, box_size, atom_coord):
    """
    Maps the coordinates of an atom to the corresponding grid position in a voxel grid.

    Parameters:
    origin (tuple): The origin of the voxel grid.
    voxel_size (tuple): The size of each voxel in the grid.
    box_size (tuple): The size of the voxel grid.
    atom_coord (tuple): The coordinates of the atom.

    Returns:
    tuple: The grid position (x, y, z) of the atom in the voxel grid.

    If the atom is outside the voxel grid, returns (0, 0, 0).
    """
    x_pos = int(round((atom_coord[0] - origin[0]) / voxel_size[0], 0))
    y_pos = int(round((atom_coord[1] - origin[1]) / voxel_size[1], 0))
    z_pos = int(round((atom_coord[2] - origin[2]) / voxel_size[2], 0))

    if (box_size[2] > x_pos >= 0) and (box_size[1] > y_pos >= 0) and (box_size[0] > z_pos >= 0):
        return x_pos, y_pos, z_pos
    else:
        return 0


def make_atom_overlay_map(origin, voxel_size, box_size, atom_list, atom_type_list):
    """
    Creates an atom overlay map based on the given parameters.

    Parameters:
    origin (tuple): The origin coordinates of the map.
    voxel_size (float): The size of each voxel in the map.
    box_size (tuple): The size of the map in each dimension.
    atom_list (list): The list of atom coordinates.
    atom_type_list (list): The list of atom types.

    Returns:
    numpy.ndarray: The atom overlay map.
    """
    map_data = np.zeros(box_size)
    for atom, atom_type in zip(atom_list, atom_type_list):
        pos = mapGridPosition(origin, voxel_size, box_size, atom)
        if pos:
            map_data[pos[2], pos[1], pos[0]] += atom_mass_dict.get(atom_type, 0.0)
    return map_data


def write_mrc_file(data, origin, voxel_size, mrc_file):
    """
    Write a data array to an MRC file.

    Args:
        data (ndarray): The data array to be written.
        origin (tuple): The origin coordinates of the data.
        voxel_size (array_like): The voxel size of the data.
        mrc_file (str): The path to the output MRC file.

    Returns:
        None
    """
    with mrcfile.new(mrc_file, overwrite=True) as mrc:
        mrc.set_data(data.astype(np.float32))
        mrc.update_header_from_data()
        mrc.voxel_size = tuple(voxel_size)
        mrc.header.origin.x = origin[0]
        mrc.header.origin.y = origin[1]
        mrc.header.origin.z = origin[2]
        mrc.update_header_stats()
        mrc.flush()


def blur_map(data, resolution, sigma_coeff):
    """
    Blurs the input data using a Gaussian filter.

    Args:
        data (ndarray): The input data to be blurred.
        resolution (float): The resolution of the data.
        sigma_coeff (float): The coefficient to determine the sigma value for the Gaussian filter.

    Returns:
        ndarray: The blurred data.

    """
    sigma = resolution * sigma_coeff
    new_data = fourier_gaussian(fftn(data), sigma)
    return np.real(ifftn(new_data))


def blur_map_real_space(data, resolution, sigma_coeff):
    """
    Blurs a map in real space using a Gaussian filter.

    Args:
        data (ndarray): The input map data.
        resolution (float): The resolution of the map.
        sigma_coeff (float): The coefficient to determine the sigma value for the Gaussian filter.

    Returns:
        ndarray: The blurred map data.
    """
    sigma = resolution * sigma_coeff
    print(sigma)
    new_data = gaussian_filter(data, sigma)
    return new_data


def normalize_map(map_data):
    """
    Normalize a map by subtracting the mean and dividing by the standard deviation.

    Parameters:
    map_data (numpy.ndarray): The input map data.

    Returns:
    numpy.ndarray: The normalized map data.
    """
    if map_data.std() != 0:
        return (map_data - map_data.mean()) / map_data.std()
    else:
        return map_data


def resample_by_box_size(data, box_size):
    """
    Resamples the given data array to match the specified box size using spline interpolation.

    Parameters:
        data (ndarray): The input data array.
        box_size (tuple): The desired box size for resampling.

    Returns:
        ndarray: The resampled data array.
    """
    # cubic spline interpolation
    zoom_factor = np.array(box_size) / np.array(data.shape)
    return zoom(data, zoom_factor, order=3)


def pdb2vol(input_pdb, output_mrc, resolution, ref_map=False, sigma_coeff=0.356, real_space=False, normalize=True):
    """
    Convert a PDB or CIF file to a volumetric map in MRC format.

    Args:
        input_pdb (str): Path to the input PDB or CIF file.
        output_mrc (str): Path to save the output MRC file.
        resolution (float): Resolution of the output map.
        ref_map (str, optional): Path to a reference map in MRC format. Defaults to False.
        sigma_coeff (float, optional): Sigma coefficient for blurring. Defaults to 0.356.
        real_space (bool, optional): Whether to perform real-space blurring. Defaults to False.
        normalize (bool, optional): Whether to normalize the output map. Defaults to True.

    Raises:
        ValueError: If the input file is not a PDB or CIF file.
        ValueError: If no atoms are found in the input file.
        ValueError: If the number of atoms and atom types do not match.

    Returns:
        None
    """

    # sigma_coeff = 1/(pi*sqrt(2*log(2)) = 0.187, makes the FT fall to half maximum at wavenumber 1/resolution
    # sigma_coeff = 1/(pi*sqrt(2)) = 0.225 is the default value used in Chimera, makes the Fourier transform (FT) of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution
    # sigma_coeff = 1/(2*sqrt(2)) =  0.356 makes the Gaussian width at 1/e maximum height equal the resolution
    # sigma_coeff = 1/(2*sqrt(2*log(2))) = 0.4247 makes the Gaussian width at half maximum height equal the resolution

    if input_pdb.split(".")[-1] not in ["pdb", "cif"]:
        raise ValueError("Input file must be a pdb or cif file")
    atoms, types = get_atom_list(input_pdb)

    if len(atoms) == 0:
        raise ValueError("No atoms found in input file")
    if len(atoms) != len(types):
        raise ValueError("Number of atoms and atom types does not match")

    if not ref_map:
        r = np.clip(resolution / 4.0, a_min=1.0, a_max=3.5)
        voxel_size = np.array([r, r, r])
        dims, origin = prot2map(atoms, types, voxel_size, resolution)
    else:
        with mrcfile.open(ref_map, permissive=True) as mrc:
            voxel_size = np.array([mrc.voxel_size.x, mrc.voxel_size.y, mrc.voxel_size.z])
            dims = mrc.data.shape
            origin = np.array([mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z])

    x_s = int(dims[2] * voxel_size[2])
    y_s = int(dims[1] * voxel_size[1])
    z_s = int(dims[0] * voxel_size[0])

    new_voxel_size = np.array([voxel_size[2] * dims[2] / x_s, voxel_size[1] * dims[1] / y_s, voxel_size[0] * dims[0] / z_s])

    map_data = make_atom_overlay_map(origin, new_voxel_size, np.array((z_s, y_s, x_s)), atoms, types)
    if real_space:
        blurred_data = blur_map_real_space(map_data, resolution, sigma_coeff)
    else:
        blurred_data = blur_map(map_data, resolution, sigma_coeff)
    blurred_data = resample_by_box_size(blurred_data, dims)
    if normalize:
        blurred_data = normalize_map(blurred_data)
    write_mrc_file(blurred_data, origin, voxel_size, output_mrc)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert a PDB or CIF file to a volumetric map in MRC format.')
    parser.add_argument('input_pdb', help='Path to the input PDB or CIF file.')
    parser.add_argument('output_mrc', help='Path to save the output MRC file.')
    parser.add_argument('resolution', type=float, help='Resolution of the output map.')
    parser.add_argument('-m', '--ref_map', help='Path to a reference map in MRC format.')
    parser.add_argument('-s', '--sigma_coeff', type=float, default=0.356, help='Sigma coefficient for blurring.')
    parser.add_argument('-r', '--real_space', action='store_true', default=False, help='Whether to perform real-space blurring.')
    parser.add_argument('-n', '--normalize', action='store_true', default=True, help='Whether to normalize the output map.')
    args = parser.parse_args()

    pdb2vol(args.input_pdb, args.output_mrc, args.resolution, args.ref_map, args.sigma_coeff, args.real_space, args.normalize)
