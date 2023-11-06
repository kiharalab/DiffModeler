import copy

import mrcfile
import numpy as np
from numba import jit

# from julia.api import Julia
# jl = Julia(compiled_modules=False)
# from julia import Main


class EMmap:
    """A mrc object that represents the density data and statistics of a given mrc file"""

    def __init__(self, path, ss_data=None):
        """
        The EMmap class takes in a path to a mrc file and reads the data and header information.
        ss_data is an optional parameter that takes in a numpy array of secondary structure data information.
        """
        # open the specified mrcfile and read the header information

        self.mrcfile_path = path
        self.ss_data = ss_data
        self.new_ss_data = None
        if self.ss_data is not None:
            self.ss_data = np.swapaxes(self.ss_data, 0, 2)  # order is coil, beta, alpha, nucleotides

        # read the mrc file
        mrc = mrcfile.open(path, permissive=True)
        data = mrc.data
        header = mrc.header

        # read and store the voxel widths and dimensions from the header
        self.xdim = int(header.nx)
        self.ydim = int(header.ny)
        self.zdim = int(header.nz)

        self.xwidth = mrc.voxel_size.x
        self.ywidth = mrc.voxel_size.y
        self.zwidth = mrc.voxel_size.z

        # set the center to be the half the dimensions
        self.cent = np.array(
            [
                self.xdim * 0.5,
                self.ydim * 0.5,
                self.zdim * 0.5,
            ]
        )

        # read and store the origin coordinate from the header
        self.orig = np.array((header.origin.x, header.origin.y, header.origin.z))

        if np.all(self.orig == 0):
            # MRC2000 format uses nxstart, nystart, nzstart instead of origin
            self.orig_idx = np.array((header.nxstart, header.nystart, header.nzstart))
            self.orig = self.orig_idx * np.array((self.xwidth, self.ywidth, self.zwidth))

        # swap the xz axes of density data array and store in self.data
        # also convert the data type to float32
        self.data = np.swapaxes(copy.deepcopy(data), 0, 2).astype(np.float32)

        # initialize the vector array to be same shape as data but will all zeros
        self.vec = np.zeros((self.xdim, self.ydim, self.zdim, 3), dtype="float32")

        # initialize all the statistics values
        self.dsum = None  # total density value
        self.Nact = None  # non-zero density voxel count
        self.ave = None  # average density value
        self.std_norm_ave = None  # L2 norm normalized with average density value
        self.std = None  # denormalize L2 norm

        # resampled map parameter and fitter related
        self.new_dim = None
        self.new_width = None
        self.new_cent = None
        self.new_orig = None
        self.new_data = None

    def set_vox_size(self, thr=0.0, voxel_size=7.0):
        """Set the voxel size according to the given threshold and granularity

        Args:
            thr (float, optional): preset threshold for density cutoff. Defaults to 0.
            voxel_size (float, optional): the granularity of the voxel in terms of angstroms. Defaults to 7.

        """

        # if th < 0 add th to all value
        if thr < 0:
            self.data = self.data - thr
            thr = 0.0

        # zero all the values less than threshold
        self.data[self.data < thr] = 0.0

        # calculate maximum distance for non-zero entries
        non_zero_index_list = np.array(np.nonzero(self.data)).T
        if len(non_zero_index_list) == 0:
            dmax = 0
        else:
            d2_list = np.linalg.norm(non_zero_index_list - np.array(self.cent), axis=1)
            dmax = np.max(d2_list)

        print()
        print("#dmax=" + str(dmax / self.xwidth))

        # set new center
        self.new_cent = self.cent * self.xwidth + self.orig

        tmp_size = 2 * dmax * self.xwidth / voxel_size

        # get the best size suitable for fft operation
        from pyfftw import pyfftw

        new_dim = pyfftw.next_fast_len(int(tmp_size))

        # set new origins
        self.new_orig = self.new_cent - 0.5 * new_dim * voxel_size
        self.new_dim = new_dim
        self.new_width = voxel_size

        print("Nvox= " + str(self.new_dim) + ", " + str(self.new_dim) + ", " + str(self.new_dim))
        print("cent= " + str(self.new_cent[0]) + ", " + str(self.new_cent[1]) + ", " + str(self.new_cent[2]))
        print("ori= " + str(self.new_orig[0]) + ", " + str(self.new_orig[1]) + ", " + str(self.new_orig[2]))

    def resample_and_vec(self, dreso=16.0, density_map=None):
        """
        The resample_and_vec function takes a map and resamples it to the desired resolution.
        It also calculates the vectors needed for calculating the DOT score.
        The input parameters are:

        :param dreso: Set the resolution of the map
        :param density_map: Pass in a density map that is used to mask the data
        :return: The resampled data, the vectors and map statistics
        """
        src_dims = np.array((self.xdim, self.ydim, self.zdim))

        res_data, res_vec, res_ss_data = do_resample_and_vec(
            self.xwidth,
            self.orig,
            src_dims,
            self.data,
            self.new_width,
            self.new_orig,
            self.new_dim,
            dreso,
            density_map,
            self.ss_data,
        )

        # Main.include("resample_and_vec.jl")
        # res_data, res_vec = Main.res_vec_jl(float(self.xwidth), self.orig, src_dims, self.data, self.new_width, self.new_orig, int(self.new_dim), dreso)

        # calculate map statistics
        density_sum = np.sum(res_data)
        count = np.count_nonzero(res_data)
        ave = np.mean(res_data[res_data > 0])
        std = np.linalg.norm(res_data[res_data > 0])
        std_norm_ave = np.linalg.norm(res_data[res_data > 0] - ave)

        print(
            "#MAP SUM={sum} COUNT={cnt} AVE={ave} STD={std} STD_norm={std_norm}".format(
                sum=density_sum, cnt=count, ave=ave, std=std, std_norm=std_norm_ave
            )
        )
        if self.ss_data is not None:
            print(
                f"#SS Coil SUM={np.sum(res_ss_data[..., 0])}, "
                f"#COUNT={np.count_nonzero(res_ss_data[..., 0])}, "
                f"#AVE={np.mean(res_ss_data[res_ss_data[..., 0] > 0])}, "
                f"#STD={np.linalg.norm(res_ss_data[res_ss_data[..., 0] > 0])}, "
                f"#STD_norm={np.linalg.norm(res_ss_data[res_ss_data[..., 0] > 0] - ave)}"
            )
            print(
                f"#SS Beta SUM={np.sum(res_ss_data[..., 1])}, "
                f"#COUNT={np.count_nonzero(res_ss_data[..., 1])}, "
                f"#AVE={np.mean(res_ss_data[res_ss_data[..., 1] > 0])}, "
                f"#STD={np.linalg.norm(res_ss_data[res_ss_data[..., 1] > 0])}, "
                f"#STD_norm={np.linalg.norm(res_ss_data[res_ss_data[..., 1] > 0] - ave)}"
            )
            print(
                f"#SS Alpha Sum={np.sum(res_ss_data[..., 2])}, "
                f"#COUNT={np.count_nonzero(res_ss_data[..., 2])}, "
                f"#AVE={np.mean(res_ss_data[res_ss_data[..., 2] > 0])}, "
                f"#STD={np.linalg.norm(res_ss_data[res_ss_data[..., 2] > 0])}, "
                f"#STD_norm={np.linalg.norm(res_ss_data[res_ss_data[..., 2] > 0] - ave)}"
            )
        else:
            res_ss_data = None
        # update the dest object with the new data and vectors
        self.new_data = res_data
        self.vec = res_vec
        self.new_ss_data = res_ss_data
        self.dsum = density_sum
        self.Nact = count
        self.ave = ave
        self.std = std
        self.std_norm_ave = std_norm_ave


def unify_dims(map_list, voxel_size):
    """
    The unify_dims function takes a list of EM maps and a voxel size as input.
    It then finds the maximum dimension among all the maps in the list, and sets
    the new_dim attribute of each map to this value. It also adjusts each map's
    new_orig attribute so that it is centered on its original center point.

    :param map_list: Pass in the list of maps to be unified
    :param voxel_size: Calculate the new origin of the map
    :return: A list of maps with the same dimensions
    """
    dims = np.array([em_map.new_dim for em_map in map_list])
    max_dim = np.max(dims)
    for em_map in map_list:
        em_map.new_dim = max_dim
        if em_map.xdim != max_dim:
            em_map.new_orig = em_map.new_cent - 0.5 * voxel_size * max_dim


@jit(nopython=True)
def calc_prob(stp, endp, pos, density_data, prob_data, fsiv):
    """Density weighted mean shift algorithm using Gaussian filter to sample data in original MRC density map"""
    dtotal = 0.0
    pos2 = np.zeros((3,))

    for xp in range(stp[0], endp[0]):
        rx = float(xp) - pos[0]
        rx = rx**2
        for yp in range(stp[1], endp[1]):
            ry = float(yp) - pos[1]
            ry = ry**2
            for zp in range(stp[2], endp[2]):
                rz = float(zp) - pos[2]
                rz = rz**2
                d2 = rx + ry + rz
                # v = density_data[xp][yp][zp] * prob_data[xp][yp][zp] *  np.exp(-1.5 * d2 * fsiv)
                v = prob_data[xp][yp][zp] * np.exp(-1.5 * d2 * fsiv)
                dtotal += v
                pos2[0] += v * xp
                pos2[1] += v * yp
                pos2[2] += v * zp

    return dtotal, pos2


@jit(nopython=True)
def calc(stp, endp, pos, data, fsiv):
    """Mean shift algorithm using Gaussian filter to sample data in original MRC density map"""
    dtotal = 0.0
    pos2 = np.zeros((3,))

    for xp in range(stp[0], endp[0]):
        rx = float(xp) - pos[0]
        rx = rx**2
        for yp in range(stp[1], endp[1]):
            ry = float(yp) - pos[1]
            ry = ry**2
            for zp in range(stp[2], endp[2]):
                rz = float(zp) - pos[2]
                rz = rz**2
                d2 = rx + ry + rz
                v = data[xp][yp][zp] * np.exp(-1.5 * d2 * fsiv)
                dtotal += v
                pos2[0] += v * xp
                pos2[1] += v * yp
                pos2[2] += v * zp

    return dtotal, pos2


@jit(nopython=True)
def do_resample_and_vec(
    src_xwidth,
    src_orig,
    src_dims,
    src_data,
    dest_xwidth,
    dest_orig,
    new_dim,
    dreso,
    density_map,
    ss_data,
):
    """
    The do_resample_and_vec function takes in the following parameters:
        src_xwidth - The width of the source map.
        src_orig - The origin of the source map.
        src_dims - The dimensions of the source map.
        src_data - A array contains the density data. The shape of the array is (xdim, ydim, zdim).

    :param src_xwidth: Calculate the grid step
    :param src_orig: Define the origin of the source volume
    :param src_dims: Determine the size of the source volume
    :param src_data: Data of the source volume
    :param dest_xwidth: Determine the size of the new grid
    :param dest_orig: Define the origin of the destination volume
    :param new_dim: Determine the size of the new volume
    :param dreso: Determine the Guassian filter size
    :param ss_data: data of the secondary structure
    :param : Determine the resolution of the new map
    :return: The density map, the unit vectors and the secondary structure data
    """

    gstep = src_xwidth
    fs = (dreso / gstep) * 0.5
    fs = fs**2
    fsiv = 1.0 / fs
    fmaxd = (dreso / gstep) * 2.0
    print("#maxd=", fmaxd)
    print("#fsiv=", fsiv)

    dest_vec = np.zeros((new_dim, new_dim, new_dim, 3), dtype="float32")
    dest_data = np.zeros((new_dim, new_dim, new_dim), dtype="float32")
    dest_ss_data = np.zeros((new_dim, new_dim, new_dim, 4), dtype="float32")

    for x in range(new_dim):
        for y in range(new_dim):
            for z in range(new_dim):

                xyz_arr = np.array((x, y, z))
                pos = (xyz_arr * dest_xwidth + dest_orig - src_orig) / src_xwidth

                # check density

                if (
                    pos[0] < 0
                    or pos[1] < 0
                    or pos[2] < 0
                    or pos[0] >= src_dims[0]
                    or pos[1] >= src_dims[1]
                    or pos[2] >= src_dims[2]
                ):
                    continue

                if src_data[int(pos[0])][int(pos[1])][int(pos[2])] == 0:
                    continue

                # Start Point
                stp = (pos - fmaxd).astype(np.int32)

                # set start and end point
                if stp[0] < 0:
                    stp[0] = 0
                if stp[1] < 0:
                    stp[1] = 0
                if stp[2] < 0:
                    stp[2] = 0

                # End Point
                endp = (pos + fmaxd + 1).astype(np.int32)

                if endp[0] >= src_dims[0]:
                    endp[0] = src_dims[0]
                if endp[1] >= src_dims[1]:
                    endp[1] = src_dims[1]
                if endp[2] >= src_dims[2]:
                    endp[2] = src_dims[2]

                # compute the total density
                if density_map is not None:
                    dtotal, pos2 = calc_prob(stp, endp, pos, density_map, src_data, fsiv)
                else:
                    dtotal, pos2 = calc(stp, endp, pos, src_data, fsiv)

                if ss_data is not None:
                    dest_ss_data[x][y][z][0], _ = calc(stp, endp, pos, ss_data[..., 0], fsiv)
                    dest_ss_data[x][y][z][1], _ = calc(stp, endp, pos, ss_data[..., 1], fsiv)
                    dest_ss_data[x][y][z][2], _ = calc(stp, endp, pos, ss_data[..., 2], fsiv)
                    dest_ss_data[x][y][z][3], _ = calc(stp, endp, pos, ss_data[..., 3], fsiv)

                if dtotal == 0:
                    continue

                dest_data[x][y][z] = dtotal

                rd = 1.0 / dtotal

                pos2 *= rd

                tmpcd = pos2 - pos

                dvec = np.sqrt(tmpcd[0] ** 2 + tmpcd[1] ** 2 + tmpcd[2] ** 2)

                if dvec == 0:
                    dvec = 1.0

                rdvec = 1.0 / dvec

                dest_vec[x][y][z] = tmpcd * rdvec

    return dest_data, dest_vec, dest_ss_data
