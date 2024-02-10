import copy

import numpy as np


def rot_mrc(orig_mrc_data, orig_mrc_vec, mtx):
    """A function to rotation the density and vector array by a specified angle.
    Args:
        orig_mrc_data (numpy.array): the data array to be rotated
        orig_mrc_vec (numpy.array): the vector array to be rotated
        mtx (numpy.array): the rotation matrix to be used
    Returns:
        new_vec_array (numpy.array): rotated vector array
        new_data_array (numpy.array): rotated data array
    """

    # set the dimension to be x dimension as all dimension are the same
    dim = orig_mrc_data.shape[0]

    # create array for the positions after rotation
    new_pos = np.array(
        np.meshgrid(
            np.arange(dim),
            np.arange(dim),
            np.arange(dim),
        )
    ).T.reshape(-1, 3)

    # set the rotation center
    cent = 0.5 * float(dim)

    # get relative new positions from center
    new_pos = new_pos - cent

    # reversely rotate the new position lists to get old positions
    old_pos = np.einsum("ij, kj->ki", mtx.T, new_pos) + cent

    # concatenate combine two position array horizontally for later filtering
    combined_arr = np.hstack((old_pos, new_pos))

    # filter out the positions that are out of the original array
    in_bound_mask = (
        (combined_arr[:, 0] >= 0)
        * (combined_arr[:, 1] >= 0)
        * (combined_arr[:, 2] >= 0)
        * (combined_arr[:, 0] < dim)
        * (combined_arr[:, 1] < dim)
        * (combined_arr[:, 2] < dim)
    )

    # init new vec and dens array
    new_vec_array = np.zeros_like(orig_mrc_vec)
    new_data_array = np.zeros_like(orig_mrc_data)

    # get the mask of all the values inside boundary
    combined_arr = combined_arr[in_bound_mask]

    # convert the index to integer
    combined_arr = combined_arr.astype(np.int32)

    # get the old index array
    index_arr = combined_arr[:, 0:3]

    # index_arr = old_pos[in_bound_mask].astype(np.int32)

    # get the index that has non-zero density by masking

    dens_mask = orig_mrc_data[index_arr[:, 0], index_arr[:, 1], index_arr[:, 2]] != 0.0
    non_zero_rot_list = combined_arr[dens_mask]

    # dens_mask = orig_mrc_data[index_arr[:, 0], index_arr[:, 1], index_arr[:, 2]] != 0.0
    # non_zero_rot_list = old_pos[dens_mask].astype(np.int32)

    # get the non-zero vec and dens values
    non_zero_vec = orig_mrc_vec[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens = orig_mrc_data[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]

    # rotate the vectors
    new_vec = np.einsum("ij, kj->ki", mtx, non_zero_vec)

    # find the new indices
    new_ind_arr = (non_zero_rot_list[:, 3:6] + cent).astype(np.int32)
    # new_ind_arr = (new_pos[dens_mask] + cent).astype(np.int32)

    # fill in the values to new vec and dens array
    new_vec_array[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = new_vec
    new_data_array[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens

    return new_vec_array, new_data_array


def rot_mrc_prob(data, vec, prob_c1, prob_c2, prob_c3, prob_c4, mtx):
    dim = data.shape[0]

    new_pos = np.array(
        np.meshgrid(
            np.arange(dim),
            np.arange(dim),
            np.arange(dim),
        )
    ).T.reshape(-1, 3)

    cent = 0.5 * float(dim)
    new_pos = new_pos - cent

    old_pos = np.einsum("ij, kj->ki", mtx.T, new_pos) + cent

    combined_arr = np.hstack((old_pos, new_pos))

    in_bound_mask = (
        (old_pos[:, 0] >= 0)
        * (old_pos[:, 1] >= 0)
        * (old_pos[:, 2] >= 0)
        * (old_pos[:, 0] < dim)
        * (old_pos[:, 1] < dim)
        * (old_pos[:, 2] < dim)
    )

    # create new array for density, vector and probability
    new_vec_array = np.zeros_like(vec)
    new_data_array = np.zeros_like(data)
    new_p1 = np.zeros_like(prob_c1)
    new_p2 = np.zeros_like(prob_c2)
    new_p3 = np.zeros_like(prob_c3)
    new_p4 = np.zeros_like(prob_c4)

    combined_arr = combined_arr[in_bound_mask]

    combined_arr = combined_arr.astype(np.int32)

    index_arr = combined_arr[:, 0:3]

    dens_mask = data[index_arr[:, 0], index_arr[:, 1], index_arr[:, 2]] != 0.0

    non_zero_rot_list = combined_arr[dens_mask]

    # get the index of the non-zero density
    non_zero_vec = vec[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens = data[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens_p1 = prob_c1[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens_p2 = prob_c2[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens_p3 = prob_c3[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]
    non_zero_dens_p4 = prob_c4[non_zero_rot_list[:, 0], non_zero_rot_list[:, 1], non_zero_rot_list[:, 2]]

    # find the new indices
    new_ind_arr = (non_zero_rot_list[:, 3:6] + cent).astype(int)

    # save the rotated data
    new_vec_array[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = np.einsum("ij, kj->ki", mtx, non_zero_vec)
    new_data_array[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens
    new_p1[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens_p1
    new_p2[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens_p2
    new_p3[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens_p3
    new_p4[new_ind_arr[:, 0], new_ind_arr[:, 1], new_ind_arr[:, 2]] = non_zero_dens_p4

    return new_vec_array, new_data_array, new_p1, new_p2, new_p3, new_p4


# @numba.jit(nopython=True)
# def calc(stp, endp, pos, mrc1_data, fsiv):
#     """Vectorized version of calc"""
#
#     xx = np.arange(stp[0], endp[0], 1)
#     yy = np.arange(stp[1], endp[1], 1)
#     zz = np.arange(stp[2], endp[2], 1)
#
#     xx = np.expand_dims(xx, axis=1)
#     xx = np.expand_dims(xx, axis=1)
#     yy = np.expand_dims(yy, axis=1)
#     yy = np.expand_dims(yy, axis=0)
#     zz = np.expand_dims(zz, axis=0)
#     zz = np.expand_dims(zz, axis=0)
#
#     # calculate the distance between the center of the voxel and the center of the particle
#     d2 = (xx - pos[0])**2 + (yy - pos[1])**2 + (zz - pos[2])**2
#
#     # calculate the density and vector in resized map using Gaussian interpolation in original MRC density map
#     d = np.exp(- 1.5 * d2 * fsiv) * mrc1_data[stp[0]:endp[0], stp[1]:endp[1], stp[2]:endp[2]]
#     dtotal = np.sum(d)
#
#     # calculate the vector
#     v = np.array([np.sum(d * xx), np.sum(d * yy), np.sum(d * zz)])
#
#     return dtotal, v

# def rot_mrc(orig_mrc_data, orig_mrc_vec, mtx, interp=interp):
#     """A function to rotation the density and vector array by a specified angle.

#     Args:
#         orig_mrc_data (numpy.array): the data array to be rotated
#         orig_mrc_vec (numpy.array): the vector array to be rotated
#         mtx (numpy.array): the rotation matrix
#         interp (str): the interpolation method

#     Returns:
#         new_vec_array (numpy.array): rotated vector array
#         new_data_array (numpy.array): rotated data array
#     """

#     Nx, Ny, Nz = orig_mrc_data.shape
#     x = np.linspace(0, Nx - 1, Nx)
#     y = np.linspace(0, Ny - 1, Ny)
#     z = np.linspace(0, Nz - 1, Nz)
#     # xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
#     zz, yy, xx = np.meshgrid(z, y, x, indexing='ij')

#     x_center = x[0] + x[-1] / 2
#     y_center = y[0] + y[-1] / 2
#     z_center = z[0] + z[-1] / 2

#     # center the coord
#     coor = np.array([xx - x_center, yy - y_center, zz - z_center])
#     # apply rotation
#     # coor_prime = np.tensordot(np.flip(mtx.T), coor, axes=((0), (1)))
#     coor_prime = np.einsum("il, ljkm->ijkm", mtx.T, coor)

#     # uncenter the coord
#     xx_prime = coor_prime[0] + x_center
#     yy_prime = coor_prime[1] + y_center
#     zz_prime = coor_prime[2] + z_center

#     # trim the values outside boundaries
#     x_valid1 = xx_prime >= 0
#     x_valid2 = xx_prime <= Nx - 1
#     y_valid1 = yy_prime >= 0
#     y_valid2 = yy_prime <= Ny - 1
#     z_valid1 = zz_prime >= 0
#     z_valid2 = zz_prime <= Nz - 1

#     # get non-zero indicies in original density
#     nonzero_dens = orig_mrc_data > 0

#     # get voxels with all valid dimensions
#     valid_voxel = x_valid1 * x_valid2 * y_valid1 * y_valid2 * z_valid1 * z_valid2 * nonzero_dens

#     # get nonzero positions
#     #x_valid_idx, y_valid_idx, z_valid_idx = np.where(valid_voxel > 0)
#     z_valid_idx, y_valid_idx, x_valid_idx = np.where(valid_voxel > 0)

#     # create new arrays to store the final result
#     new_data_array = np.zeros_like(orig_mrc_data)
#     new_vec_array = np.zeros_like(orig_mrc_vec)

#     # gather points to be interpolated
#     # interp_points = np.array(
#     #     [
#     #         xx_prime[x_valid_idx, y_valid_idx, z_valid_idx],
#     #         yy_prime[x_valid_idx, y_valid_idx, z_valid_idx],
#     #         zz_prime[x_valid_idx, y_valid_idx, z_valid_idx],
#     #     ]
#     # ).T

#     interp_points = np.array(
#         [
#             zz_prime[z_valid_idx, y_valid_idx, x_valid_idx],
#             yy_prime[z_valid_idx, y_valid_idx, x_valid_idx],
#             xx_prime[z_valid_idx, y_valid_idx, x_valid_idx],
#         ]
#     ).T

#     if interp is not None:
#         # create grid interpolator
#         # data_w_coor = RegularGridInterpolator((x, y, z), orig_mrc_data, method=interp)
#         # vec_w_coor = RegularGridInterpolator((x, y, z), orig_mrc_vec, method=interp)

#         data_w_coor = RegularGridInterpolator((z, y, x), orig_mrc_data, method=interp)
#         vec_w_coor = RegularGridInterpolator((z, y, x), orig_mrc_vec, method=interp)

#         # do interpolation
#         interp_result = data_w_coor(interp_points)
#         vec_result = vec_w_coor(interp_points)

#     else:
#         # no interpolation
#         # interp_result = orig_mrc_data[interp_points[:, 0].astype(np.int32),
#         #                               interp_points[:, 1].astype(np.int32),
#         #                               interp_points[:, 2].astype(np.int32)]
#         # vec_result = orig_mrc_vec[interp_points[:, 0].astype(np.int32),
#         #                           interp_points[:, 1].astype(np.int32),
#         #                           interp_points[:, 2].astype(np.int32)]

#         interp_result = orig_mrc_data[interp_points[:, 2].astype(np.int32),
#                                       interp_points[:, 1].astype(np.int32),
#                                       interp_points[:, 0].astype(np.int32)]
#         vec_result = orig_mrc_vec[interp_points[:, 2].astype(np.int32),
#                                   interp_points[:, 1].astype(np.int32),
#                                   interp_points[:, 0].astype(np.int32)]

#     # save interpolated data

#     new_data_array[z_valid_idx, y_valid_idx, x_valid_idx] = interp_result
#     new_vec_array[z_valid_idx, y_valid_idx, x_valid_idx] = np.einsum("ij, kj->ki", mtx, vec_result)

#     # new_data_array[x_valid_idx, y_valid_idx, z_valid_idx] = interp_result
#     # new_vec_array[x_valid_idx, y_valid_idx, z_valid_idx] = np.einsum("ij, kj->ki", mtx, vec_result)

#     return new_vec_array, new_data_array


# @numba.jit(nopython=True)
# def laplacian_filter(arr):
#     """A simple laplacian filter applied to an array with the kernel [[0, 1, 0], [1, -6, 1], [0, 1, 0]].
#
#     Args:
#         arr (numpy.array): the array to be filtered
#
#     Returns:
#         new_arr (numpy.array): the filtered array
#     """
#     xdim = arr.shape[0]
#     ydim = arr.shape[1]
#     zdim = arr.shape[2]
#     new_arr = np.zeros_like(arr)
#     for x in range(xdim):
#         for y in range(ydim):
#             for z in range(zdim):
#                 if arr[x][y][z] > 0:
#                     new_arr[x][y][z] = -6.0 * arr[x][y][z]
#                     if x + 1 < xdim:
#                         new_arr[x][y][z] += arr[x + 1][y][z]
#                     if x - 1 >= 0:
#                         new_arr[x][y][z] += arr[x - 1][y][z]
#                     if y + 1 < ydim:
#                         new_arr[x][y][z] += arr[x][y + 1][z]
#                     if y - 1 >= 0:
#                         new_arr[x][y][z] += arr[x][y - 1][z]
#                     if z + 1 < zdim:
#                         new_arr[x][y][z] += arr[x][y][z + 1]
#                     if z - 1 >= 0:
#                         new_arr[x][y][z] += arr[x][y][z - 1]
#     return new_arr

def next_fast_fft_len(size):
    a = 2
    while 1:
        if a > size:
            break
        a *= 2

    b = 3
    while 1:
        if b > size:
            break
        b *= 2

    b = 9
    while 1:
        if b > size:
            break
        b *= 2
    if a > b:
        return b
    else:
        return a

# @jit(nopython=True)
# def interpolate(data):
#     found = False
#     for i in range(data.shape[0]):
#         for j in range(data.shape[1]):
#             for k in range(data.shape[2]):
#                 if (data[i, j, k, 0] > 0):
#                     found = True
#                     if i % 2 == 0:
#                         ref_x = 0
#                     else:
#                         ref_x = 1
#
#                     if j % 2 == 0:
#                         ref_y = 0
#                     else:
#                         ref_y = 1
#
#                     if k % 2 == 0:
#                         ref_z = 0
#                     else:
#                         ref_z = 1
#
#                     #print(i, j, k)
#
#                 if found:
#                     break
#             if found:
#                 break
#         if found:
#             break
#
#     #print(ref_x, ref_y, ref_z)
#     #print((ref_x + 1) % 2, (ref_y + 1) % 2, (ref_z + 1) % 2)
#
#     #Change a single into the opposite of the reference at a time
#     for i in range((ref_x + 1) % 2, data.shape[0]-1, 2):
#         for j in range(ref_y, data.shape[1]-1, 2):
#             for k in range(ref_z, data.shape[2]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i+1, j, k, l]
#                         avg += data[i-1, j, k, l]
#                         avg = avg / float(2)
#                         data[i, j, k, l] = avg
#
#
#     for j in range((ref_y + 1) % 2, data.shape[1]-1, 2):
#         for i in range(ref_x, data.shape[0]-1, 2):
#             for k in range(ref_z, data.shape[2]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i, j+1, k, l]
#                         avg += data[i, j-1, k, l]
#                         avg = avg / float(2)
#                         data[i, j, k, l] = avg
#
#
#     for k in range((ref_z + 1) % 2, data.shape[2]-1, 2):
#         for i in range(ref_x, data.shape[0]-1, 2):
#             for j in range(ref_y, data.shape[1]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i, j, k+1, l]
#                         avg += data[i, j, k-1, l]
#                         avg = avg / float(2)
#                         data[i, j, k, l] = avg
#
#
#     #Change two indices into the opposite of the reference at a time
#     for i in range((ref_x + 1) % 2, data.shape[0]-1, 2):
#         for j in range((ref_y + 1) % 2, data.shape[1]-1, 2):
#             for k in range(ref_z, data.shape[2]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i+1, j, k, l]
#                         avg += data[i-1, j, k, l]
#                         avg += data[i, j+1, k, l]
#                         avg += data[i, j-1, k, l]
#                         avg = avg / float(4)
#                         data[i, j, k, l] = avg
#
#
#     for i in range((ref_x + 1) % 2, data.shape[0]-1, 2):
#         for k in range((ref_z + 1) % 2, data.shape[2]-1, 2):
#             for j in range(ref_y, data.shape[1]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i+1, j, k, l]
#                         avg += data[i-1, j, k, l]
#                         avg += data[i, j, k+1, l]
#                         avg += data[i, j, k-1, l]
#                         avg = avg / float(4)
#                         data[i, j, k, l] = avg
#
#
#     for j in range((ref_y + 1) % 2, data.shape[1]-1, 2):
#         for k in range((ref_z + 1) % 2, data.shape[2]-1, 2):
#             for i in range(ref_x, data.shape[0]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i, j+1, k, l]
#                         avg += data[i, j-1, k, l]
#                         avg += data[i, j, k+1, l]
#                         avg += data[i, j, k-1, l]
#                         avg = avg / float(4)
#                         data[i, j, k, l] = avg
#
#
#     #Change all indices into the opposite of the reference
#     for i in range((ref_x + 1) % 2, data.shape[0]-1, 2):
#         for j in range((ref_y + 1) % 2, data.shape[1]-1, 2):
#             for k in range((ref_z + 1) % 2, data.shape[2]-1, 2):
#                 for l in range(data.shape[3]):
#                     if (data[i, j, k, l] == 0):
#                         avg = 0
#                         avg += data[i+1, j, k, l]
#                         avg += data[i-1, j, k, l]
#                         avg += data[i, j+1, k, l]
#                         avg += data[i, j-1, k, l]
#                         avg += data[i, j, k+1, l]
#                         avg += data[i, j, k-1, l]
#                         avg = avg / float(6)
#                         data[i, j, k, l] = avg

# def conv_prediction_to_mrc(arr,reffile,keyidx,outmapfile):
#     with mrcfile.open(reffile) as mrc:
#         nx,ny,nz = arr.shape[0], arr.shape[1], arr.shape[2]
#         mx,my,mz,cella = mrc.header.mx,mrc.header.my,mrc.header.mz, mrc.header.cella
#
#
#         mrc_new = mrcfile.new(outmapfile,overwrite=True)
#
#         mrc_new.set_data(np.zeros((nz, ny, nx), dtype=np.float32))
#         mrc_new.header.nxstart=arr.shape[0]
#         mrc_new.header.nystart=arr.shape[1]
#         mrc_new.header.nzstart=arr.shape[2]
#         mrc_new.header.origin.x = mrc.header.origin.x
#         mrc_new.header.origin.y = mrc.header.origin.y
#         mrc_new.header.origin.z = mrc.header.origin.z
#         mrc_new.header.cella['x'] = cella['x']
#         mrc_new.header.cella['y'] = cella['y']
#         mrc_new.header.cella['z'] = cella['z']
#
#         print(nx, ny, nz)
#         for i in range(nx - 2):
#             for j in range(ny - 2):
#                 for k in range(nz - 2):
#                     mrc_new.data[k, j, i] = arr[i, j, k, keyidx]
#
#         vsize=mrc_new.voxel_size
#         vsize.flags.writeable = True
#         mrc_new.voxel_size=vsize
#
#         mrc_new.update_header_stats()
#
#         print("original", mrc.voxel_size)
#         mrc.print_header()
#
#         print()
#
#         print("new", mrc_new.voxel_size)
#         mrc_new.print_header()
#
#         print()
#
#         mrc_new.close()

# def find_best_trans_mixed(vec_fft_results, prob_fft_results, alpha, vstd, vave, pstd, pave):
#     """
#     It takes the sum of the two arrays, normalizes them, mixes them, and then finds the best translation
#     :param vec_fft_results: the results of the FFT on the vectorized image
#     :param prob_fft_results: the FFT of the probability map
#     :param alpha: the weight of the probability map
#     :param vstd: standard deviation of the vector fft results
#     :param vave: the mean of the vector fft results
#     :param pstd: standard deviation of the secondary structure matching score fft results
#     :param pave: the mean of the secondary structure matching score fft results
#     :return: The best score and the translation that produced it.
#     """
#     sum_arr_v = vec_fft_results[0] + vec_fft_results[1] + vec_fft_results[2]
#
#     # sum_arr_v = sum(vec_fft_results)
#     # sum_arr_p = sum(prob_fft_results)
#     sum_arr_p = prob_fft_results[0] + prob_fft_results[1] + prob_fft_results[2]
#
#     # z-score normalization
#     sum_arr_v = (sum_arr_v - vave) / vstd
#     sum_arr_p = (sum_arr_p - pave) / pstd
#
#     # mix the two arrays
#     sum_arr_mixed = (1 - alpha) * sum_arr_v + alpha * sum_arr_p
#
#     # find the best translation
#     best_score = sum_arr_mixed.max()
#     best_trans = np.unravel_index(sum_arr_mixed.argmax(), sum_arr_mixed.shape)
#
#     return best_score, best_trans

# import copy
#
# import numpy as np
# from pyfftw import pyfftw
#
# from search import fft_search_best_dot
# from utils import new_rot_mrc_prob, new_rot_mrc
#
# from scipy.spatial.transform import Rotation as R
#
#
# def eval_score_orig(mrc_target, mrc_search, angle, trans, dot_score_ave, dot_score_std):
#     trans = [int(i) for i in trans]
#
#     # Function to evaluate the DOT score for input rotation angle and translation
#
#     # init rotation grid
#     search_pos_grid = (
#         np.mgrid[
#             0 : mrc_search.data.shape[0],
#             0 : mrc_search.data.shape[0],
#             0 : mrc_search.data.shape[0],
#         ]
#         .reshape(3, -1)
#         .T
#     )
#
#     # init the target map vectors
#     x1 = copy.deepcopy(mrc_target.vec[:, :, :, 0])
#     y1 = copy.deepcopy(mrc_target.vec[:, :, :, 1])
#     z1 = copy.deepcopy(mrc_target.vec[:, :, :, 2])
#
#     rd3 = 1.0 / mrc_target.data.size
#
#     # init fft transformation for the target map
#
#     X1 = np.fft.rfftn(x1)
#     X1 = np.conj(X1)
#     Y1 = np.fft.rfftn(y1)
#     Y1 = np.conj(Y1)
#     Z1 = np.fft.rfftn(z1)
#     Z1 = np.conj(Z1)
#     target_list = [X1, Y1, Z1]
#
#     # fftw plans initialization
#     a = pyfftw.empty_aligned(mrc_search.vec[..., 0].shape, dtype="float32")
#     b = pyfftw.empty_aligned((a.shape[0], a.shape[1], a.shape[2] // 2 + 1), dtype="complex64")
#     c = pyfftw.empty_aligned(mrc_search.vec[..., 0].shape, dtype="float32")
#
#     fft_object = pyfftw.FFTW(a, b, axes=(0, 1, 2))
#     ifft_object = pyfftw.FFTW(b, c, direction="FFTW_BACKWARD", axes=(0, 1, 2), normalise_idft=False)
#
#     # init the rotation matrix by euler angle
#     rot_mtx = R.from_euler("xyz", angle, degrees=True).as_matrix()
#
#     new_vec, new_data = new_rot_mrc(mrc_search.data, mrc_search.vec, rot_mtx, search_pos_grid)
#
#     x2 = new_vec[..., 0]
#     y2 = new_vec[..., 1]
#     z2 = new_vec[..., 2]
#
#     query_list_vec = [x2, y2, z2]
#
#     fft_result_list_vec = fft_search_best_dot(target_list[:3], query_list_vec, a, b, c, fft_object, ifft_object)
#
#     sum_arr = np.zeros_like(fft_result_list_vec[0])
#     for arr in fft_result_list_vec:
#         sum_arr = sum_arr + arr
#
#     dot_score = sum_arr[trans[0]][trans[1]][trans[2]] * rd3
#     dot_score = (dot_score - dot_score_ave) / dot_score_std
#
#     return dot_score
#
#
# def eval_score_mix(
#     mrc_target,
#     mrc_input,
#     mrc_P1,
#     mrc_P2,
#     mrc_P3,
#     mrc_P4,
#     mrc_search_p1,
#     mrc_search_p2,
#     mrc_search_p3,
#     mrc_search_p4,
#     angle_list,
#     trans_list,
#     vstd,
#     vave,
#     pstd,
#     pave,
#     mix_score_ave,
#     mix_score_std,
#     alpha,
# ):
#     # convert the translation list to integer
#     for trans in trans_list:
#         trans = [int(i) for i in trans]
#
#     # init the target map vectors
#     x1 = copy.deepcopy(mrc_target.vec[:, :, :, 0])
#     y1 = copy.deepcopy(mrc_target.vec[:, :, :, 1])
#     z1 = copy.deepcopy(mrc_target.vec[:, :, :, 2])
#
#     p1 = copy.deepcopy(mrc_P1.data)
#     p2 = copy.deepcopy(mrc_P2.data)
#     p3 = copy.deepcopy(mrc_P3.data)
#     p4 = copy.deepcopy(mrc_P4.data)
#
#     # Score normalization constant
#
#     rd3 = 1.0 / (mrc_target.xdim**3)
#
#     # Calculate the FFT results for target map
#
#     X1 = np.fft.rfftn(x1)
#     X1 = np.conj(X1)
#     P1 = np.fft.rfftn(p1)
#     P1 = np.conj(P1)
#     P2 = np.fft.rfftn(p2)
#     P2 = np.conj(P2)
#     P3 = np.fft.rfftn(p3)
#     P3 = np.conj(P3)
#     P4 = np.fft.rfftn(p4)
#     P4 = np.conj(P4)
#
#     Y1 = np.fft.rfftn(y1)
#     Y1 = np.conj(Y1)
#     Z1 = np.fft.rfftn(z1)
#     Z1 = np.conj(Z1)
#
#     # Compose target result list
#
#     target_list = [X1, Y1, Z1, P1, P2, P3, P4]
#
#     # fftw plans initialization
#     a = pyfftw.empty_aligned(mrc_search_p1.data.shape, dtype="float32")
#     b = pyfftw.empty_aligned((a.shape[0], a.shape[1], a.shape[2] // 2 + 1), dtype="complex64")
#     c = pyfftw.empty_aligned(mrc_search_p1.data.shape, dtype="float32")
#
#     fft_object = pyfftw.FFTW(a, b, axes=(0, 1, 2))
#     ifft_object = pyfftw.FFTW(b, c, direction="FFTW_BACKWARD", axes=(0, 1, 2), normalise_idft=False)
#
#     # init rotation grid
#     search_pos_grid = (
#         np.mgrid[
#             0 : mrc_input.data.shape[0],
#             0 : mrc_input.data.shape[0],
#             0 : mrc_input.data.shape[0],
#         ]
#         .reshape(3, -1)
#         .T
#     )
#
#     mix_score_list = []
#
#     for angle, trans in zip(angle_list, trans_list):
#         trans = [int(i) for i in trans]
#
#         rot_mtx = R.from_euler("xyz", angle, degrees=True).as_matrix()
#
#         r_vec, r_data, rp1, rp2, rp3, rp4 = new_rot_mrc_prob(
#             mrc_input.data,
#             mrc_input.vec,
#             mrc_search_p1.data,
#             mrc_search_p2.data,
#             mrc_search_p3.data,
#             mrc_search_p4.data,
#             rot_mtx,
#             search_pos_grid,
#         )
#
#         # extract XYZ components from vector representation
#         x2 = r_vec[..., 0]
#         y2 = r_vec[..., 1]
#         z2 = r_vec[..., 2]
#
#         p21 = rp1
#         p22 = rp2
#         p23 = rp3
#         p24 = rp4
#
#         query_list_vec = [x2, y2, z2]
#         query_list_prob = [p21, p22, p23, p24]
#
#         fft_result_list_vec = fft_search_best_dot(target_list[:3], query_list_vec, a, b, c, fft_object, ifft_object)
#         fft_result_list_prob = fft_search_best_dot(target_list[3:], query_list_prob, a, b, c, fft_object, ifft_object)
#
#         sum_arr_v = fft_result_list_vec[0] + fft_result_list_vec[1] + fft_result_list_vec[2]
#
#         sum_arr_p = fft_result_list_prob[0] + fft_result_list_prob[1] + fft_result_list_prob[2]
#
#         sum_arr_v = (sum_arr_v - vave) / vstd
#         sum_arr_p = (sum_arr_p - pave) / pstd
#
#         sum_arr_mixed = (1 - alpha) * sum_arr_v + alpha * sum_arr_p
#
#         best_score = sum_arr_mixed.max()
#         best_trans = np.unravel_index(sum_arr_mixed.argmax(), sum_arr_mixed.shape)
#
#         # print(best_score, best_trans)
#
#         mix_score = sum_arr_mixed[trans[0]][trans[1]][trans[2]]
#
#         # print(mix_score)
#
#         mix_score = (mix_score - mix_score_ave) / mix_score_std
#
#         mix_score_list.append(mix_score)
#
#         # print(f"vave={vave}, vstd={vstd}, pave={pave}, pstd={pstd}, mix_score_ave={mix_score_ave}, mix_score_std={mix_score_std}")
#
#     return mix_score_list


# @jit(nopython=True, nogil=True)
# def do_resample_and_vec(
#         src_xwidth,
#         src_orig,
#         src_dims,
#         src_data,
#         dest_xwidth,
#         dest_orig,
#         new_dim,
#         dreso,
#         ss_data,
#         old_pos_list,
#         new_pos_list,
# ):
#     """
#     The do_resample_and_vec function takes in the following parameters:
#         src_xwidth - The width of the source map.
#         src_orig - The origin of the source map.
#         src_dims - The dimensions of the source map.
#         src_data - A array contains the density data. The shape of the array is (xdim, ydim, zdim).
#
#     :param src_xwidth: Calculate the grid step
#     :param src_orig: Define the origin of the source volume
#     :param src_dims: Determine the size of the source volume
#     :param src_data: Data of the source volume
#     :param dest_xwidth: Determine the size of the new grid
#     :param dest_orig: Define the origin of the destination volume
#     :param new_dim: Determine the size of the new volume
#     :param dreso: Determine the Guassian filter size
#     :param ss_data: data of the secondary structure
#     :param : Determine the resolution of the new map
#     :return: The density map, the unit vectors and the secondary structure data
#     """
#
#     gstep = src_xwidth
#     fs = (dreso / gstep) * 0.5
#     fs = fs ** 2
#     fsiv = 1.0 / fs
#     fmaxd = (dreso / gstep) * 2.0
#     print("#maxd=", fmaxd)
#     print("#fsiv=", fsiv)
#
#     dest_vec = np.zeros((new_dim, new_dim, new_dim, 3), dtype="float32")
#     dest_data = np.zeros((new_dim, new_dim, new_dim), dtype="float32")
#     dest_ss_data = np.zeros((new_dim, new_dim, new_dim, 4), dtype="float32")
#
#     for new_pos, old_pos in zip(new_pos_list, old_pos_list):
#         stp = np.maximum(old_pos - fmaxd, 0).astype(np.int32)
#         endp = np.minimum(old_pos + fmaxd + 1, src_dims).astype(np.int32)
#
#         # compute the total density
#         dtotal, v = calc(stp, endp, old_pos, src_data, fsiv)
#
#         if ss_data is not None:
#             dest_ss_data[new_pos[0], new_pos[1], new_pos[2], 0], _ = calc(stp, endp, old_pos, ss_data[..., 0], fsiv, False)
#             dest_ss_data[new_pos[0], new_pos[1], new_pos[2], 1], _ = calc(stp, endp, old_pos, ss_data[..., 1], fsiv, False)
#             dest_ss_data[new_pos[0], new_pos[1], new_pos[2], 2], _ = calc(stp, endp, old_pos, ss_data[..., 2], fsiv, False)
#             dest_ss_data[new_pos[0], new_pos[1], new_pos[2], 3], _ = calc(stp, endp, old_pos, ss_data[..., 3], fsiv, False)
#
#         if np.isclose(dtotal, 0.0):
#             continue
#
#         dest_data[new_pos[0], new_pos[1], new_pos[2]] = dtotal
#         # tmpcd = pos2 / dtotal - old_pos
#         # dvec = np.sqrt(np.sum(tmpcd ** 2))
#         # if dvec == 0:
#         #     dvec = 1.0
#         # dest_vec[new_pos[0], new_pos[1], new_pos[2]] = tmpcd / dvec
#         dest_vec[new_pos[0], new_pos[1], new_pos[2]] = v
#
#     return dest_data, dest_vec, dest_ss_data

# @jit(nopython=True, nogil=True)
# def calc(stp, endp, pos, mrc1_data, fsiv, calc_vector=True):
#     """Vectorized version of calc"""
#
#     xx = np.arange(stp[0], endp[0], 1)
#     yy = np.arange(stp[1], endp[1], 1)
#     zz = np.arange(stp[2], endp[2], 1)
#
#     xx = np.expand_dims(xx, axis=1)
#     xx = np.expand_dims(xx, axis=1)
#     yy = np.expand_dims(yy, axis=1)
#     yy = np.expand_dims(yy, axis=0)
#     zz = np.expand_dims(zz, axis=0)
#     zz = np.expand_dims(zz, axis=0)
#
#     # calculate the distance between the center of the voxel and the center of the particle
#     d2 = (xx - pos[0])**2 + (yy - pos[1])**2 + (zz - pos[2])**2
#
#     # calculate the density and vector in resized map using Gaussian interpolation in original MRC density map
#     d = np.exp(-1.5 * d2 * fsiv) * mrc1_data[stp[0]:endp[0], stp[1]:endp[1], stp[2]:endp[2]]
#     dtotal = np.sum(d)
#
#     if calc_vector:
#         # calculate the vector
#         pos2 = np.array([np.sum(d * xx), np.sum(d * yy), np.sum(d * zz)])
#
#         pos2 = pos2 / dtotal
#
#         v = (pos2 - pos) / np.sqrt(np.sum((pos2 - pos)**2))
#     else:
#         v = None
#
#     return dtotal, v