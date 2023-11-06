import numpy as np


def format_score_result(result, ave, std):
    return (
        f"Rotation {result['angle']}, DOT Score: {result['vec_score']}, DOT Trans: {result['vec_trans']}, "
        + f"Prob Score: {result['prob_score']}, Prob Trans: {result['prob_trans']}, "
        + f"MIX Score: {result['mixed_score']}, MIX Trans: {result['mixed_trans']}, "
        + f"Normalized Mix Score: {(result['mixed_score'] - ave) / std}"
    )


def euler_to_mtx(ang):
    # extrinsic rotation in xyz order or intrinsic rotation in zyx order
    # ang: [alpha, beta, gamma]

    import torch

    cos_a = torch.cos(ang[0])
    sin_a = torch.sin(ang[0])
    cos_b = torch.cos(ang[1])
    sin_b = torch.sin(ang[1])
    cos_c = torch.cos(ang[2])
    sin_c = torch.sin(ang[2])
    one = torch.ones_like(ang[0])
    zero = torch.zeros_like(ang[0])

    mtx_z = torch.stack((cos_c, -sin_c, zero, sin_c, cos_c, zero, zero, zero, one)).reshape(3, 3)
    mtx_y = torch.stack((cos_b, zero, sin_b, zero, one, zero, -sin_b, zero, cos_b)).reshape(3, 3)
    mtx_x = torch.stack((one, zero, zero, zero, cos_a, -sin_a, zero, sin_a, cos_a)).reshape(3, 3)

    return mtx_z @ mtx_y @ mtx_x


def get_score(ref_map, tgt_map_data, tgt_map_vec, trans):

    ave1 = ref_map.ave
    ave2 = np.mean(tgt_map_data[tgt_map_data > 0])

    std1 = ref_map.std
    std2 = np.linalg.norm(tgt_map_data[tgt_map_data > 0])

    pstd1 = ref_map.std_norm_ave
    pstd2 = np.linalg.norm(tgt_map_data[tgt_map_data > 0] - ave2)

    ref_map_data = ref_map.new_data
    ref_map_vec = ref_map.vec

    dim = ref_map_data.shape[0]
    total = 0

    t = np.array(trans)
    if trans[0] > 0.5 * dim:
        t[0] -= dim
    if trans[1] > 0.5 * dim:
        t[1] -= dim
    if trans[2] > 0.5 * dim:
        t[2] -= dim

    target_pos = np.array(
        np.meshgrid(
            np.arange(dim),
            np.arange(dim),
            np.arange(dim),
        )
    ).T.reshape(-1, 3)

    search_pos = target_pos + t

    total += np.count_nonzero(ref_map_data[target_pos[:, 0], target_pos[:, 1], target_pos[:, 2]])

    combined_arr = np.hstack((target_pos, search_pos))

    combined_arr = combined_arr[
        (combined_arr[:, 3] >= 0)
        & (combined_arr[:, 4] >= 0)
        & (combined_arr[:, 5] >= 0)
        & (combined_arr[:, 3] < dim)
        & (combined_arr[:, 4] < dim)
        & (combined_arr[:, 5] < dim)
    ]

    target_pos = combined_arr[:, 0:3]
    search_pos = combined_arr[:, 3:6]

    d1 = ref_map_data[target_pos[:, 0], target_pos[:, 1], target_pos[:, 2]]
    d2 = tgt_map_data[search_pos[:, 0], search_pos[:, 1], search_pos[:, 2]]

    d1 = np.where(d1 <= 0, 0.0, d1)  # trim negative values
    d2 = np.where(d2 <= 0, 0.0, d1)  # trim negative values

    pd1 = np.where(d1 <= 0, 0.0, d1 - ave1)  # trim negative values
    pd2 = np.where(d2 <= 0, 0.0, d2 - ave2)  # trim negative values

    cc = np.sum(np.multiply(d1, d2))  # cross correlation
    pcc = np.sum(np.multiply(pd1, pd2))  # Pearson cross correlation

    target_zero_mask = ref_map_data[target_pos[:, 0], target_pos[:, 1], target_pos[:, 2]] == 0
    target_non_zero_mask = ref_map_data[target_pos[:, 0], target_pos[:, 1], target_pos[:, 2]] > 0
    search_non_zero_mask = tgt_map_data[search_pos[:, 0], search_pos[:, 1], search_pos[:, 2]] > 0
    search_non_zero_count = np.count_nonzero(np.multiply(target_zero_mask, search_non_zero_mask))

    trimmed_target_vec = ref_map_vec[target_pos[:, 0], target_pos[:, 1], target_pos[:, 2]]
    trimmed_search_vec = tgt_map_vec[search_pos[:, 0], search_pos[:, 1], search_pos[:, 2]]

    total += search_non_zero_count

    sco_arr = np.zeros_like(tgt_map_data)
    sco = np.einsum("ij,ij->i", trimmed_target_vec, trimmed_search_vec)
    sco_arr[search_pos[:, 0], search_pos[:, 1], search_pos[:, 2]] = sco
    sco_sum = np.sum(sco_arr)
    Nm = np.count_nonzero(np.multiply(target_non_zero_mask, search_non_zero_mask))

    overlap = float(Nm) / float(total)
    cc = cc / (std1 * std2)  # cross correlation
    pcc = pcc / (pstd1 * pstd2)  # Pearson cross correlation
    dot = sco_sum

    return sco_arr, overlap, cc, pcc, Nm, total, dot
