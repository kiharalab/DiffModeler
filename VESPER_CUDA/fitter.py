import concurrent.futures
import os
from datetime import datetime
from itertools import product

import numpy as np
import torch
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import laplace
from tqdm import tqdm

from utils.fileio import save_rotated_pdb, save_rotated_mrc, save_vec_as_pdb
from utils.utils import euler_to_mtx, get_score


class MapFitter:
    def __init__(
        self,
        ref,
        tgt,
        ang_interval,
        mode,
        remove_dup,
        ldp_path,
        backbone_path,
        input_pdb,
        threads,
        gpu,
        device,
        topn=10,
        outdir=None,
        save_mrc=False,
        alpha=None,
        confine_angles=None,
        save_vec=False,
    ):

        print("###Initializing fitter###")

        self.ref_map = ref
        self.tgt_map = tgt
        self.ang_interval = ang_interval
        self.mode = mode
        self.remove_dup = remove_dup
        self.ldp_path = ldp_path
        self.backbone_path = backbone_path
        self.input_pdb = input_pdb
        self.threads = threads
        self.gpu = gpu
        self.device = device
        self.topn = topn
        self.outdir = outdir
        self.save_mrc = save_mrc
        self.save_vec = save_vec
        self.angle_comb = []

        self.result_list = None
        self.refined_list = None
        self.final_list = None
        self.score_ave = None
        self.score_std = None
        self.alpha = alpha
        self.confine_angles = confine_angles

        if self.outdir is None and self.input_pdb is not None:
            self.outdir = os.path.join("./outputs", "VESPER_RUN_" + datetime.now().strftime("%m%d_%H%M%S"))
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir, exist_ok=True)

        self.ldp_recall_mode = (ldp_path is not None) and (backbone_path is not None)
        self.ss_mix_score_mode = (
            (alpha is not None) and (self.ref_map.new_ss_data is not None) and (self.tgt_map.new_ss_data is not None)
        )

        # init rotation grid
        self.search_pos_grid = (
            np.mgrid[
                0 : self.tgt_map.new_dim,
                0 : self.tgt_map.new_dim,
                0 : self.tgt_map.new_dim,
            ]
            .reshape(3, -1)
            .T
        )

        # calculate combination of rotation angles
        self._calc_angle_comb()
        # self._calc_angle_comb_quat()
        self.total_rotations = len(self.angle_comb)

        # init the target map vectors
        ref_x_real = self.ref_map.vec[:, :, :, 0]

        # Postprocessing for other modes
        if mode == "Overlap":
            ref_x_real = np.where(self.ref_map.new_data > 0, 1.0, 0.0)
        elif mode == "CC":
            ref_x_real = np.where(self.ref_map.new_data > 0, self.ref_map.new_data, 0.0)
        elif mode == "PCC":
            ref_x_real = np.where(self.ref_map.new_data > 0, self.ref_map.new_data - self.ref_map.ave, 0.0)
        elif mode == "Laplacian":
            ref_x_real = laplace(self.ref_map.new_data, mode="constant", cval=0.0)

        # init fft transformation for the target map
        ref_x_fourier = np.fft.rfftn(ref_x_real)
        ref_x_fourier = np.conj(ref_x_fourier)

        self.ref_map_fft_list = [ref_x_fourier]

        # init fft transformation for the target map
        if mode == "VecProduct":
            ref_y_real = self.ref_map.vec[:, :, :, 1]
            ref_z_real = self.ref_map.vec[:, :, :, 2]
            ref_y_fourier = np.fft.rfftn(ref_y_real)
            ref_y_fourier = np.conj(ref_y_fourier)
            ref_z_fourier = np.fft.rfftn(ref_z_real)
            ref_z_fourier = np.conj(ref_z_fourier)
            self.ref_map_fft_list = [ref_x_fourier, ref_y_fourier, ref_z_fourier]
            if self.ss_mix_score_mode:
                ref_ss_c_real = self.ref_map.new_ss_data[..., 0]  # coil
                ref_ss_b_real = self.ref_map.new_ss_data[..., 1]  # beta
                ref_ss_a_real = self.ref_map.new_ss_data[..., 2]  # alpha
                res_ss_n_real = self.ref_map.new_ss_data[..., 3]  # nucleotide
                ref_ss_c_fourier = np.fft.rfftn(ref_ss_c_real)
                ref_ss_c_fourier = np.conj(ref_ss_c_fourier)
                ref_ss_b_fourier = np.fft.rfftn(ref_ss_b_real)
                ref_ss_b_fourier = np.conj(ref_ss_b_fourier)
                ref_ss_a_fourier = np.fft.rfftn(ref_ss_a_real)
                ref_ss_a_fourier = np.conj(ref_ss_a_fourier)
                res_ss_n_fourier = np.fft.rfftn(res_ss_n_real)
                res_ss_n_fourier = np.conj(res_ss_n_fourier)
                self.ref_map_fft_list.extend([ref_ss_c_fourier, ref_ss_b_fourier, ref_ss_a_fourier, res_ss_n_fourier])

        # ldp recall mode init
        if self.ldp_recall_mode:
            assert self.gpu, "LDP recall mode only works with GPU"
            import torch

            # get atom coords from ldp
            ldp_atoms = []
            with open(ldp_path) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        ldp_atoms.append(np.array((float(line[30:38]), float(line[38:46]), float(line[46:54]))))

            assert len(ldp_atoms) > 0, "No points found in LDP file."
            self.ldp_atoms = torch.from_numpy(np.array(ldp_atoms)).to(self.device)

            # get ca atoms from backbone
            backbone_ca = []
            with open(backbone_path) as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":  # only CA atoms
                        # if tokens[0] == "ATOM": # all atoms
                        backbone_ca.append(np.array((float(line[30:38]), float(line[38:46]), float(line[46:54]))))

            assert len(backbone_ca) > 0, "No CA atoms found in backbone file."
            self.backbone_ca = torch.from_numpy(np.array(backbone_ca)).to(self.device)

        if self.gpu:
            # move everything to GPU
            import torch

            self.tgt_map.new_data_gpu = torch.from_numpy(self.tgt_map.new_data).to(self.device).share_memory_()
            self.tgt_map.vec_gpu = torch.from_numpy(self.tgt_map.vec).to(self.device).share_memory_()
            self.tgt_map.search_pos_grid_gpu = torch.from_numpy(self.search_pos_grid).to(self.device).share_memory_()

            if self.ss_mix_score_mode:
                self.tgt_map.new_ss_data_gpu = torch.from_numpy(self.tgt_map.new_ss_data).to(self.device).share_memory_()

            self.ref_map_fft_list_gpu = [
                torch.from_numpy(fft_arr).to(self.device).share_memory_() for fft_arr in self.ref_map_fft_list
            ]
        else:
            import pyfftw.config

            pyfftw.config.PLANNER_EFFORT = "FFTW_MEASURE"
            pyfftw.config.NUM_THREADS = max(os.cpu_count() - 2, 2)  # Maybe the CPU is sweating too much?

    def fit_ss(self):
        print("###Start Searching###")
        self.result_list = []

        with tqdm(total=len(self.angle_comb)) as pbar:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
                futures = {
                    executor.submit(
                        self._rot_and_search_fft_ss,
                        rot_ang,
                        False,
                        None,
                        None,
                        None,
                        None,
                    ): rot_ang
                    for rot_ang in self.angle_comb
                }
                for future in concurrent.futures.as_completed(futures):
                    rot_ang = futures[future]
                    result = future.result()
                    pbar.update(1)
                    self.result_list.append(
                        {
                            "angle": rot_ang,
                            "vec_score": result[0] / self.ref_map.new_data.size,
                            "vec_vox_trans": result[1],
                            "ss_score": result[2] / self.ref_map.new_data.size,
                            "ss_vox_trans": result[3],
                        }
                    )

        # gather stats
        vec_score_arr = np.array([item["vec_score"] for item in self.result_list])
        self.vec_ave = np.mean(vec_score_arr)
        self.vec_std = np.std(vec_score_arr)
        ss_score_arr = np.array([item["ss_score"] for item in self.result_list])
        self.ss_ave = np.mean(ss_score_arr)
        self.ss_std = np.std(ss_score_arr)
        print(f"Vec Score Mean: {self.vec_ave}, StDev: {self.vec_std}")
        print(f"SS Score Mean: {self.ss_ave}, StDev: {self.ss_std}")
        print()

        self.result_list = []

        with tqdm(total=len(self.angle_comb)) as pbar:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
                futures = {
                    executor.submit(
                        self._rot_and_search_fft_ss,
                        rot_ang,
                        False,
                        self.vec_ave,
                        self.vec_std,
                        self.ss_ave,
                        self.ss_std,
                    ): rot_ang
                    for rot_ang in self.angle_comb
                }
                for future in concurrent.futures.as_completed(futures):
                    rot_ang = futures[future]
                    result = future.result()
                    pbar.update(1)
                    self.result_list.append(
                        {
                            "angle": rot_ang,
                            "vec_score": result[0] / self.ref_map.new_data.size,
                            "vec_vox_trans": result[1],
                            "ss_score": result[2] / self.ref_map.new_data.size,
                            "ss_vox_trans": result[3],
                            "mix_score": result[4] / self.ref_map.new_data.size,
                            "vox_trans": result[5],
                            "real_trans": result[6],
                        }
                    )

        # rank by mix score

        self.result_list = sorted(self.result_list, key=lambda x: x["mix_score"], reverse=True)

        # remove duplicates

        if self.remove_dup:
            self._remove_dup_results()
            print("#Non-duplicate count: " + str(len(self.result_list)))
            print()

        # refine

        if self.ang_interval >= 5:
            self.refine_ss(2)

        if self.refined_list:
            self.final_list = self.refined_list[: self.topn]
        else:
            self.final_list = self.result_list[: self.topn]

        for i, item in enumerate(self.final_list):
            print(
                "\n#" + str(i),
                "Rotation=",
                "(" + str(item["angle"][0]),
                str(item["angle"][1]),
                str(item["angle"][2]) + ")",
                "Translation=",
                "(" + "{:.3f}".format(item["real_trans"][0]),
                "{:.3f}".format(item["real_trans"][1]),
                "{:.3f}".format(item["real_trans"][2]) + ")",
            )
            print("Mix Score=", "{:.6f}".format(item["mix_score"]))
            print("Vec Score=", "{:.6f}".format(item["vec_score"]))
            print("SS Score=", "{:.6f}".format(item["ss_score"]))

        if self.input_pdb is not None:
            self._save_topn_pdb()

    def fit(self):
        print("###Start Searching###")
        self.result_list = []
        with tqdm(total=len(self.angle_comb)) as pbar:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
                futures = {
                    executor.submit(
                        self._rot_and_search_fft,
                        rot_ang,
                        False,
                    ): rot_ang
                    for rot_ang in self.angle_comb
                }
                for future in concurrent.futures.as_completed(futures):
                    rot_ang = futures[future]
                    result = future.result()
                    pbar.update(1)
                    self.result_list.append(
                        {
                            "angle": rot_ang,
                            "score": result[0] / self.ref_map.new_data.size,
                            "vox_trans": result[1],
                        }
                    )
                for result in self.result_list:
                    result["real_trans"] = self._convert_trans(result["angle"], result["vox_trans"])

        # sort the result list
        self.result_list.sort(key=lambda x: x["score"], reverse=True)

        print("###Preliminary Results Summary###")
        self.score_ave, self.score_std = self._print_result_stats(self.result_list, return_stats=True)
        print()
        # for i, item in enumerate(self.result_list[:10]):
        #     self._print_result_item(item, i)

        # remove duplicates
        if self.remove_dup:
            self._remove_dup_results()
            print("#Non-duplicate count: " + str(len(self.result_list)))
            print()

        # calculate ldp recall if specified
        if self.ldp_recall_mode:
            self._calc_ldp_recall(self.result_list, sort=True)
            print()

        if self.ang_interval >= 5:
            self.refine(2, self.topn, sort_by_ldp_recall=self.ldp_recall_mode)

        if self.refined_list:
            self.final_list = self.refined_list[: self.topn]
        else:
            self.final_list = self.result_list[: self.topn]

        print("###Final Results###")
        self._print_result_stats(self.final_list)

        # print individual stats
        for i, item in enumerate(self.final_list):
            self._print_result_item(item, i)

        # save rotated pdb structure for visualization
        if self.input_pdb is not None:
            self._save_topn_pdb()

        # save vectors as pdb for visualization
        if self.save_vec:
            self._retrieve_data(self.final_list)
            self._save_topn_vec_as_pdb()

        if self.save_mrc:
            self._save_topn_mrc()

    def refine(self, ang_interval, top_n, sort_by_ldp_recall=False):
        print("###Start Refining###")
        self.refined_list = []
        top_n_list = self.result_list[:top_n]
        for result in tqdm(top_n_list, desc="Refining Top N", position=0):
            curr_result_list = []

            # compose angle list using the interval
            x_list = range(int(result["angle"][0]) - 5, int(result["angle"][0]) + 6, ang_interval)
            y_list = range(int(result["angle"][1]) - 5, int(result["angle"][1]) + 6, ang_interval)
            z_list = range(int(result["angle"][2]) - 5, int(result["angle"][2]) + 6, ang_interval)
            curr_refine_ang_list = np.array(list(product(x_list, y_list, z_list))).astype(np.float32)

            # make sure the angles are in the range of 0-360
            curr_refine_ang_list[curr_refine_ang_list < 0] += 360
            curr_refine_ang_list[curr_refine_ang_list > 360] -= 360

            executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.threads)
            with tqdm(total=len(curr_refine_ang_list), position=1, leave=False) as pbar:
                futures = {
                    executor.submit(
                        self._rot_and_search_fft,
                        rot_ang,
                        False,
                    ): rot_ang
                    for rot_ang in curr_refine_ang_list
                }
                for future in concurrent.futures.as_completed(futures):
                    rot_ang = futures[future]
                    result = future.result()
                    pbar.update(1)
                    curr_result_list.append(
                        {
                            "angle": rot_ang,
                            "score": result[0] / self.ref_map.new_data.size,
                            "vox_trans": result[1],
                        }
                    )
                for result in curr_result_list:
                    result["real_trans"] = self._convert_trans(result["angle"], result["vox_trans"])
                if self.ldp_recall_mode:
                    self._calc_ldp_recall(curr_result_list, progress_bar=False)
            if sort_by_ldp_recall and self.ldp_recall_mode:
                self.refined_list.append(max(curr_result_list, key=lambda x: x["ldp_recall"]))
            else:
                self.refined_list.append(max(curr_result_list, key=lambda x: x["score"]))
            # close the executor
            executor.shutdown(wait=True)
        # sort the refined list
        if sort_by_ldp_recall and self.ldp_recall_mode:
            self.refined_list.sort(key=lambda x: x["ldp_recall"], reverse=True)
        else:
            self.refined_list.sort(key=lambda x: x["score"], reverse=True)

    def refine_ss(self, ang_interval):
        print("###Start Refining###")
        self.refined_list = []
        top_n_list = self.result_list[: self.topn]
        for result in tqdm(top_n_list, desc="Refining Top N", position=0):
            curr_result_list = []

            # compose angle list using the interval
            x_list = range(int(result["angle"][0]) - 5, int(result["angle"][0]) + 6, ang_interval)
            y_list = range(int(result["angle"][1]) - 5, int(result["angle"][1]) + 6, ang_interval)
            z_list = range(int(result["angle"][2]) - 5, int(result["angle"][2]) + 6, ang_interval)
            curr_refine_ang_list = np.array(list(product(x_list, y_list, z_list))).astype(np.float32)

            # make sure the angles are in the range of 0-360
            curr_refine_ang_list[curr_refine_ang_list < 0] += 360
            curr_refine_ang_list[curr_refine_ang_list > 360] -= 360
            with tqdm(total=len(curr_refine_ang_list), position=1, leave=False) as pbar:
                for rot_ang in curr_refine_ang_list:
                    result = self._rot_and_search_fft_ss(rot_ang, False, self.vec_ave, self.vec_std, self.ss_ave, self.ss_std)
                    curr_result_list.append(
                        {
                            "angle": rot_ang,
                            "vec_score": result[0] / self.ref_map.new_data.size,
                            "vec_vox_trans": result[1],
                            "ss_score": result[2] / self.ref_map.new_data.size,
                            "ss_vox_trans": result[3],
                            "mix_score": result[4] / self.ref_map.new_data.size,
                            "vox_trans": result[5],
                            "real_trans": result[6],
                        }
                    )
                    pbar.update(1)
            self.refined_list.append(max(curr_result_list, key=lambda x: x["mix_score"]))
        self.refined_list.sort(key=lambda x: x["mix_score"], reverse=True)

    def _save_topn_pdb(self):
        os.makedirs(os.path.join(self.outdir, "PDB"), exist_ok=True)
        for i, item in enumerate(self.final_list):
            rot_mtx = R.from_euler("xyz", item["angle"], degrees=True).inv().as_matrix()
            angle_str = f"rx{int(item['angle'][0])}_ry{int(item['angle'][1])}_rz{int(item['angle'][2])}"
            trans_str = f"tx{item['real_trans'][0]:.3f}_ty{item['real_trans'][1]:.3f}_tz{item['real_trans'][2]:.3f}"
            filename = f"#{i}_{angle_str}_{trans_str}.pdb"
            save_rotated_pdb(self.input_pdb, rot_mtx, item["real_trans"], os.path.join(self.outdir, "PDB", filename), i)

    def _save_topn_vec_as_pdb(self):
        os.makedirs(os.path.join(self.outdir, "VEC"), exist_ok=True)
        for i, item in enumerate(self.final_list):
            angle_str = f"rx{int(item['angle'][0])}_ry{int(item['angle'][1])}_rz{int(item['angle'][2])}"
            trans_str = f"tx{item['real_trans'][0]:.3f}_ty{item['real_trans'][1]:.3f}_tz{item['real_trans'][2]:.3f}"
            filename = f"#{i}_{angle_str}_{trans_str}_VEC.pdb"
            save_vec_as_pdb(
                self.ref_map.new_orig,
                item["vec"],
                item["data"],
                item["score_arr"],
                item["score"],
                self.ref_map.new_width,
                item["vox_trans"],
                os.path.join(self.outdir, "VEC", filename),
                i,
            )

    def _save_topn_mrc(self):
        os.makedirs(os.path.join(self.outdir, "MRC"), exist_ok=True)
        for i, item in enumerate(self.final_list):
            angle_str = f"rx{int(item['angle'][0])}_ry{int(item['angle'][1])}_rz{int(item['angle'][2])}"
            trans_str = f"tx{item['real_trans'][0]:.3f}_ty{item['real_trans'][1]:.3f}_tz{item['real_trans'][2]:.3f}"
            filename = f"#{i}_{angle_str}_{trans_str}.mrc"
            save_rotated_mrc(
                self.tgt_map.mrcfile_path, item["angle"], item["real_trans"], os.path.join(self.outdir, "MRC", filename)
            )

    def _retrieve_data(self, result_list):
        for result in result_list:
            ret_dict = self._rot_and_search_fft(
                result["angle"],
                True,
            )
            result["data"] = ret_dict[2]
            result["vec"] = ret_dict[3]
            sco_arr, overlap, cc, pcc, Nm, total, dot = get_score(
                self.ref_map, result["data"], result["vec"], result["vox_trans"]
            )
            result["score_arr"] = sco_arr
            result["overlap"] = overlap
            result["cc"] = cc
            result["pcc"] = pcc
            result["Nm"] = Nm
            result["total"] = total

    @staticmethod
    def _calc_ldp_recall_score_item(ldp_arr, ca_arr, rot_mtx, trans):
        """
        Calculate the recall score of LDP points given a rotation matrix and translation vector
        All arguments have to be torch tensors on GPU
        # ldp_arr: torch tensor of shape (N, 3)
        # ca_arr: torch tensor of shape (N, 3)
        # rot_mtx: torch tensor of shape (3, 3)
        # trans: torch tensor of shape (3, )
        """
        import torch

        # rotated backbone CA
        rot_backbone_ca = torch.matmul(ca_arr, rot_mtx) + trans

        # calculate all pairwise distances
        dist_mtx = torch.cdist(rot_backbone_ca, ldp_arr, p=2)

        # get distance from the closest LDP point for each CA atom
        min_dist = torch.min(dist_mtx, dim=1).values

        # count the coverage of CA atoms within 3.0 angstrom of LDP points in the total amount of CA atoms
        return (min_dist < 3.0).sum().item() / len(rot_backbone_ca)

    def _calc_ldp_recall(self, results, sort=False, progress_bar=True):
        import torch

        if not progress_bar:
            iter_results = results  # calculate for each rotation
        else:
            iter_results = tqdm(results, desc="Calculating LDP Recall Score")
        for result in iter_results:
            r = R.from_euler("xyz", result["angle"], degrees=True)
            rot_mtx = (r.as_matrix()).T
            rot_mtx = torch.from_numpy(rot_mtx).to(self.device)
            # rot_mtx = euler_to_mtx(torch.tensor(result["angle"], device=self.device)).t()
            result["ldp_recall"] = self._calc_ldp_recall_score_item(
                self.ldp_atoms, self.backbone_ca, rot_mtx, torch.from_numpy(result["real_trans"]).to(self.device)
            )

        # sort by LDP recall
        if sort:
            results.sort(key=lambda x: x["ldp_recall"], reverse=True)

    def _remove_dup_results(self):
        no_dup_results = []

        print("###Start Duplicate Removal###")

        # duplicate removal
        hash_angs = {}

        # non_dup_count = 0

        # at least 30 degrees apart
        n_angles_apart = 30 // self.ang_interval  # could be directly specified
        ang_range = n_angles_apart * int(self.ang_interval)
        ang_range = int(ang_range)

        for result in tqdm(self.result_list, desc="Removing Duplicates"):
            # duplicate removal
            if tuple(result["angle"]) in hash_angs:
                # print(f"Duplicate: {result_mrc['angle']}")
                trans = hash_angs[tuple(result["angle"])]
                # manhattan distance
                if np.sum(np.abs(trans - result["vox_trans"])) < self.tgt_map.new_dim:
                    # result_mrc["vec_score"] = 0
                    continue

            # add to hash
            hash_angs[tuple(result["angle"])] = np.array(result["vox_trans"])

            ang_x, ang_y, ang_z = int(result["angle"][0]), int(result["angle"][1]), int(result["angle"][2])

            # add surrounding angles to hash
            for xx in range(ang_x - ang_range, ang_x + ang_range + 1, int(self.ang_interval)):
                for yy in range(ang_y - ang_range, ang_y + ang_range + 1, int(self.ang_interval)):
                    for zz in range(ang_z - ang_range, ang_z + ang_range + 1, int(self.ang_interval)):
                        x_positive = xx % 360
                        y_positive = yy % 360
                        z_positive = zz % 180

                        x_positive = x_positive + 360 if x_positive < 0 else x_positive
                        y_positive = y_positive + 360 if y_positive < 0 else y_positive
                        z_positive = z_positive + 180 if z_positive < 0 else z_positive

                        curr_trans = np.array([x_positive, y_positive, z_positive]).astype(np.float64)
                        # insert into hash
                        hash_angs[tuple(curr_trans)] = np.array(result["vox_trans"])

            # non_dup_count += 1
            no_dup_results.append(result)

        self.result_list = no_dup_results

    @staticmethod
    def _print_result_stats(results, return_stats=False):
        score_arr = np.array([item["score"] for item in results])
        ave = np.mean(score_arr)
        std = np.std(score_arr)
        print(f"Score Mean: {ave}, StDev: {std}")
        if return_stats:
            return ave, std

    def _print_result_item(self, item, i):

        print(
            "\n#" + str(i),
            "Rotation=",
            "(" + str(item["angle"][0]),
            str(item["angle"][1]),
            str(item["angle"][2]) + ")",
            "Translation=",
            "(" + "{:.3f}".format(item["real_trans"][0]),
            "{:.3f}".format(item["real_trans"][1]),
            "{:.3f}".format(item["real_trans"][2]) + ")",
        )

        print("Score=", "{:.6f}".format(item["score"]))
        norm_score = (item["score"] - self.score_ave) / self.score_std
        print(f"Voxel Trans= {item['vox_trans']}, Normalized Score= {norm_score:.6f}")

        if self.ldp_recall_mode:
            print("LDP Recall Score:", "{:.6f}".format(item["ldp_recall"]))

    def _convert_trans(self, rot_ang, vox_trans):

        r = R.from_euler("xyz", rot_ang, degrees=True)
        trans = np.array(vox_trans)

        if trans[0] > 0.5 * self.tgt_map.new_dim:
            trans[0] -= self.tgt_map.new_dim
        if trans[1] > 0.5 * self.tgt_map.new_dim:
            trans[1] -= self.tgt_map.new_dim
        if trans[2] > 0.5 * self.tgt_map.new_dim:
            trans[2] -= self.tgt_map.new_dim

        tgt_new_cent = r.apply(self.tgt_map.new_cent)  # rotate the center
        real_trans = self.ref_map.new_cent - (tgt_new_cent + trans * self.tgt_map.new_width)  # calculate new translation
        return real_trans

    def _rot_and_search_fft_ss(self, rot_ang, return_data, v_ave, v_std, ss_ave, ss_std):

        if self.gpu:
            import torch

            rot_mtx = euler_to_mtx(torch.tensor(np.radians(rot_ang), device=self.device, dtype=torch.float32))
            new_vec, new_data, new_ss_data = self._gpu_rot_map(
                self.tgt_map.new_data_gpu,
                self.tgt_map.vec_gpu,
                rot_mtx,
                self.tgt_map.search_pos_grid_gpu,
                self.device,
                rot_vec=True,
                ss_mix_score_mode=True,
                tgt_map_ss_data=self.tgt_map.new_ss_data_gpu,
            )
        else:
            rot_mtx = R.from_euler("xyz", rot_ang, degrees=True).as_matrix().astype(np.float32)
            new_vec, new_data, new_ss_data = self._rot_map(
                self.tgt_map.new_data,
                self.tgt_map.vec,
                rot_mtx,
                self.search_pos_grid,
                rot_vec=True,
                ss_mix_score_mode=True,
                tgt_map_ss_data=self.tgt_map.new_ss_data,
            )

        x2, y2, z2 = new_vec[..., 0], new_vec[..., 1], new_vec[..., 2]
        tgt_map_pre_fft_list = [x2, y2, z2]
        tgt_map_pre_fft_list.extend([new_ss_data[..., i] for i in range(4)])

        if self.gpu:
            fft_result_list = self._fft_get_prod_list(self.ref_map_fft_list_gpu, tgt_map_pre_fft_list)
        else:
            fft_result_list = self._fft_get_prod_list(self.ref_map_fft_list, tgt_map_pre_fft_list)

        # convert back to numpy if needed
        if self.gpu and return_data:
            new_data = new_data.cpu().numpy()
            new_vec = new_vec.cpu().numpy()
            new_ss_data = new_ss_data.cpu().numpy()

        # Search for best translation using FFT
        vec_score, vec_vox_trans = self._find_best_trans_by_fft_list(fft_result_list[:3], self.gpu)
        ss_score, ss_vox_trans = self._find_best_trans_by_fft_list(fft_result_list[3:-1], self.gpu)

        if v_ave is not None and v_std is not None and ss_ave is not None and ss_std is not None:
            mix_score, mix_vox_trans = self._find_best_trans_by_fft_list_ss(
                fft_result_list, self.alpha, v_ave, v_std, ss_ave, ss_std, self.gpu
            )
        else:
            mix_score, mix_vox_trans = None, None

        if mix_vox_trans is not None:
            mix_real_trans = self._convert_trans(rot_ang, mix_vox_trans)
        else:
            mix_real_trans = None

        if return_data:
            return (
                vec_score,
                vec_vox_trans,
                ss_score,
                ss_vox_trans,
                mix_score,
                mix_vox_trans,
                mix_real_trans,
                new_data,
                new_vec,
                new_ss_data,
            )
        else:
            return vec_score, vec_vox_trans, ss_score, ss_vox_trans, mix_score, mix_vox_trans, mix_real_trans

    def _rot_and_search_fft(self, rot_ang, return_data):
        if self.gpu:
            import torch

            rot_mtx = R.from_euler("xyz", rot_ang, degrees=True).as_matrix().astype(np.float32)
            rot_mtx = torch.from_numpy(rot_mtx).to(self.device)
            # rot_mtx = euler_to_mtx(torch.tensor(np.radians(rot_ang), device=self.device, dtype=torch.float32))
            new_vec, new_data, _ = self._gpu_rot_map(
                self.tgt_map.new_data_gpu,
                self.tgt_map.vec_gpu,
                rot_mtx,
                self.tgt_map.search_pos_grid_gpu,
                self.device,
                rot_vec=(self.mode == "VecProduct"),
            )
        else:
            rot_mtx = R.from_euler("xyz", rot_ang, degrees=True).as_matrix().astype(np.float32)
            new_vec, new_data, _ = self._rot_map(
                self.tgt_map.new_data,
                self.tgt_map.vec,
                rot_mtx,
                self.search_pos_grid,
                rot_vec=(self.mode == "VecProduct"),
            )

        if self.mode == "VecProduct":
            # compose query list for vector product mode (x, y, z)
            x2, y2, z2 = new_vec[..., 0], new_vec[..., 1], new_vec[..., 2]
            tgt_map_pre_fft_list = [x2, y2, z2]
        else:
            # 1D mode
            if not self.gpu:
                x2 = new_data
                if self.mode == "Overlap":
                    x2 = np.where(x2 > 0, 1.0, 0.0)
                elif self.mode == "CC":
                    x2 = np.where(x2 > 0, x2, 0.0)
                elif self.mode == "PCC":
                    x2 = np.where(x2 > 0, x2 - self.tgt_map.ave, 0.0)
                elif self.mode == "Laplacian":
                    x2 = laplace(x2, mode="constant", cval=0.0)
                tgt_map_pre_fft_list = [x2]
            else:
                x2 = new_data
                zero_tensor = torch.tensor([0], device=self.device, dtype=torch.float32)
                one_tensor = torch.tensor([1], device=self.device, dtype=torch.float32)
                if self.mode == "Overlap":
                    x2 = torch.where(x2 > zero_tensor, one_tensor, zero_tensor)  # binaries
                elif self.mode == "CC":
                    x2 = torch.where(x2 > zero_tensor, x2, zero_tensor)  # zero mean
                elif self.mode == "PCC":
                    x2 = torch.where(x2 > zero_tensor, x2 - self.tgt_map.ave, zero_tensor)  # zero mean
                elif self.mode == "Laplacian":
                    x2 = laplace(x2.cpu().numpy(), mode="constant", cval=0.0)
                    x2 = torch.from_numpy(x2).to(self.device)
                tgt_map_pre_fft_list = [x2]

        if self.gpu:
            fft_result_list = self._fft_get_prod_list(self.ref_map_fft_list_gpu, tgt_map_pre_fft_list)
        else:
            fft_result_list = self._fft_get_prod_list(self.ref_map_fft_list, tgt_map_pre_fft_list)

        # convert back to numpy if needed
        if self.gpu and return_data:
            new_data = new_data.cpu().numpy()
            if new_vec is not None:
                new_vec = new_vec.cpu().numpy()

        # Search for best translation using FFT
        score, vox_trans = self._find_best_trans_by_fft_list(fft_result_list, gpu=self.gpu)

        if self.mode == "CC":
            score = score / (self.ref_map.std**2)
        if self.mode == "PCC":
            score = score / (self.ref_map.std_norm_ave**2)

        # return data if specified
        if return_data:
            return score, vox_trans, new_data, new_vec
        else:
            return score, vox_trans

    @staticmethod
    def _gpu_rot_map(data, vec, mtx, new_pos_grid, device, rot_vec=True, ss_mix_score_mode=False, tgt_map_ss_data=None):
        import torch

        with torch.no_grad():
            # set the dimension to be x dimension as all dimension are the same
            dim = data.shape[0]

            # set the rotation center
            cent = 0.5 * float(dim)
            cent = torch.tensor(cent, device=device, dtype=torch.float32, requires_grad=False)

            # get relative new positions from center
            new_pos = new_pos_grid - cent

            # reversely rotate the new position lists to get old positions
            # old_pos = torch.einsum("ij, kj->ki", mtx.T, new_pos) + cent
            old_pos = new_pos @ mtx + cent

            # round old positions to nearest integer
            old_pos.round_()

            # init new vec and dens array
            new_data_array = torch.zeros_like(data, device=device, dtype=torch.float32, requires_grad=False)

            in_bound_mask = torch.all((old_pos >= 0) & (old_pos < dim), axis=1)

            # get valid old positions in bound
            valid_old_pos = (old_pos[in_bound_mask]).long()

            # get nonzero density positions in the map
            # non_zero_mask = data[valid_old_pos[:, 0], valid_old_pos[:, 1], valid_old_pos[:, 2]] > 0

            # apply nonzero mask to valid positions
            # non_zero_old_pos = valid_old_pos[non_zero_mask]

            # get corresponding new positions
            # new_pos = (new_pos[in_bound_mask][non_zero_mask] + cent).long()
            new_pos = new_pos[in_bound_mask].add_(cent).long()

            # fill new density entries
            new_data_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = data[
                valid_old_pos[:, 0], valid_old_pos[:, 1], valid_old_pos[:, 2]
            ]

            if ss_mix_score_mode:
                new_ss_array = torch.zeros_like(tgt_map_ss_data, device=device, dtype=torch.float32, requires_grad=False)
                new_ss_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = tgt_map_ss_data[
                    valid_old_pos[:, 0], valid_old_pos[:, 1], valid_old_pos[:, 2]
                ]
            else:
                new_ss_array = None

            if rot_vec:
                new_vec_array = torch.zeros_like(vec, device=device, dtype=torch.float32, requires_grad=False)
                # fetch and rotate the vectors
                non_zero_vecs = vec[valid_old_pos[:, 0], valid_old_pos[:, 1], valid_old_pos[:, 2]]

                new_vec = non_zero_vecs @ mtx.T
                # new_vec = torch.einsum("ij, kj->ki", mtx, non_zero_vecs)

                # fill new vector entries
                new_vec_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = new_vec
            else:
                new_vec_array = None

        return new_vec_array, new_data_array, new_ss_array

    @staticmethod
    def _rot_map(data, vec, mtx, new_pos_grid, rot_vec=True, ss_mix_score_mode=False, tgt_map_ss_data=None):
        # set the dimension to be x dimension as all dimension are the same
        dim = data.shape[0]

        # set the rotation center
        cent = 0.5 * float(dim)

        # get relative new positions from center
        new_pos = new_pos_grid - cent

        # reversely rotate the new position lists to get old positions
        old_pos = np.einsum("ij, kj->ki", mtx.T, new_pos) + cent
        # old_pos = new_pos @ mtx + 0.5 * float(dim)
        old_pos = np.round(old_pos)  # round old positions to nearest integer (nearest voxel)

        # init new vec and dens array
        new_data_array = np.zeros_like(data)

        in_bound_mask = np.all((old_pos >= 0) & (old_pos < dim), axis=1)

        # get valid old positions in bound
        valid_old_pos = (old_pos[in_bound_mask]).astype(int)

        # get nonzero density positions in the map
        non_zero_mask = data[valid_old_pos[:, 0], valid_old_pos[:, 1], valid_old_pos[:, 2]] > 0

        # apply nonzero mask to valid positions
        non_zero_old_pos = valid_old_pos[non_zero_mask]

        # get corresponding new positions
        new_pos = (new_pos[in_bound_mask][non_zero_mask] + cent).astype(np.int32)

        # fill new density entries
        new_data_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = data[
            non_zero_old_pos[:, 0], non_zero_old_pos[:, 1], non_zero_old_pos[:, 2]
        ]

        if rot_vec:
            new_vec_array = np.zeros_like(vec)
            # fetch and rotate the vectors
            non_zero_vecs = vec[non_zero_old_pos[:, 0], non_zero_old_pos[:, 1], non_zero_old_pos[:, 2]]

            # new_vec = non_zero_vecs @ mtx.T
            new_vec = np.einsum("ij, kj->ki", mtx, non_zero_vecs)

            # fill new vector entries
            new_vec_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = new_vec
        else:
            new_vec_array = None

        if ss_mix_score_mode:
            new_ss_array = np.zeros_like(tgt_map_ss_data)
            new_ss_array[new_pos[:, 0], new_pos[:, 1], new_pos[:, 2]] = tgt_map_ss_data[
                non_zero_old_pos[:, 0], non_zero_old_pos[:, 1], non_zero_old_pos[:, 2]
            ]
        else:
            new_ss_array = None

        return new_vec_array, new_data_array, new_ss_array

    def _fft_get_prod_list(self, ref_map_fft_list, tgt_map_fft_list):
        dot_product_list = []
        if self.gpu:
            import torch

            for ref_fourier, tgt_real in zip(ref_map_fft_list, tgt_map_fft_list):
                tgt_fourier = torch.fft.rfftn(tgt_real)
                dot_fourier = ref_fourier * tgt_fourier
                dot_real = torch.fft.irfftn(dot_fourier, norm="ortho")
                dot_product_list.append(dot_real)
        else:
            from pyfftw.interfaces import numpy_fft

            for ref_fourier, tgt_real in zip(ref_map_fft_list, tgt_map_fft_list):
                tgt_fourier = numpy_fft.rfftn(tgt_real)
                dot_fourier = ref_fourier * tgt_fourier
                dot_real = numpy_fft.irfftn(dot_fourier, norm="ortho")
                dot_product_list.append(dot_real)

        return dot_product_list

    @staticmethod
    def _find_best_trans_by_fft_list(fft_result_list, gpu=False):
        if gpu:
            import torch

            sum_arr = torch.zeros_like(fft_result_list[0])
            for item in fft_result_list:
                sum_arr += item
            best_score = torch.amax(sum_arr).cpu().numpy()
            trans = torch.argmax(sum_arr).cpu().numpy()
            trans = np.unravel_index(trans, sum_arr.shape)
            return best_score, trans
        else:
            sum_arr = np.sum(fft_result_list, axis=0)
            best_score = np.amax(sum_arr)
            trans = np.unravel_index(sum_arr.argmax(), sum_arr.shape)
            return best_score, trans

    @staticmethod
    def _find_best_trans_by_fft_list_ss(
        fft_result_list, alpha, vec_score_mean, vec_score_std, ss_score_mean, ss_score_std, gpu=False
    ):
        if gpu:
            sum_arr_v = torch.stack(fft_result_list[:3]).sum(dim=0)
            sum_arr_ss = torch.stack(fft_result_list[4:-1]).sum(dim=0)  # do not include nucleotide score

            # z-score normalization
            sum_arr_v = (sum_arr_v - vec_score_mean) / vec_score_std
            sum_arr_ss = (sum_arr_ss - ss_score_mean) / ss_score_std

            sum_arr_mixed = (1 - alpha) * sum_arr_v + alpha * sum_arr_ss
            best_score = torch.amax(sum_arr_mixed).cpu().numpy()
            best_trans = np.unravel_index(sum_arr_mixed.cpu().numpy().argmax(), sum_arr_mixed.shape)

        else:
            sum_arr_v = np.sum(fft_result_list[:3], axis=0)
            sum_arr_ss = np.sum(fft_result_list[3:-1], axis=0)  # do not include nucleotide score

            # z-score normalization
            sum_arr_v = (sum_arr_v - vec_score_mean) / vec_score_std
            sum_arr_ss = (sum_arr_ss - ss_score_mean) / ss_score_std

            sum_arr_mixed = (1 - alpha) * sum_arr_v + alpha * sum_arr_ss
            best_score = sum_arr_mixed.max()
            best_trans = np.unravel_index(sum_arr_mixed.argmax(), sum_arr_mixed.shape)

        return best_score, best_trans

    def _calc_angle_comb(self):
        """Calculate the all the possible combination of angles given the interval in degrees"""

        if self.confine_angles is not None:
            x_angle = y_angle = z_angle = np.arange(-self.confine_angles, self.confine_angles + 1, self.ang_interval)
        else:
            xy_limit = 360
            z_limit = 180

            x_angle = np.arange(0, xy_limit, self.ang_interval)
            y_angle = np.arange(0, xy_limit, self.ang_interval)
            z_angle = np.arange(0, z_limit + 1, self.ang_interval)

        # make sure positive angles are in the range of 0-360
        x_angle[x_angle < 0] += 360
        y_angle[y_angle < 0] += 360
        z_angle[z_angle < 0] += 180

        angle_comb = np.array(np.meshgrid(x_angle, y_angle, z_angle)).T.reshape(-1, 3)

        # compose angle_comb
        seen = set()
        for ang in angle_comb:
            r = R.from_euler("xyz", ang, degrees=True)
            quat = tuple(np.round(r.as_quat(), 4))
            if quat not in seen:
                seen.add(quat)
                self.angle_comb.append(ang)
