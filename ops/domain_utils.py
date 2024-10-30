import os
import pickle
import shutil

# from ops.fasta_to_similar_pdb import print_usage
from ops.io_utils import load_json
from ops.pdb_utils import filter_pdb_by_residues


# call sword2 to get different domains and then run diffmodeler with domains
def prepare_domain_input(fitting_dict, save_dir, num_cpu=8):
    """
    :param fitting_dict:    [pdb_path]:[list of chain_id]
    :param save_dir:
    return:
    fitting_dict
    """
    domain_fitting_dict_path = os.path.join(save_dir, "domain_fitting_dict.pkl")
    if os.path.exists(domain_fitting_dict_path):
        fitting_dict = pickle.load(open(domain_fitting_dict_path, "rb"))
        # check every path exists
        verify_fitting_dict = True
        for pdb_path in fitting_dict:
            if not os.path.exists(pdb_path):
                verify_fitting_dict = False
                break
        if verify_fitting_dict:
            return fitting_dict

    current_dir = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.path.join(current_dir, "SWORD2")
    sword_script_path = os.path.join(current_dir, "SWORD2.py")
    # clean the save_dir
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    domain_fitting_dict = {}

    # clean the binaries
    dssp_bin = os.path.join(current_dir, "bin/SWORD/bin/SWORD/bin/Dssp/dsspcmbi")
    peeling_bin = os.path.join(current_dir, "bin/SWORD/bin/SWORD/bin/Peeling_omp")
    scoring_bin = os.path.join(current_dir, "bin/mypmfs-master/scoring_omp")
    if os.path.exists(dssp_bin):
        os.remove(dssp_bin)
    if os.path.exists(peeling_bin):
        os.remove(peeling_bin)
    if os.path.exists(scoring_bin):
        os.remove(scoring_bin)

    # compile the SWORD2 C++ binaries
    os.environ["MKL_ENABLE_INSTRUCTIONS"] = "SSE4_2"
    os.system(f"{os.path.join(current_dir, 'bin/SWORD/bin/SWORD/SWORD')} > /dev/null")
    os.system(f"make -C {os.path.join(current_dir, 'bin/mypmfs-master/')} scoring_omp > /dev/null")
    os.system(f"make -C {os.path.join(current_dir, 'bin/SWORD/bin/SWORD/bin')} Peeling_omp > /dev/null")

    for pdb_path in fitting_dict:
        chain_ids = fitting_dict[pdb_path]
        pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
        pdb_dir = os.path.join(save_dir, pdb_name)
        os.makedirs(pdb_dir, exist_ok=True)
        print("Launching SWORD2 for %s" % pdb_path)
        if num_cpu > 0:
            # split pdb into different domains
            command_line = f"python {sword_script_path} -i {pdb_path} -o {pdb_dir} -c A --cpu {num_cpu} --disable-energies --disable-plots"
        else:
            # 0 indicates use all cpu available
            command_line = f"python {sword_script_path} -i {pdb_path} -o {pdb_dir} -c A --disable-energies --disable-plots"
        os.system(command_line)
        cur_output_dir = os.path.join(pdb_dir, f"{pdb_name}_A")
        if not os.path.exists(cur_output_dir):
            print(f"Error: SWORD2 failed to split {pdb_path}, output dir does not exist.")
            continue
        cur_output_json = os.path.join(cur_output_dir, "SWORD2_summary.json")
        if not os.path.exists(cur_output_json):
            print(f"Error: SWORD2 failed to split {pdb_path}, summary json does not exist.")
            domain_pdb_path = os.path.join(save_dir, f"{pdb_name}.pdb")
            shutil.copy(pdb_path, domain_pdb_path)
            domain_fitting_dict[domain_pdb_path] = chain_ids
            continue
        try:
            domain_info = load_json(cur_output_json)
            boundary_info = domain_info["Optimal partition"]["Domains"]
            print(f"Domain info for {pdb_path}: {boundary_info}")
        except:
            print(f"Error: SWORD2 failed to split {pdb_path}, boundary info missing")
            domain_pdb_path = os.path.join(save_dir, f"{pdb_name}.pdb")
            shutil.copy(pdb_path, domain_pdb_path)
            domain_fitting_dict[domain_pdb_path] = chain_ids
            continue
        for domain_name, current_domain in boundary_info.items():
            # filter all residues
            domain_name = domain_name.replace(" ", "_")
            residue_index_list = []
            dom = current_domain.get("PUs", {})
            for tmp_range in dom:
                res_range = tmp_range.split("-")
                residue_index_list.extend(list(range(int(res_range[0]), int(res_range[1]) + 1)))
            print("Domain name:", domain_name, "Residue indices:", residue_index_list)
            # output the domain pdb
            #
            pdb_path_renumed = os.path.join(cur_output_dir, f"{pdb_name}_A")
            print("Renumbered PDB path:", pdb_path_renumed)
            domain_pdb_path = os.path.join(save_dir, f"{pdb_name}_{domain_name}.pdb")
            filter_pdb_by_residues(pdb_path_renumed, domain_pdb_path, residue_index_list)
            current_chain_ids = [x + "_" + str(domain_name) for x in chain_ids]
            domain_fitting_dict[domain_pdb_path] = current_chain_ids
            print(f"Domain {domain_name} for {pdb_path} is saved at {domain_pdb_path}")
    pickle.dump(domain_fitting_dict, open(domain_fitting_dict_path, "wb"))

    return domain_fitting_dict
