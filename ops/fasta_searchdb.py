import os
from ops.os_operation import mkdir, functToDeleteItems
from ops.io_utils import download_file
from ops.pdb_utils import count_atom_line, filter_chain_cif, cif2pdb, filter_chain_pdb, count_residues
from ops.fasta_utils import read_fasta, write_all_fasta
from ops.io_utils import write_pickle, load_pickle
from collections import defaultdict
import time
import requests


def parse_blast_output(blast_file):
    match_dict = defaultdict(list)
    read_flag = False
    with open(blast_file, 'r') as rfile:
        all_lines = rfile.readlines()
        for line_index, line in enumerate(all_lines):
            if line.startswith("Query="):
                line = line.strip("\n")
                line = line.replace("Query=", "")
                current_id = line.replace(" ", "")
                # check if following line still have some thing
                next_line_index = line_index + 1
                while len(all_lines[next_line_index].strip("\n").replace(" ", "")) != 0:
                    new_line = all_lines[next_line_index].strip("\n").replace(" ", "")
                    current_id += new_line
                    next_line_index += 1
                print("parse blast id:", current_id)
            elif line.startswith("Sequences producing significant"):
                read_flag = True
                continue
            if read_flag:
                line = line.strip("\n")
                if len(line) > 0:
                    line = line.split()
                    match_id = line[0]
                    score = float(line[-1])
                    match_dict[current_id].append([match_id, score])
                    for k in range(1, 500):
                        new_line = all_lines[k + line_index]
                        new_line = new_line.strip("\n").split()
                        try:
                            match_id = new_line[0]
                            score = float(new_line[-1])
                            match_dict[current_id].append([match_id, score])
                        except:
                            break
                    read_flag = False
    return match_dict


def download_pdb(pdb_id, current_chain_dir, final_pdb_path):
    pdb = pdb_id.split("_")[0]
    chain_id = pdb_id.split("_")[1]
    download_link = "https://files.rcsb.org/download/%s.cif" % pdb
    cif_file = os.path.join(current_chain_dir, "input.cif")
    download_file(download_link, cif_file)
    if os.path.exists(cif_file) and count_atom_line(cif_file) >= 50:
        # segment the specific chains
        chain_cif = os.path.join(current_chain_dir, "input_%s.cif" % chain_id)
        filter_chain_cif(cif_file, chain_id, chain_cif)

        # then convert cif file format to pdb
        cif2pdb(chain_cif, final_pdb_path)
    else:
        # download the pdb if the old one did not exist
        download_link = "https://files.rcsb.org/download/%s.pdb" % pdb
        pdb_file = os.path.join(current_chain_dir, "input.pdb")
        download_file(download_link, pdb_file)
        filter_chain_pdb(pdb_file, chain_id, final_pdb_path)
    return final_pdb_path


def get_metadata(pdb_id, max_retry, type):
    if type == "pdb":
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    elif type == "afdb":
        url = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{pdb_id}.json"
    while max_retry > 0:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response.json()
        except:
            max_retry -= 1  # try again
            time.sleep(30)
        finally:
            max_retry -= 1  # try again after 30 seconds
    return None  # explicitly return None


def fasta_searchdb(params, save_path):
    max_file_name_limit = 250  # default is 255, to be safe
    single_chain_pdb_dir = os.path.join(save_path, "single_chain_pdb")
    mkdir(single_chain_pdb_dir)
    fitting_pickle_path = os.path.join(save_path, "fitting_dict.pkl")
    if os.path.exists(fitting_pickle_path) and os.path.getsize(fitting_pickle_path) > 10:
        fitting_dict = load_pickle(fitting_pickle_path)
        return fitting_dict
    # reorganize data
    fasta_path = os.path.abspath(params['P'])
    chain_dict = read_fasta(fasta_path)
    fasta_path = os.path.join(single_chain_pdb_dir, "input.fasta")
    write_all_fasta(chain_dict, fasta_path)

    # first do experimental database search
    if not params['af_only']:
        output_path = os.path.join(single_chain_pdb_dir, "exp_search.out")
        search_command = "blastp -query %s -db %s -out %s -num_threads %d" % (fasta_path, params['db_exp_path'],
                                                                              output_path, params['search_thread'])
        if not os.path.exists(output_path) or os.path.getsize(output_path) < 1000:
            os.system(search_command)
        exp_match_dict = parse_blast_output(output_path)
    else:
        exp_match_dict = {}  # empty dict
    matched_dict = {}  # [key]: chain id list, [value]: the structure id
    fitting_dict = {}
    # parse the information

    # merge only identical chains, otherwise, rely on combine search
    for key in exp_match_dict:
        current_match_list = exp_match_dict[key]
        print("match candidate:", key, current_match_list)
        for k in range(len(current_match_list)):
            match_id, evalue = current_match_list[k]
            if evalue == 0:
                # fetch the pdb to see if the protein size really reasonable
                expected_seq_length = len(chain_dict[key]) * params['search']['length_ratio']
                chain_name_list = key.replace(",", "-")
                current_chain_dir = os.path.join(single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit])
                mkdir(current_chain_dir)
                final_pdb_path = os.path.join(single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb")
                download_pdb(match_id, current_chain_dir, final_pdb_path)
                actual_structure_length = count_residues(final_pdb_path)
                if actual_structure_length >= expected_seq_length and actual_structure_length <= len(chain_dict[key]):
                    matched_dict[key] = "PDB:" + match_id
                    final_chain_list = chain_name_list.split("-")
                    fitting_dict[final_pdb_path] = final_chain_list
                    break
                else:
                    os.remove(final_pdb_path)
                    functToDeleteItems(current_chain_dir)
            if evalue > 0:
                break
    matched_keys = matched_dict.keys()
    # write a new fasta to search

    remain_chain_dict = {k: v for k, v in chain_dict.items() if k not in matched_keys}
    print("experimental db search finished, get %s/%s matched single chain structure" % (
        len(matched_keys), len(chain_dict)))
    if len(remain_chain_dict) > 0:
        print("continue PDB+AFDB search for remained %d chains" % len(remain_chain_dict))
        remain_fasta_path = os.path.join(single_chain_pdb_dir, "remain_search.fasta")
        write_all_fasta(remain_chain_dict, remain_fasta_path)
        output_path = os.path.join(single_chain_pdb_dir, "expaf_search.out")
        search_command = "blastp -query %s -db %s -out %s -num_threads %d" % (remain_fasta_path, params['db_path'],
                                                                              output_path, params['search_thread'])
        if not os.path.exists(output_path) or os.path.getsize(output_path) < 1000:
            os.system(search_command)
        expaf_match_dict = parse_blast_output(output_path)
        if len(expaf_match_dict) != len(remain_chain_dict):
            for tmp_key in remain_chain_dict:
                if tmp_key not in expaf_match_dict:
                    print("warning, query %s can not find any templates in database" % tmp_key)

        for key in expaf_match_dict:
            closest_choice = None
            num_res_difference = 999999999
            current_match_list = expaf_match_dict[key]
            print("match candidate:", key, current_match_list)
            chain_name_list = key.replace(",", "-")
            current_chain_dir = os.path.join(single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit])
            mkdir(current_chain_dir)
            final_pdb_path = os.path.join(single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb")
            expected_seq_length = len(chain_dict[key]) * params['search']['length_ratio']
            for k in range(len(current_match_list)):
                match_id, evalue = current_match_list[k]
                if params['af_only']:
                    if "AFDB" not in match_id:
                        continue
                if "AFDB" in match_id:
                    split_info = match_id.split(":")
                    database = split_info[0]
                    pdb_id = split_info[1]

                    metadata = get_metadata(pdb_id.split("-")[1], 3, "afdb")

                    try:
                        actual_structure_length = metadata["uniprot_entry"]["sequence_length"]
                    except:
                        print("get metadata failed for %s" % pdb_id)
                        actual_structure_length = 0
                else:
                    metadata = get_metadata(match_id.split("_")[0], 3, "pdb")
                    try:
                        actual_structure_length = metadata["rcsb_entry_info"]["deposited_polymer_monomer_count"]
                    except:
                        print("get metadata failed for %s" % match_id)
                        actual_structure_length = 0

                if k == 0:
                    # save the top 1 in case of no length match
                    top1_structure_length = actual_structure_length

                curr_seq_length_diff = abs(expected_seq_length - actual_structure_length)

                if expected_seq_length <= actual_structure_length <= len(chain_dict[key]):
                    if "AFDB" in match_id:
                        matched_dict[key] = match_id
                    else:
                        matched_dict[key] = "PDB:" + match_id
                    final_chain_list = chain_name_list.split("-")
                    fitting_dict[final_pdb_path] = final_chain_list
                    closest_choice = match_id
                    break  # already in match dict, no need to check anymore

                if curr_seq_length_diff < num_res_difference:
                    num_res_difference = curr_seq_length_diff
                    closest_choice = match_id

            # print("Closest choice is %s with %d difference" % (closest_choice))

            if key not in matched_dict:  # no match under length condition
                if params['af_only'] and closest_choice is None:
                    # we can only pick the top 1 candidate
                    match_id, evalue = current_match_list[0]
                    closest_choice = match_id
                    num_res_difference = abs(expected_seq_length - top1_structure_length)

                if "AFDB" in closest_choice:
                    matched_dict[key] = closest_choice
                else:
                    matched_dict[key] = "PDB:" + closest_choice
                print("we have no better choice but pick %s with %d residues differences" % (
                    closest_choice, num_res_difference))

            # download the closest choice as final
            if "AFDB" in closest_choice:
                download_link = "https://alphafold.ebi.ac.uk/files/%s-model_v4.pdb" % closest_choice.split(":")[1]
                download_flag = download_file(download_link, final_pdb_path)
                if download_flag is False:
                    time.sleep(60)
                    continue
            else:
                download_pdb(closest_choice, current_chain_dir, final_pdb_path)

    print("DB search finished! Match relationship ", matched_dict)
    # get the matched dict

    for chain_name_list in chain_dict:
        if chain_name_list not in matched_dict:
            print("*" * 100)
            print("WARNING! query %s can not find any templates in database" % chain_name_list)
            print("*" * 100)
            continue
        matched_id = matched_dict[chain_name_list]
        chain_name_list = chain_name_list.replace(",", "-")
        current_chain_dir = os.path.join(single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit])
        mkdir(current_chain_dir)
        final_pdb_path = os.path.join(single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb")
        if final_pdb_path in fitting_dict:
            continue
        if os.path.exists(final_pdb_path) and count_atom_line(final_pdb_path) >= 50:
            final_chain_list = chain_name_list.split("-")
            fitting_dict[final_pdb_path] = final_chain_list
            continue
        split_info = matched_id.split(":")
        database = split_info[0]
        pdb_id = split_info[1]
        if database == "AFDB":
            # alphafold db
            download_link = "https://alphafold.ebi.ac.uk/files/%s-model_v4.pdb" % pdb_id
            download_flag = download_file(download_link, final_pdb_path)
            if download_flag is False:
                time.sleep(60)
                continue
        else:
            download_pdb(pdb_id, current_chain_dir, final_pdb_path)
        final_chain_list = chain_name_list.split("-")
        fitting_dict[final_pdb_path] = final_chain_list
    print("collecting finish: fitting dict: ", fitting_dict)
    write_pickle(fitting_dict, fitting_pickle_path)
    return fitting_dict
