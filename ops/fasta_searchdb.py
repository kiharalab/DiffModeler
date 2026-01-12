from collections import defaultdict
import os
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

from ops.os_operation import mkdir, functToDeleteItems
from ops.pdb_utils import (
    count_atom_line,
    filter_chain_cif,
    cif2pdb,
    filter_chain_pdb,
    count_residues,
)
from ops.fasta_utils import read_fasta, write_all_fasta
from ops.io_utils import write_pickle, load_pickle, download_file


def get_afdb_pdb_url(model_id: str, max_retry: int = 2):
    """
    get available download url from AFDB api.
    input should be a form of 'AF-A0A1B0GTW7-F1'

    Returns the pdbUrl if found, None otherwise.
    """
    # uniprot ID format
    # 1. XXXXXX (single string/number)
    # 2. XXXXXX-1 (with isoform ID)
    # API is only searchable without isoform ID, but it returns all isoforms.

    uniprot_id = model_id.split("-")[1]
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    while max_retry > 0:
        try:
            r = requests.get(url)
            if r.status_code == 200:
                break
            elif r.status_code == 404:
                print(f"AlphaFold prediction not found (404) for UniProt ID: {uniprot_id}")
                return None
        except:
            max_retry = max_retry - 1
            print(f"request failed, retring in 5 seconds... {url}")
            time.sleep(5)

    response_data = r.json()

    # Check if response is empty
    if not response_data:
        print(f"No AlphaFold prediction found for UniProt ID: {uniprot_id} (model: {model_id})")
        return None

    # Try to find exact match first
    for d in response_data:
        if d["modelEntityId"] == model_id:
            return d["pdbUrl"]

    # Fallback to first entry if no exact match
    return response_data[0]["pdbUrl"]


def parse_blast_output(blast_file):
    match_dict = defaultdict(list)
    read_flag = False
    with open(blast_file, "r") as rfile:
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


def download_pdb(pdb_id, current_chain_dir, final_pdb_path, first_model_only=True):
    pdb = pdb_id.split("_")[0]
    chain_id = pdb_id.split("_")[1]
    download_link = "https://files.rcsb.org/download/%s.cif" % pdb
    cif_file = os.path.join(current_chain_dir, "input.cif")
    download_file(download_link, cif_file)
    if os.path.exists(cif_file) and count_atom_line(cif_file) >= 50:
        # segment the specific chains
        chain_cif = os.path.join(current_chain_dir, "input_%s.cif" % chain_id)
        filter_chain_cif(cif_file, chain_id, chain_cif, first_model_only=first_model_only)

        # then convert cif file format to pdb
        cif2pdb(chain_cif, final_pdb_path)
    else:
        # download the pdb if the old one did not exist
        download_link = "https://files.rcsb.org/download/%s.pdb" % pdb
        pdb_file = os.path.join(current_chain_dir, "input.pdb")
        download_file(download_link, pdb_file)
        filter_chain_pdb(pdb_file, chain_id, final_pdb_path, first_model_only=first_model_only)
    return final_pdb_path


def get_metadata(pdb_id, max_retry, type, chain_id=None):
    if type == "pdb":
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/{pdb_id}/chain/{chain_id}"
    elif type == "afdb":
        url = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{pdb_id}.json"
    # print("fetching metadata from %s" % url)
    while max_retry > 0:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response.json()
            if response.status_code == 404:
                print(f"Fetching metadata failed (404): {url} for {pdb_id} {chain_id}")
                return None
        except:
            max_retry -= 1  # try again
            print("Fetching metadata failed, retrying after 5 seconds %s" % url)
            time.sleep(5)
        finally:
            max_retry -= 1  # try again after 30 seconds
    return None  # explicitly return None


def fetch_afdb_candidate(match_id, pdb_id, uniprot_id):
    """Fetch metadata and download URL for an AFDB candidate."""
    try:
        # Prefetch metadata
        metadata = get_metadata(uniprot_id, 3, "afdb")
        if metadata is None:
            return None
        
        # Prefetch download URL
        download_link = get_afdb_pdb_url(pdb_id)
        if download_link is None:
            return None
        
        try:
            actual_structure_length = metadata["uniprot_entry"]["sequence_length"]
        except:
            return None
        
        return {
            "match_id": match_id,
            "actual_structure_length": actual_structure_length,
            "download_link": download_link,
            "metadata": metadata
        }
    except Exception as e:
        print(f"Error fetching AFDB candidate {pdb_id}: {e}")
        return None


def fetch_pdb_candidate(match_id):
    """Fetch metadata for a PDB candidate."""
    try:
        pdb_id = match_id.split("_")[0]
        chain_id = match_id.split("_")[1]
        metadata = get_metadata(pdb_id, 3, "pdb", chain_id)
        if metadata is None:
            return None
        
        try:
            actual_structure_length = 0
            segments = metadata[pdb_id.lower()]["molecules"][0]["chains"][0]["observed"]
            for segment in segments:
                actual_structure_length += (
                    segment["end"]["residue_number"]
                    - segment["start"]["residue_number"]
                    + 1
                )
        except:
            return None
        
        return {
            "match_id": match_id,
            "actual_structure_length": actual_structure_length,
            "download_link": None,  # PDB downloads handled separately
            "metadata": metadata
        }
    except Exception as e:
        print(f"Error fetching PDB candidate {match_id}: {e}")
        return None


def fasta_searchdb(params, save_path):
    max_file_name_limit = 250  # default is 255, to be safe
    single_chain_pdb_dir = os.path.join(save_path, "single_chain_pdb")
    mkdir(single_chain_pdb_dir)
    fitting_pickle_path = os.path.join(save_path, "fitting_dict.pkl")
    if (
        os.path.exists(fitting_pickle_path)
        and os.path.getsize(fitting_pickle_path) > 10
    ):
        fitting_dict = load_pickle(fitting_pickle_path)
        # check all the path exists, otherwise do not allow skip
        use_flag = True
        for key in fitting_dict:
            if not os.path.exists(key):
                use_flag = False
                break
        if use_flag:
            return fitting_dict
        else:
            print("some of the pdb files in the dict are missing, re-searching")
    # reorganize data
    fasta_path = os.path.abspath(params["P"])
    chain_dict = read_fasta(fasta_path)
    if len(chain_dict) == 0:
        return {}
    fasta_path = os.path.join(single_chain_pdb_dir, "input.fasta")
    write_all_fasta(chain_dict, fasta_path)

    # first do experimental database search
    if not params["af_only"]:
        output_path = os.path.join(single_chain_pdb_dir, "exp_search.out")
        search_command = "blastp -query %s -db %s -out %s -num_threads %d" % (
            fasta_path,
            params["db_exp_path"],
            output_path,
            params["search_thread"],
        )
        if not os.path.exists(output_path) or os.path.getsize(output_path) < 1000:
            os.system(search_command)
        exp_match_dict = parse_blast_output(output_path)
    else:
        exp_match_dict = {}  # empty dict
    matched_dict = {}  # [key]: chain id list, [value]: the structure id
    fitting_dict = {}
    best_candidate_info = {}  # Store best candidate info for each chain
    # parse the information

    # merge only identical chains, otherwise, rely on combine search
    for key in exp_match_dict:
        current_match_list = exp_match_dict[key]
        print("match candidate:", key, current_match_list)
        for k in range(len(current_match_list)):
            match_id, evalue = current_match_list[k]
            if evalue == 0:
                # fetch the pdb to see if the protein size really reasonable
                expected_seq_length = (
                    len(chain_dict[key]) * params["search"]["length_ratio"]
                )
                max_allow_length = (
                    len(chain_dict[key]) * params["search"]["max_length_ratio"]
                )
                chain_name_list = key.replace(",", "-")
                current_chain_dir = os.path.join(
                    single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit]
                )
                mkdir(current_chain_dir)
                final_pdb_path = os.path.join(
                    single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb"
                )
                download_pdb(match_id, current_chain_dir, final_pdb_path)
                actual_structure_length = count_residues(final_pdb_path)
                if (
                    actual_structure_length >= expected_seq_length
                    and actual_structure_length <= max_allow_length
                ):
                    matched_dict[key] = "PDB:" + match_id
                    final_chain_list = chain_name_list.split("-")
                    fitting_dict[final_pdb_path] = final_chain_list
                    # Track best candidate info
                    best_candidate_info[key] = {
                        "chain_id": key,
                        "match_id": match_id,
                        "matched_id": matched_dict[key],
                        "actual_structure_length": actual_structure_length,
                        "expected_seq_length": expected_seq_length,
                        "max_allow_length": max_allow_length,
                        "num_res_difference": 0,  # Perfect match
                        "evalue": evalue,
                        "download_link": None,
                        "database": "PDB",
                        "is_good_match": True
                    }
                    break
                else:
                    os.remove(final_pdb_path)
                    functToDeleteItems(current_chain_dir)
            if evalue > 0:
                break
    matched_keys = matched_dict.keys()
    # write a new fasta to search

    remain_chain_dict = {k: v for k, v in chain_dict.items() if k not in matched_keys}
    print(
        "experimental db search finished, get %s/%s matched single chain structure"
        % (len(matched_keys), len(chain_dict))
    )
    if len(remain_chain_dict) > 0:
        print(
            "continue PDB+AFDB search for remained %d chains" % len(remain_chain_dict)
        )
        remain_fasta_path = os.path.join(single_chain_pdb_dir, "remain_search.fasta")
        write_all_fasta(remain_chain_dict, remain_fasta_path)
        output_path = os.path.join(single_chain_pdb_dir, "expaf_search.out")
        search_command = "blastp -query %s -db %s -out %s -num_threads %d" % (
            remain_fasta_path,
            params["db_path"],
            output_path,
            params["search_thread"],
        )
        if not os.path.exists(output_path) or os.path.getsize(output_path) < 1000:
            os.system(search_command)
        expaf_match_dict = parse_blast_output(output_path)
        if len(expaf_match_dict) != len(remain_chain_dict):
            for tmp_key in remain_chain_dict:
                if tmp_key not in expaf_match_dict:
                    print(
                        "warning, query %s can not find any templates in database"
                        % tmp_key
                    )

        for key in expaf_match_dict:
            closest_choice = None
            num_res_difference = 999999999
            current_match_list = expaf_match_dict[key]
            print("match candidate:", key, current_match_list)
            chain_name_list = key.replace(",", "-")
            current_chain_dir = os.path.join(
                single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit]
            )
            mkdir(current_chain_dir)
            final_pdb_path = os.path.join(
                single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb"
            )
            expected_seq_length = (
                len(chain_dict[key]) * params["search"]["length_ratio"]
            )
            max_allow_length = (
                len(chain_dict[key]) * params["search"]["max_length_ratio"]
            )
            
            # Prefetch metadata and download URLs for all candidates
            valid_candidates = []
            max_workers = min(len(current_match_list), params.get("search_thread", 16))
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                candidate_info = {}  # Map futures to (match_id, evalue, index)
                
                for k in range(len(current_match_list)):
                    match_id, evalue = current_match_list[k]
                    if params["af_only"]:
                        if "AFDB" not in match_id:
                            continue
                    
                    if "AFDB" in match_id:
                        split_info = match_id.split(":")
                        database = split_info[0]
                        pdb_id = split_info[1]
                        uniprot_id = pdb_id.split("-")[1]
                        
                        future = executor.submit(fetch_afdb_candidate, match_id, pdb_id, uniprot_id)
                        futures.append(future)
                        candidate_info[future] = (match_id, evalue, k)
                    else:
                        future = executor.submit(fetch_pdb_candidate, match_id)
                        futures.append(future)
                        candidate_info[future] = (match_id, evalue, k)
                
                # Collect results as they complete
                for future in as_completed(futures):
                    match_id, evalue, k = candidate_info[future]
                    try:
                        result = future.result()
                        if result is not None:
                            result["evalue"] = evalue
                            result["index"] = k
                            valid_candidates.append(result)
                        else:
                            if "AFDB" in match_id:
                                pdb_id = match_id.split(":")[1]
                                print("get metadata/download URL failed for %s, skipping" % pdb_id)
                            else:
                                print("get metadata failed for %s, skipping" % match_id)
                    except Exception as e:
                        print(f"Error processing candidate {match_id}: {e}")
            
            # Sort valid candidates by their original index to maintain order
            valid_candidates.sort(key=lambda x: x["index"])
            # Remove index key after sorting
            for candidate in valid_candidates:
                candidate.pop("index", None)
            
            if len(valid_candidates) == 0:
                print("No valid candidates with metadata/download URL for %s" % key)
                continue
            
            # Process valid candidates
            top1_structure_length = None
            for k, candidate in enumerate(valid_candidates):
                match_id = candidate["match_id"]
                actual_structure_length = candidate["actual_structure_length"]

                if k == 0:
                    # save the top 1 in case of no length match
                    top1_structure_length = actual_structure_length

                curr_seq_length_diff = abs(
                    expected_seq_length - actual_structure_length
                )

                if expected_seq_length <= actual_structure_length <= max_allow_length:
                    if "AFDB" in match_id:
                        matched_dict[key] = match_id
                    else:
                        matched_dict[key] = "PDB:" + match_id
                    final_chain_list = chain_name_list.split("-")
                    fitting_dict[final_pdb_path] = final_chain_list
                    closest_choice = candidate
                    # Track best candidate info
                    best_candidate_info[key] = {
                        "chain_id": key,
                        "match_id": match_id,
                        "matched_id": matched_dict[key],
                        "actual_structure_length": actual_structure_length,
                        "expected_seq_length": expected_seq_length,
                        "max_allow_length": max_allow_length,
                        "num_res_difference": 0,  # Perfect match
                        "evalue": candidate["evalue"],
                        "download_link": candidate.get("download_link"),
                        "database": "AFDB" if "AFDB" in match_id else "PDB",
                        "is_good_match": True
                    }
                    break  # already in match dict, no need to check anymore

                if curr_seq_length_diff < num_res_difference:
                    num_res_difference = curr_seq_length_diff
                    closest_choice = candidate
            
            if closest_choice is not None and key not in matched_dict:
                # Recalculate difference to ensure it's accurate
                num_res_difference = abs(
                    expected_seq_length - closest_choice["actual_structure_length"]
                )
                print(
                    "Closest choice is %s with %d difference"
                    % (closest_choice["match_id"], num_res_difference)
                )

            if key not in matched_dict:  # no match under length condition
                if params["af_only"] and closest_choice is None:
                    # we can only pick the top 1 candidate
                    if len(valid_candidates) > 0:
                        closest_choice = valid_candidates[0]
                        num_res_difference = abs(
                            expected_seq_length - closest_choice["actual_structure_length"]
                        )
                
                if closest_choice is not None:
                    match_id = closest_choice["match_id"]
                    # Ensure difference is recalculated
                    num_res_difference = abs(
                        expected_seq_length - closest_choice["actual_structure_length"]
                    )
                    if "AFDB" in match_id:
                        matched_dict[key] = match_id
                    else:
                        matched_dict[key] = "PDB:" + match_id
                    # Track best candidate info
                    best_candidate_info[key] = {
                        "chain_id": key,
                        "match_id": match_id,
                        "matched_id": matched_dict[key],
                        "actual_structure_length": closest_choice["actual_structure_length"],
                        "expected_seq_length": expected_seq_length,
                        "max_allow_length": max_allow_length,
                        "num_res_difference": num_res_difference,
                        "evalue": closest_choice["evalue"],
                        "download_link": closest_choice.get("download_link"),
                        "database": "AFDB" if "AFDB" in match_id else "PDB",
                        "is_good_match": False
                    }
                    print(
                        "we have no better choice but pick %s with %d residues differences"
                        % (match_id, num_res_difference)
                    )
                else:
                    print("No valid candidate found for %s" % key)
                    continue

            # download the closest choice as final
            if closest_choice is not None:
                match_id = closest_choice["match_id"]
                if "AFDB" in match_id:
                    download_link = closest_choice["download_link"]
                    download_flag = download_file(download_link, final_pdb_path)
                    if download_flag is False:
                        print("Download failed for %s, removing from matched_dict" % match_id)
                        if key in matched_dict:
                            del matched_dict[key]
                        if final_pdb_path in fitting_dict:
                            del fitting_dict[final_pdb_path]
                        if key in best_candidate_info:
                            del best_candidate_info[key]
                        time.sleep(60)
                        continue
                else:
                    download_pdb(match_id, current_chain_dir, final_pdb_path)

    print("DB search finished! Match relationship ", matched_dict)
    # get the matched dict

    for chain_name_list in chain_dict:
        if chain_name_list not in matched_dict:
            print("*" * 100)
            print(
                "WARNING! query %s can not find any templates in database"
                % chain_name_list
            )
            print("*" * 100)
            continue
        matched_id = matched_dict[chain_name_list]
        chain_name_list = chain_name_list.replace(",", "-")
        current_chain_dir = os.path.join(
            single_chain_pdb_dir, str(chain_name_list)[:max_file_name_limit]
        )
        mkdir(current_chain_dir)
        final_pdb_path = os.path.join(
            single_chain_pdb_dir, chain_name_list[:max_file_name_limit] + ".pdb"
        )
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
            download_link = get_afdb_pdb_url(pdb_id)
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
    
    # Save best candidate info to JSON
    best_candidate_json_path = os.path.join(save_path, "best_candidates.json")
    with open(best_candidate_json_path, "w") as f:
        json.dump(best_candidate_info, f, indent=2)
    print(f"Best candidate info saved to {best_candidate_json_path}")
    
    return fitting_dict
