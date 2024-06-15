#This code is from https://github.com/DSIMB/SWORD2/blob/main/SWORD2.py
#Some functions are removed for faster processing

#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import logging
import math
import multiprocessing
import os
import random
import re
import shlex
import shutil
import plotly.express as px
import subprocess
import sys
import textwrap
import time
import json
from copy import copy
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib import patches
from prody import *
from multiprocessing import cpu_count
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def check_parsing_pdb(uniprot_id, mgnify_id, pdb_id, pdb_chain, model, input_file):
    """
    This function tries to fetch and parse the input PDB submitted,
    either PDB code and chain or a whole PDB file downloaded by the user.

    Args:
        - uniprot_id (str): AlphaFold Uniprot Accession Id
        - mgnify_id (str): MGnifyID for the ESM Metagenomic Atlas
        - pdb_id (str): PDB code to fetch and parse
        - pdb_chain (str): PDB chain to fetch and parse
        - model (int): Structure model to parse
        - input_file (str): The path to the file downloaded by the user


    Returns:
        - prot (ProDy Protein object): if fetched and parsed correctly, else False
    """
    # Custom user file
    if input_file:
        logging.info("Try to parse user input structure file")
        file_ext = os.path.splitext(input_file)[1]
        try:
            if file_ext in [".cif", ".mmcif"]:
                prot = parseMMCIF(input_file, chain=pdb_chain, model=model)
            else:
                prot = parsePDB(input_file, chain=pdb_chain, model=model)
        except Exception as e:
            sys.exit(str(e))
        if prot is None:
            # Try again without specifying the chain. If it works and
            # if the protein has no chain, we set it to "A" by default.
            # It is necessary because the scoring program
            # accepts only PDB files which contain a chain
            if file_ext in [".cif", ".mmcif"]:
                prot = parseMMCIF(input_file, model=model)
            else:
                prot = parsePDB(input_file, model=model)
            if prot is not None and prot.getChids()[0] == ' ':
                prot.setChids([pdb_chain for i in range(prot.numAtoms())])
            else:
                sys.exit(f"Atomic data could not be parsed for chain {pdb_chain}. Please check the input PDB file")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
    # Download the AlphaFold model
    elif uniprot_id:
        logging.info(f"Download AlphaFold Uniprot Accession ID: {uniprot_id}")
        ok, response = download_af_model(uniprot_id)
        if not ok:
            sys.exit(f"Error: {response}. Please try again.")
        else:
            input_file = response
        # Redirect useful error of ProDy
        try:  # Try to parse PDB file
            prot = parsePDB(input_file, chain=pdb_chain, model=model)
        except Exception as e:
            sys.exit(str(e))
        # A parsed PDB by ProDy returns an AtomGroup, else it could be an EMD file, which we don't want...
        if type(prot) is not AtomGroup:
            sys.exit(f"Error: Something went wrong with the AlphaFold Uniprot Accession Id {uniprot_id}. Please make sure chain {pdb_chain} is valid and that it is a valid ID referenced in the AlphaFold database (https://alphafold.ebi.ac.uk/)")
        if prot is None:
            sys.exit(f"Error: Atomic data could not be parsed. Please check the PDB file corresponding to the code {pdb_id}. Also check that it actually contains the chain {pdb_chain}.")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
        input_file = os.path.basename(input_file)
    # Download the ESMFold model
    elif mgnify_id:
        logging.info(f"Download the ESM Metagenomic Atlas ID: {mgnify_id}")
        ok, response = download_esm_model(mgnify_id)
        if not ok:
            sys.exit(f"Error: {response}. Please try again.")
        else:
            input_file = response
        # Redirect useful error of ProDy
        try:  # Try to parse PDB file
            prot = parsePDB(input_file, chain=pdb_chain, model=model)
        except Exception as e:
            sys.exit(str(e))
        # A parsed PDB by ProDy returns an AtomGroup, else it could be an EMD file, which we don't want...
        if type(prot) is not AtomGroup:
            sys.exit(f"Error: Something went wrong with the ESM Metagenomic Atlas ID (MGnifyID) {mgnify_id}. Please make sure chain {pdb_chain} is valid and that it is a valid ID referenced in the ESM Metagenomic Atlas database (https://esmatlas.com/)")
        if prot is None:
            sys.exit(f"Error: Atomic data could not be parsed. Please check the PDB file corresponding to the code {pdb_id}. Also check that it actually contains the chain {pdb_chain}.")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
        input_file = os.path.basename(input_file)
    # Fetch and parse a PDB from a given PDB code and chain
    else:
        logging.info(f"Fetch PDB ID: {pdb_id}")
        # Redirect useful error of ProDy
        try:  # Try to parse fetched PDB file
            prot = parsePDB(pdb_id, chain=pdb_chain, model=model, folder=RESULTS_DIR, compressed=False)
        except Exception as e:
            sys.exit(str(e))
        # A parsed PDB by ProDy returns an AtomGroup, else it could be an EMD file, which we don't want...
        if type(prot) is not AtomGroup:
            sys.exit(f"Error: No PDB file could be parsed. Please check that the PDB code {pdb_id} exists in the PDB RCSB database (https://www.rcsb.org) with a legacy PDB format file available, and contains the chain {pdb_chain}. Careful, 'A' is different than 'a'. Please note that mmCIF files are not yet supported.")
        if prot is None:
            sys.exit(f"Error: Atomic data could not be parsed. Please check the PDB file corresponding to the code {pdb_id}. Also check that it actually contains the chain {pdb_chain}.")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
    return prot

def requests_retry_session(retries=3,
                           backoff_factor=0.3,
                           status_forcelist=(500, 502, 504),
                           session=None):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def download_af_model(id):
    """
    Download the Alphafold2 model corresponding to the Uniprot Id given by user
    https://alphafold.ebi.ac.uk/

    Returns:
        - File path (string): Path of the downloaded PDB file
        or
        False if wrong id
        "DOWNLOAD ERROR" if could not download
    """
    name = f"AF-{id}-F1-model_v4"
    url = f"https://alphafold.ebi.ac.uk/files/{name}.pdb"
    try:
        response = requests_retry_session().get(url)
    except Exception as x:
        return (False, x)
    with open(f"{RESULTS_DIR}/{name}.pdb", "w") as f:
        f.write(response.text)
    return (True, f"{RESULTS_DIR}/{name}.pdb")

def download_esm_model(id):
    """
    Download the ESM-2 model corresponding to the MGnify Id given by user
    https://esmatlas.com

    Returns:
        - File path (string): Path of the downloaded PDB file
        or
        False if wrong id
        "DOWNLOAD ERROR" if could not download
    """
    url = "https://api.esmatlas.com/fetchPredictedStructure/"+id
    try:
        response = requests_retry_session().get(url)
    except Exception as x:
        return (False, x)
    with open(f"{RESULTS_DIR}/{id}.pdb", "w") as f:
        f.write(response.text)
    return (True, f"{RESULTS_DIR}/{id}.pdb")

def parse_sword(output):
    """
    Parse the SWORD output and return a dictionary.

    Args:
        - output (list of str): The output of SWORD as a list of strings

    Returns:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
    """
    logging.info("Parse SWORD output")
    # Index of the alternative partitionings
    nb_alt = 0
    # Dict containing all the informations given by SWORD2 output
    sword_results = {"DOMAINS": {}, "AMBIGUITY": "n/a"}
    for line in output:
        amb_found = re.search(r"^A-index = (\++)$", line)
        ass_found = re.search(r"\d{1,}\s+\|", line)
        # Found Ambiguity-index line
        if amb_found:
            sword_results["AMBIGUITY"] = amb_found.group(1)
        # Found a domain assignment
        elif ass_found:
            splitted_ass = re.split(r"\s{0,}\|\s{0,}", line)
            sword_results["DOMAINS"][nb_alt] = {}
            sword_results["DOMAINS"][nb_alt]["NB_DOMAINS"] = int(splitted_ass[0].strip())
            sword_results["DOMAINS"][nb_alt]["MIN_SIZE"] = int(splitted_ass[1].strip())
            boundaries = re.split(r"\s", splitted_ass[2].strip())
            sword_results["DOMAINS"][nb_alt]["BOUNDARIES"] = {}
            for i, boundary in enumerate(boundaries):
                mult_boundaries = re.split(r";", boundary)
                sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i] = []
                for mb in mult_boundaries:
                    start_pu, end_pu = re.split(r"-", mb)
                    sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i].append((int(start_pu), int(end_pu)))
            sword_results["DOMAINS"][nb_alt]["AVERAGE K"] = float(splitted_ass[3].strip())
            sword_results["DOMAINS"][nb_alt]["QUALITY"] = splitted_ass[4].strip()
            nb_alt += 1
    return sword_results


def get_quality_as_nb_bars(quality):
    """
    Transform "*****" or "+++" into number.
    """
    return len(quality) if quality != "n/a" else 0


def write_partitionings_json(sword_results, energies):
    """
    Write the partitionnings into JSON formatted file

    Args:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
    Returns:
        - None
    """
    partitioning = os.path.join(RESULTS_DIR, "SWORD2_summary.json")
    with open(partitioning, "w") as f:
        json_results = {}
        # AMBIGUITY INDEX
        nb_bars = get_quality_as_nb_bars(sword_results["AMBIGUITY"])
        json_results["Ambiguity index"] = "*" * nb_bars
        for nb_alt_part, alt_part in sword_results["DOMAINS"].items():
            alt_part_json = {}
            # OPTIMAL PARTITIONING
            if nb_alt_part == 0:
                alt_part_json["Partition"] = "Optimal partition"
            # ALTERNATIVE PARTITIONINGS
            else:
                alt_part_json["Partition"] = f"Alternative partition {nb_alt_part}"
            nb_bars = get_quality_as_nb_bars(alt_part["QUALITY"])
            alt_part_json["Quality"] = "*" * nb_bars
            alt_part_json["Nb. domains"] = len(alt_part["BOUNDARIES"])
            alt_part_json['BOUNDARIES']=alt_part["BOUNDARIES"]
            # domains_json = {}
            # for i, dom in alt_part["BOUNDARIES"].items():
            #     domain_json = []
            #     domain_json.append(("AUL", int((1-(1/(energies[(nb_alt_part, i)][1])**2))*100) if abs(energies[(nb_alt_part, i)][1]) >= 1 else 0))
            #     domain_json.append(("Z-score", round(energies[(nb_alt_part, i)][1], 1)))
            #     p_unit_json = {}
            #     for j, (start_pu, end_pu) in enumerate(dom):
            #         p_unit_json[f"{str(start_pu)}-{str(end_pu)}"] = {"AUL": int((1-(1/(energies[(nb_alt_part, i, start_pu, end_pu)][1])**2))*100) if abs(energies[(nb_alt_part, i, start_pu, end_pu)][1]) >= 1 else 0, "Z-score": round(energies[(nb_alt_part, i, start_pu, end_pu)][1], 1)}
            #     domain_json.append(("PUs", p_unit_json))
            #     domains_json[f"Domain {i+1}"] = dict(domain_json)
            # alt_part_json["Domains"] = domains_json
            json_results[alt_part_json["Partition"]] = alt_part_json
        f.write(json.dumps(json_results, indent=4))


def define_colors(sword_results):
    """
        Set visually distinct colors for Domains (dark palette) and PUs (ligh palette)

        Return:
            - pus_colors (dict): key=(start_pu, end_pu) --> value=(r, g, b)
            - dom_colors (dict): key=domain_id --> value=(r, g, b)
    """
    pus_colors = {}
    dom_colors = {}
    color_domain_cnt = 0
    color_pu_cnt = 0
    dom_id = None
    # Set of 20 visually distinct colors for PUs (light palette)
    # https://medialab.github.io/iwanthue/
    colors_for_pus = ["#baeae5", "#e1c65b", "#b4bcf7", "#d0e47b", "#f0a8e5", "#6de4ac", "#d8c0e4", "#a5e18d", "#68d1f1", "#f3b175", "#63e3d8", "#ebbaba", "#c3d28c", "#aac5e2", "#e8da92", "#bcdbec", "#e1c298", "#98c7c6", "#abddb4", "#d4d8bb"]

    # Set of 25 visually distinct colors for domains (dark palette)
    # https://medialab.github.io/iwanthue/
    colors_for_domains = ["#27a3b4", "#c08423", "#d83e7c", "#986a35", "#686fdf", "#559a3b", "#763da6", "#8d8d36", "#ce61c7", "#406021", "#a42c88", "#3d956b","#341d79", "#cf3a44", "#3c8cc9", "#cf6430", "#4d4b92", "#7d3119", "#7b81cd", "#cf6c61", "#401d56", "#c95f7a", "#792e65", "#82263a", "#b86da8"]
    for i, part in sword_results["DOMAINS"].items():
        # COLOR ALL PUS OF A SWORD ALTERNATIVE DOMAIN WITH A DIFFERENT COLOR
        for j, dom in part["BOUNDARIES"].items():
            # Consider that a domain is a list of PUs delineation sorted by 1st delimitation of PUs
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            if dom_id not in dom_colors:
                # Pick a new color
                hex_val = colors_for_domains[color_domain_cnt].lstrip("#")
                (r, g, b) = tuple(int(hex_val[k:k+2], 16) for k in (0, 2, 4))
                dom_colors[dom_id] = (r, g, b)
                color_domain_cnt += 1
        for j, dom in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            r, g, b = dom_colors[dom_id]
        for j, dom in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            for start_pu, end_pu in dom:
                if (start_pu, end_pu) not in pus_colors:
                    # Pick a new color
                    hex_val = colors_for_pus[color_pu_cnt].lstrip("#")
                    (r, g, b) = tuple(int(hex_val[k:k+2], 16) for k in (0, 2, 4))
                    pus_colors[(start_pu, end_pu)] = (r, g, b)
                    color_pu_cnt += 1
    return pus_colors, dom_colors


def write_domains_histogram(sword_results, domains_colors):
    """
    Generate histogram of SWORD domains consistency.
    """
    logging.info("Generate histogram of SWORD2 domains consistency")
    histogram = os.path.join(RESULTS_DIR, "domains_histogram.png")
    domains = {}
    for i, part in sword_results["DOMAINS"].items():
        for j, domain in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(domain, key=lambda x: x[0]))
            if dom_id not in domains:
                domains[dom_id] = {}
            if "nb" not in domains[dom_id]:
                domains[dom_id]["nb"] = 1
            else:
                domains[dom_id]["nb"] += 1
    sorted_domains_by_decreasing_count = sorted(domains.items(), key=lambda x: x[1]["nb"], reverse=True)
    x = []
    for dom_part, infos in sorted_domains_by_decreasing_count:
        x.append(str(list(dom_part)).strip("[]"))
    y = []
    for dom_part, infos in sorted_domains_by_decreasing_count:
        y.append(infos["nb"])
    colors = {}
    for dom_part, infos in sorted_domains_by_decreasing_count:
        r, g, b = domains_colors[dom_part]
        colors[dom_part] = f"rgb({r}, {g}, {b})"
    
    # Calling DataFrame constructor after zipping
    # both lists, with columns specified
    df = pd.DataFrame(list(zip(x, y)),
                columns =['SWORD Domains', 'Count'])
    fig = px.bar(df, x="SWORD Domains", y="Count", color=colors, title="Consistency of domains determined by SWORD")
    fig.update_layout(showlegend=False)
    fig.write_image(histogram, scale=4)


def predict_time(prot):
    """ Time predicted in function of len of protein in minutes"""
    return str(int(133*math.exp(2.06*10**-3*len(set(prot.getResnums())))/60))


def get_energy_and_z_score(BIN_DIR, pdb, res_list=None):
    """
    Calculate pseudo-energy and z-score of a protein or specified residue list
    of the input pdb.
    List is res+chain: 10A,11A,12A,13A
    """
    if res_list:
        cmd_args = f"{BIN_DIR}/mypmfs-master/scoring -i {pdb} -d {BIN_DIR}/mypmfs-master/025_30_100_potential -q {res_list} -z -s 2000"
    else:
        cmd_args = f"{BIN_DIR}/mypmfs-master/scoring -i {pdb} -d {BIN_DIR}/mypmfs-master/025_30_100_potential -z -s 2000"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8")
    output = output.split("\n")
    energy = None
    z_score = None
    for line in output:
        pseudo_e_found = re.search(r"^Pseudo-energy = (.+)$", line)
        z_score_found = re.search(r"^Z-score = (.+)$", line)
        if pseudo_e_found:
            energy = float(pseudo_e_found.group(1))
        if z_score_found:
            z_score = float(z_score_found.group(1))
    return energy, z_score


def multiprocess_get_energy(i, pdb_chain, pdb_id_chain, RESULTS_DIR, BIN_DIR, energies, dom_bounds):
    """
    Calculate the energy and Z-score of PUs and Domains.

    Args:
        - i: index of partitionning
        - energies: dictionary to hold the results
        - dom_bounds: Boundaries of domains and PUs to calculate
    """
    j, domain = dom_bounds
    # Residues of the domain gathered little by little
    dom_residues = ""
    # Python multiprocessing is shit and cannot handle dict of dict.
    # So I have to use a tmp list for this little buddy...
    tmp_list = []
    for start_pu, end_pu in domain:
        pu_residues = ""
        dom_residues += ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)]) + ","
        pu_residues += ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)])
        pu_energy, pu_z_score = get_energy_and_z_score(BIN_DIR, f"{RESULTS_DIR}/{pdb_id_chain}", pu_residues)
        energies[(i, j, start_pu, end_pu)] = []
        tmp_list = energies[(i, j, start_pu, end_pu)]
        tmp_list.extend([pu_energy, pu_z_score])
        energies[(i, j, start_pu, end_pu)] = tmp_list
    dom_energy, dom_z_score = get_energy_and_z_score(BIN_DIR, f"{RESULTS_DIR}/{pdb_id_chain}", dom_residues)
    energies[(i, j)] = []
    tmp_list = energies[(i, j)]
    tmp_list.extend([dom_energy, dom_z_score])
    energies[(i, j)] = tmp_list


def write_peeling_results():
    """
    Parse Protein Peeling 3 results and calculate pseudo energy and AUL for
    all Protein Units.

    Args:
        - energies (dict): 
    """
    # Get old resnums
    peeling_num = os.path.join(RESULTS_DIR, "PDBs_Clean", pdb_id_chain, f"{pdb_id_chain}.num")
    with open(peeling_num, "r") as f:
        ori_resnums = [int(resnum) for resnum in f.readline().split()]
    # Parse Peeling.log
    peeling_results = {}
    peeling_log = os.path.join(RESULTS_DIR, "PDBs_Clean", pdb_id_chain, "Peeling", "Peeling.log")
    with open(peeling_log, "r") as f:
        # Index of the alternative partitionings
        nb_lvl = 1
        # Dict containing all the informations given by Peeling output
        for line in f:
            if not line.startswith("#") and line != "\n":
                line = line.split()
                peeling_results[nb_lvl] = {}
                peeling_results[nb_lvl]["i/e"] = float(line[0])
                peeling_results[nb_lvl]["i/i+e"] = float(line[1])
                peeling_results[nb_lvl]["R2"] = float(line[2])
                peeling_results[nb_lvl]["CI"] = float(line[3])
                peeling_results[nb_lvl]["N"] = int(line[4])
                # Retrieve only the PUs
                # store them as a list of tuples
                peeling_results[nb_lvl]["PUs"] = [(ori_resnums[int(line[5+i]) - 1], ori_resnums[int(line[5+i+1]) - 1]) for i in range(0, len(line[5:])-1, 2)]
                # Sorte by inceasing boundaries
                peeling_results[nb_lvl]["PUs"] = sorted(peeling_results[nb_lvl]["PUs"], key=lambda x: x[0]) 
                nb_lvl += 1
    
    # WRITE THE RESULTS
    logging.info("Calculate pseudo-energies of PUs")
    peeling = os.path.join(RESULTS_DIR, "PEELING_summary.txt")
    with open(peeling, "w") as f:
        for lvl, data in peeling_results.items():
            f.write(f"""Peeling level {lvl}\n    Number of Protein Units: {data["N"]}\n    Compaction Index: {round(data["CI"], 2)}\n""")
            for start_pu, end_pu in data["PUs"]:
                pu_residues = ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)])
                pu_energy, pu_z_score = get_energy_and_z_score(BIN_DIR, f"{RESULTS_DIR}/{pdb_id_chain}", pu_residues)
                f.write(f"""    {str(start_pu)+"-"+str(end_pu):>7}: AUL={int((1-(1/(pu_z_score)**2))*100) if abs(pu_z_score) >= 1 else 0:3}% Z-score={round(pu_z_score, 1)}\n""")
    logging.info("Write the Peeling results")

if __name__ == '__main__':

    #################
    # Parse arguments
    #################

    def check_cpu(nb_cpu):
        """
        Check if the user input CPU nb is valid
        """
        try:
            nb_cpu = int(nb_cpu)
        except ValueError as e:
            print("Unable to cast ", type(nb_cpu), " into integer (int)")
            raise argparse.ArgumentTypeError(
                "Error option -c/--cpu: please input an integer or string integer"
            ) from e
        if 0 <= nb_cpu <= cpu_count():
            return nb_cpu
        raise argparse.ArgumentTypeError(
            f"Error option -c/--cpu: nb_cpu should be 0 <= nb_cpu <= {cpu_count()}")


    def check_model(model):
        """
        Check if the user input model nb is valid
        """
        try:
            model = int(model)
        except ValueError as e:
            print("Unable to cast ", type(model), " into integer (int)")
            raise argparse.ArgumentTypeError(
                "Error option -c/--cpu: please input an integer or string integer"
            ) from e
        if model >= 1:
            return model
        raise argparse.ArgumentTypeError(
                f"Error option -d/--model: model should be >= 1")

    parser = argparse.ArgumentParser(
                description=textwrap.dedent('''\
                    SWORD2: SWift and Optimized Recognition of protein Domains.
                    The SWORD2 partitioning algorithm produces multiple alternative 
                    domain assignments for a given protein structure. 
                    This unique approach handles ambiguous protein structure partitioning, 
                    admitting several solutions. The decomposition of the protein structure
                    into domains is achieved through the hierarchical clustering of Protein Units, 
                    evolutionarily preserved structural descriptors at the interface between 
                    secondary structures and domains.'''),
                formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    optional = parser.add_argument_group('optional arguments')
    required = parser.add_argument_group('required arguments')
    parser.add_argument('--version', action='version', version='SWORD2 1.1.2')
    group.add_argument("-u", "--uniprot-id",
                        help=textwrap.dedent('''\
                                                AlphaFold Uniprot Accession Id.
                                                The corresponding predicted structure will be downloaded from the AlphaFold database.'''), type=str)
    group.add_argument("-m", "--mgnify-id",
                        help=textwrap.dedent('''\
                                                MGnify Id.
                                                The corresponding predicted structure will be downloaded from the ESM Metagenomic Atlas database.'''), type=str)
    group.add_argument("-p", "--pdb-id",
                        help=textwrap.dedent('''\
                                                PDB id.
                                                The corresponding structure will be downloaded from the PDB database.'''), type=str)
    group.add_argument("-i", "--input-file", help="Path to an input PDB or mmCIF file.", type=str)
    optional.add_argument("-c", "--pdb-chain", help="PDB chain. Default is A.", type=str, required=False, default="A")
    optional.add_argument("-d", "--model", help="Model to parse. Especially usefull for NMR files which contain several models. Default is 1.", type=check_cpu, required=False, default=1)
    optional.add_argument("-x", "--cpu", help=textwrap.dedent(f'''\
                                                                How many CPUs to use.
                                                                Default all (0). 
                                                                Max on this computer is: {cpu_count()}'''),
                        default=0, type=check_cpu, required=False)
    required.add_argument("-o", "--output", help=textwrap.dedent('''\
                                                                    Output directory.
                                                                    Results will be generated inside in a dedicated directory 
                                                                    named after OUTPUT/PDBID_CHAIN/'''), 
                        type=str, required=True)

    args = parser.parse_args()

    uniprot_id = args.uniprot_id
    mgnify_id = args.mgnify_id
    pdb_id = args.pdb_id
    input_file = args.input_file
    pdb_chain = args.pdb_chain
    model = args.model
    output_dir = args.output
    nb_cpu = args.cpu
    if nb_cpu == 0:
        nb_cpu = cpu_count()

    if input_file:
        if os.path.exists(input_file):
            pdb_id_chain = os.path.basename(os.path.splitext(input_file)[0]) + "_" + pdb_chain
        else:
            sys.exit("Unable to open file: " + input_file)
    elif uniprot_id:
        pdb_id_chain = uniprot_id + "_" + pdb_chain
    elif mgnify_id:
        pdb_id_chain = mgnify_id + "_" + pdb_chain
    else:
        pdb_id_chain = pdb_id + "_" + pdb_chain

    confProDy(verbosity="none")

    # Configure the logger to redirect to job log file
    logging.basicConfig(level=logging.INFO,
                        handlers=[
                            logging.StreamHandler()
                        ],
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt="%Y/%m/%d %H:%M:%S")

    # Be conservative, if results directory already exists, create another one with suffix
    RESULTS_DIR = os.path.join(output_dir, pdb_id_chain)
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
    else:
        name_rep = time.strftime("_%d_%m_%Y_") + ''.join(random.choice(list(map(str, list(range(0, 10))))) for i in range(5))
        RESULTS_DIR += name_rep
        os.makedirs(RESULTS_DIR)
        # Configure file log here
        fh = logging.FileHandler(os.path.join(RESULTS_DIR, "sword2.log"))
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(filename)s %(funcName)s() %(lineno)s - %(message)s")
        fh.setFormatter(formatter)
        logging.getLogger('').addHandler(fh)
        logging.warning(f"Results dir '{os.path.join(output_dir, pdb_id_chain)}' already exists. We created '{RESULTS_DIR}' instead.")


    ######################################
    # Define paths, variables, and logging
    ######################################

    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    BIN_DIR = os.path.join(BASE_DIR, "bin")
    SWORD_DIR = os.path.join(BIN_DIR, "SWORD/bin/SWORD")
    SWORD = os.path.join(SWORD_DIR, "SWORD")
    DISPLAY_SWORD2 = os.path.join(BIN_DIR, "display_SWORD2_output.pl")

    # CHECK ENTRIES
    prot = check_parsing_pdb(uniprot_id, mgnify_id, pdb_id, pdb_chain, model, input_file)

    # ESTIMATE THE RUNTIME
    est_time_in_minutes = int(predict_time(prot))
    if est_time_in_minutes < 1:
        est_time_in_minutes = 1
    logging.info(f"")
    logging.info(f">>>   Estimated runtime: {est_time_in_minutes} minutes")
    logging.info(f">>>   Using {nb_cpu} cpus")
    logging.info(f"")

    ######################################################
    # Write the specific chain as a new PDB file for SWORD
    ######################################################

    pdb_chain_file = os.path.join(RESULTS_DIR, pdb_id_chain + ".pdb")
    logging.info("Write a clean version of the PDB: remove non standard residues")
    # Remove residues which have insertion codes
    res_to_remove = " "
    hv = prot.getHierView()
    for residue in hv.iterResidues():
        if residue.getIcode():
            res_to_remove += f"{residue.getResnum()} "
    # Keep only what constitutes the protein, remove the non standard residues
    if res_to_remove != " ":
        prot = prot.select('protein and not nonstdaa and not hetatm and not resnum ' + res_to_remove)
    else:
        prot = prot.select('protein and not nonstdaa and not hetatm')
    # In case there are negative residue numbers, we renumber
    # for convenience
    hv = prot.getHierView()
    SEQ_RESNUMS = [residue.getResnum() for residue in hv.iterResidues()]
    NEW_SEQ_RESNUMS = [i for i, _ in enumerate(SEQ_RESNUMS, start=1)]
    # Renumber the residues from 0 -> len(prot)
    [residue.setResnum(NEW_SEQ_RESNUMS[i]) for i, residue in enumerate(hv.iterResidues())]
    # Write the cleaned protein
    file_written = writePDB(pdb_chain_file, prot)
    if file_written is None:
        logging.warning(f"PDB file {pdb_chain_file} could not be written using function writePDB of ProDy")
        sys.exit(1)
    # Remove the ".pdb" extension for SWORD
    os.rename(file_written, os.path.splitext(file_written)[0])

    # Now that the protein is clean, get some infos
    hv = prot.getHierView()
    ch = hv.getChain(pdb_chain)
    # Get the aa sequence
    seq = ch.getSequence()
    PROT_LEN = len(seq)
    # Map the original resnums of the PDB to new numerotation from 0 -> N
    # this is used for the selection/coloration of domains and PUs in the sequence
    # section of the 3D viewer
    SEQ_RESNUMS_NORM = {num: i for i, num in enumerate(SEQ_RESNUMS)}
    SEQ_RESNUMS_NORM_INV = {i: num for i, num in enumerate(SEQ_RESNUMS)}

    ##############
    # Launch SWORD
    ##############

    logging.info("Launch SWORD")
    # First run compiles necessary dependency for current arch
    if not os.path.exists(SWORD_DIR+"/bin/Dssp/dsspcmbi"):
        subprocess.run(SWORD, capture_output=True)
    cmd_args = f"{DISPLAY_SWORD2} '{SWORD} -i {pdb_id_chain} --dir {RESULTS_DIR} -max 9'"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8")
    output = output.split("\n")
    with open(f"{RESULTS_DIR}/sword.txt", "w") as f:
        for line in output:
            f.write(line + "\n")

    ####################
    # Parse SWORD output
    ####################

    sword_results = parse_sword(output)

    #####################################################
    # Calculate the energy and Z-score of PUs and Domains
    #####################################################

    logging.info("Calculate pseudo-energies of Domains")
    manager = multiprocessing.Manager()
    energies = manager.dict()
    for i, part in sword_results["DOMAINS"].items():
        with multiprocessing.Pool(processes=nb_cpu) as p:
            FUNC = partial(multiprocess_get_energy, i, pdb_chain, pdb_id_chain, RESULTS_DIR, BIN_DIR, energies)
            p.imap_unordered(FUNC, list(part["BOUNDARIES"].items()))
            p.close()
            p.join()

    ############################
    # Write  SWORD partitionings
    ############################

   
    write_partitionings_json(sword_results, energies)

    