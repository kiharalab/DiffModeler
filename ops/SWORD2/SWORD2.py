#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import json
import math
import multiprocessing
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
import time
from copy import copy
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import requests
from matplotlib import patches
from prody import *
import logging as log
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
log.disable(log.CRITICAL)
confProDy(verbosity="info")
log.disable(log.NOTSET)
from prody import LOGGER as logging
# Create a new formatter with the desired format
formatter = log.Formatter("%(asctime)s %(levelname)-8s %(message)s", datefmt="%Y/%m/%d %H:%M:%S")
# Get the current handler and set the new formatter
for handler in logging.getHandlers():
    handler.setFormatter(formatter)

def check_parsing_pdb(uniprot_id, mgnify_id, pdb_id, pdb_chain, model, input_file, output_dir):
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
        - output_dir (str): The output directory to save fetched files

    Returns:
        - prot (ProDy Protein object): if fetched and parsed correctly, else exits
        - pdb_chain (str): the chain that was used
    """
    prot = None
    # Custom user file
    if input_file:
        logging.info("Try to parse user input structure file")
        prot, pdb_chain = parse_structure(input_file, pdb_chain, model)
        prot = prot.select("protein and not nonstdaa")
        if prot is None:
            sys.exit(
                "Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
            )
    # Download the AlphaFold model
    elif uniprot_id:
        logging.info(f"Download AlphaFold Uniprot Accession ID: {uniprot_id}")
        ok, response = download_af_model(uniprot_id, output_dir)
        if not ok:
            sys.exit(f"Error: {response}. Please try again.")
        else:
            input_file = response
        prot, pdb_chain = parse_structure(input_file, pdb_chain, model)
        prot = prot.select("protein and not nonstdaa")
        if prot is None:
            sys.exit(
                "Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
            )
        input_file = os.path.basename(input_file)
    # Download the ESMFold model
    elif mgnify_id:
        logging.info(f"Download the ESM Metagenomic Atlas ID: {mgnify_id}")
        ok, response = download_esm_model(mgnify_id, output_dir)
        if not ok:
            sys.exit(f"Error: {response}. Please try again.")
        else:
            input_file = response
        prot, pdb_chain = parse_structure(input_file, pdb_chain, model)
        prot = prot.select("protein and not nonstdaa")
        if prot is None:
            sys.exit(
                "Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
            )
        input_file = os.path.basename(input_file)
    # Fetch and parse a PDB from a given PDB code and chain
    else:
        logging.info(f"Fetch PDB ID: {pdb_id}")
        try:
            # First fetch the PDB file
            pdb_file = fetchPDB(pdb_id, folder=output_dir, compressed=False)
            if not pdb_file:
                sys.exit(f"Error: Unable to fetch PDB ID {pdb_id}.")
            prot, pdb_chain = parse_structure(pdb_file, pdb_chain, model)
            prot = prot.select("protein and not nonstdaa")
            if prot is None:
                sys.exit(
                    "Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
                )
        except Exception as e:
            sys.exit(str(e))
    return prot, pdb_chain


def parse_structure(input_file, pdb_chain, model):
    """
    Parses a PDB or mmCIF file and selects the specified chain.
    If no chain is specified, the first chain in the file is used.

    Args:
        - input_file (str): Path to the PDB or mmCIF file
        - pdb_chain (str): Chain identifier (can be None)
        - model (int): Model number to parse

    Returns:
        - prot (ProDy AtomGroup): Parsed protein structure
        - pdb_chain (str): The chain identifier used
    """
    file_ext = os.path.splitext(input_file)[1]
    try:
        if file_ext in [".cif", ".mmcif"]:
            prot = parseMMCIF(input_file, model=model)
        else:
            prot = parsePDB(input_file, model=model)
    except Exception as e:
        sys.exit(str(e))
    if prot is None:
        sys.exit("Atomic data could not be parsed. Please check the input file.")
    chain_ids = np.unique(prot.getChids())
    if len(chain_ids) == 0:
        sys.exit("No chains found in the PDB file.")
    if pdb_chain is None:
        pdb_chain = chain_ids[0]
        logging.info(f"No chain specified. Using first chain '{pdb_chain}' in the PDB file.")
    elif pdb_chain not in chain_ids:
        logging.info(f"Chain {pdb_chain} not found in PDB file. Available chains: {', '.join(chain_ids)}")
        sys.exit(f"Chain {pdb_chain} not found in PDB file.")
    prot = prot.select('chain ' + pdb_chain)
    if prot is None:
        sys.exit(f"Error selecting chain {pdb_chain}. Please check the input PDB file.")
    return prot, pdb_chain


def requests_retry_session(
    retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504), session=None
):
    """Creates a requests session with retry logic."""
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def download_af_model(uniprot_id, output_dir):
    """
    Download the Alphafold2 model corresponding to the Uniprot Id given by user
    https://alphafold.ebi.ac.uk/

    Returns:
        - (bool, str): (True, file path) if successful, else (False, error message)
    """
    name = f"AF-{uniprot_id}-F1-model_v4"
    url = f"https://alphafold.ebi.ac.uk/files/{name}.pdb"
    try:
        response = requests_retry_session().get(url)
        response.raise_for_status()
    except Exception as x:
        return (False, str(x))
    file_path = f"{output_dir}/{name}.pdb"
    with open(file_path, "w") as f:
        f.write(response.text)
    return (True, file_path)


def download_esm_model(mgnify_id, output_dir):
    """
    Download the ESM-2 model corresponding to the MGnify Id given by user
    https://esmatlas.com

    Returns:
        - (bool, str): (True, file path) if successful, else (False, error message)
    """
    url = f"https://api.esmatlas.com/fetchPredictedStructure/{mgnify_id}"
    try:
        response = requests_retry_session().get(url)
        response.raise_for_status()
    except Exception as x:
        return (False, str(x))
    file_path = f"{output_dir}/{mgnify_id}.pdb"
    with open(file_path, "w") as f:
        f.write(response.text)
    return (True, file_path)


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
            sword_results["DOMAINS"][nb_alt]["NB_DOMAINS"] = int(
                splitted_ass[0].strip()
            )
            sword_results["DOMAINS"][nb_alt]["MIN_SIZE"] = int(splitted_ass[1].strip())
            boundaries = re.split(r"\s", splitted_ass[2].strip())
            sword_results["DOMAINS"][nb_alt]["BOUNDARIES"] = {}
            for i, boundary in enumerate(boundaries):
                mult_boundaries = re.split(r";", boundary)
                sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i] = []
                for mb in mult_boundaries:
                    start_pu, end_pu = re.split(r"-", mb)
                    sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i].append(
                        (int(start_pu), int(end_pu))
                    )
            sword_results["DOMAINS"][nb_alt]["AVERAGE K"] = float(
                splitted_ass[3].strip()
            )
            sword_results["DOMAINS"][nb_alt]["QUALITY"] = splitted_ass[4].strip()
            nb_alt += 1
    return sword_results


def get_quality_as_nb_bars(quality):
    """
    Transform "*****" or "+++" into number.
    """
    return len(quality) if quality != "n/a" else 0


def write_partitionings(sword_results, energies, disable_energies):
    """
    Write the partitionings into text file.

    Args:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
        - energies (dict): Dictionary of energies
        - disable_energies (bool): Whether to include energies in the output
    Returns:
        - None
    """
    logging.info("Write the SWORD results")
    partitioning = os.path.join(RESULTS_DIR, "SWORD2_summary.txt")
    with open(partitioning, "w") as f:
        # AMBIGUITY INDEX
        nb_bars = get_quality_as_nb_bars(sword_results["AMBIGUITY"])
        f.write("Ambiguity index: " + "*" * nb_bars + "\n")
        for nb_alt_part, alt_part in sword_results["DOMAINS"].items():
            f.write("-----------------------\n")
            # OPTIMAL PARTITIONING
            if nb_alt_part == 0:
                f.write("Optimal partition\n")
            # ALTERNATIVE PARTITIONINGS
            else:
                f.write(f"Alternative partition {nb_alt_part}\n")
            nb_bars = get_quality_as_nb_bars(alt_part["QUALITY"])
            f.write("Quality: " + "*" * nb_bars + "\n")
            f.write(f"Nb. domains: {len(alt_part['BOUNDARIES'])}\n")
            for i, dom in alt_part["BOUNDARIES"].items():
                if not disable_energies:
                    dom_energy = energies.get((nb_alt_part, i), [None, None])
                    dom_aul = (
                        int((1 - (1 / (dom_energy[1]) ** 2)) * 100)
                        if dom_energy[1] and abs(dom_energy[1]) >= 1
                        else 0
                    )
                    dom_z_score = round(dom_energy[1], 1) if dom_energy[1] else "n/a"
                    f.write(
                        f"Domain:{i+1}       AUL={dom_aul:3}% Z-score={dom_z_score}\n"
                    )
                else:
                    f.write(f"Domain:{i+1}\n")
                for start_pu, end_pu in dom:
                    if not disable_energies:
                        pu_energy = energies.get(
                            (nb_alt_part, i, start_pu, end_pu), [None, None]
                        )
                        pu_aul = (
                            int((1 - (1 / (pu_energy[1]) ** 2)) * 100)
                            if pu_energy[1] and abs(pu_energy[1]) >= 1
                            else 0
                        )
                        pu_z_score = round(pu_energy[1], 1) if pu_energy[1] else "n/a"
                        f.write(
                            f"    PU:{str(start_pu)+'-'+str(end_pu):>7} AUL={pu_aul:3}% Z-score={pu_z_score}\n"
                        )
                    else:
                        f.write(f"    PU:{str(start_pu)+'-'+str(end_pu):>7}\n")


def write_partitionings_json(sword_results, energies, disable_energies):
    """
    Write the partitionings into JSON formatted file.

    Args:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
        - energies (dict): Dictionary of energies
        - disable_energies (bool): Whether to include energies in the output
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
            domains_json = {}
            for i, dom in alt_part["BOUNDARIES"].items():
                domain_json = {}
                if not disable_energies:
                    dom_energy = energies.get((nb_alt_part, i), [None, None])
                    dom_aul = (
                        int((1 - (1 / (dom_energy[1]) ** 2)) * 100)
                        if dom_energy[1] and abs(dom_energy[1]) >= 1
                        else 0
                    )
                    dom_z_score = round(dom_energy[1], 1) if dom_energy[1] else "n/a"
                    domain_json["AUL"] = dom_aul
                    domain_json["Z-score"] = dom_z_score
                p_unit_json = {}
                for start_pu, end_pu in dom:
                    pu_key = f"{start_pu}-{end_pu}"
                    if not disable_energies:
                        pu_energy = energies.get(
                            (nb_alt_part, i, start_pu, end_pu), [None, None]
                        )
                        pu_aul = (
                            int((1 - (1 / (pu_energy[1]) ** 2)) * 100)
                            if pu_energy[1] and abs(pu_energy[1]) >= 1
                            else 0
                        )
                        pu_z_score = round(pu_energy[1], 1) if pu_energy[1] else "n/a"
                        p_unit_json[pu_key] = {"AUL": pu_aul, "Z-score": pu_z_score}
                    else:
                        p_unit_json[pu_key] = {}
                domain_json["PUs"] = p_unit_json
                domains_json[f"Domain {i+1}"] = domain_json
            alt_part_json["Domains"] = domains_json
            json_results[alt_part_json["Partition"]] = alt_part_json
        f.write(json.dumps(json_results, indent=4))


def define_colors(sword_results):
    """
    Set visually distinct colors for Domains and PUs.

    Return:
        - pus_colors (dict): key=(start_pu, end_pu) --> value=(r, g, b)
        - dom_colors (dict): key=domain_id --> value=(r, g, b)
    """
    pus_colors = {}
    dom_colors = {}
    color_domain_cnt = 0
    color_pu_cnt = 0
    colors_for_pus = [
        "#baeae5",
        "#e1c65b",
        "#b4bcf7",
        "#d0e47b",
        "#f0a8e5",
        "#6de4ac",
        "#d8c0e4",
        "#a5e18d",
        "#68d1f1",
        "#f3b175",
        "#63e3d8",
        "#ebbaba",
        "#c3d28c",
        "#aac5e2",
        "#e8da92",
        "#bcdbec",
        "#e1c298",
        "#98c7c6",
        "#abddb4",
        "#d4d8bb",
    ]
    colors_for_domains = [
        "#27a3b4",
        "#c08423",
        "#d83e7c",
        "#986a35",
        "#686fdf",
        "#559a3b",
        "#763da6",
        "#8d8d36",
        "#ce61c7",
        "#406021",
        "#a42c88",
        "#3d956b",
        "#341d79",
        "#cf3a44",
        "#3c8cc9",
        "#cf6430",
        "#4d4b92",
        "#7d3119",
        "#7b81cd",
        "#cf6c61",
        "#401d56",
        "#c95f7a",
        "#792e65",
        "#82263a",
        "#b86da8",
    ]
    for i, part in sword_results["DOMAINS"].items():
        # COLOR ALL PUS OF A SWORD ALTERNATIVE DOMAIN WITH A DIFFERENT COLOR
        for j, dom in part["BOUNDARIES"].items():
            # Consider that a domain is a list of PUs delineation sorted by 1st delimitation of PUs
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            if dom_id not in dom_colors:
                # Pick a new color
                hex_val = colors_for_domains[color_domain_cnt].lstrip("#")
                (r, g, b) = tuple(int(hex_val[k : k + 2], 16) for k in (0, 2, 4))
                dom_colors[dom_id] = (r, g, b)
                color_domain_cnt += 1
        for j, dom in part["BOUNDARIES"].items():
            for start_pu, end_pu in dom:
                if (start_pu, end_pu) not in pus_colors:
                    # Pick a new color
                    hex_val = colors_for_pus[color_pu_cnt].lstrip("#")
                    (r, g, b) = tuple(int(hex_val[k : k + 2], 16) for k in (0, 2, 4))
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
    sorted_domains = sorted(domains.items(), key=lambda x: x[1]["nb"], reverse=True)
    x = [str(list(dom_part)).strip("[]") for dom_part, _ in sorted_domains]
    y = [infos["nb"] for _, infos in sorted_domains]
    colors = [
        f"rgb({r}, {g}, {b})"
        for dom_part, _ in sorted_domains
        for r, g, b in [domains_colors[dom_part]]
    ]
    df = pd.DataFrame(list(zip(x, y)), columns=["SWORD Domains", "Count"])
    fig = px.bar(
        df,
        x="SWORD Domains",
        y="Count",
        color=colors,
        title="Consistency of domains determined by SWORD",
    )
    fig.update_layout(showlegend=False)
    fig.write_image(histogram, scale=4)


def predict_time_full(prot):
    """Predict time in seconds based on protein length when user runs SWORD2 completely."""
    x = len(set(prot.getResnums()))
    return int(16.8 + 0.163 * x - 4.3e-5 * x**2)

def predict_time_no_energy_no_plots(prot):
    """Predict time in seconds based on protein length when user runs SWORD2 with options -e and -l
    meaning without calculation of pseudo-energies and plots"""
    x = len(set(prot.getResnums()))
    return int(4.23 + 0.0412 * x - 8.42e-6 * x**2)

def get_energy_and_z_score(bin_dir, pdb, res_list=None):
    """
    Calculate pseudo-energy and z-score of a protein or specified residue list.

    Args:
        - bin_dir (str): Path to the binary directory
        - pdb (str): Path to the PDB file
        - res_list (str): Comma-separated list of residues

    Returns:
        - (float, float): Energy and Z-score
    """
    if res_list:
        cmd_args = f"{bin_dir}/mypmfs-master/scoring_omp -i {pdb} -d {bin_dir}/mypmfs-master/025_30_100_potential -q {res_list} -z -s 2000"
    else:
        cmd_args = f"{bin_dir}/mypmfs-master/scoring_omp -i {pdb} -d {bin_dir}/mypmfs-master/025_30_100_potential -z -s 2000"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8").split("\n")
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


def multiprocess_get_energy(
    i, pdb_chain, pdb_id_chain, results_dir, bin_dir, energies, dom_bounds
):
    """
    Calculate the energy and Z-score of PUs and Domains.

    Args:
        - i (int): Index of partitioning
        - pdb_chain (str): PDB chain identifier
        - pdb_id_chain (str): PDB ID and chain
        - results_dir (str): Results directory
        - bin_dir (str): Binary directory
        - energies (dict): Shared dictionary to store energies
        - dom_bounds (tuple): Boundaries of domains and PUs to calculate
    """
    j, domain = dom_bounds
    dom_residues = ""
    for start_pu, end_pu in domain:
        pu_residues = ",".join(
            [f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu + 1)]
        )
        dom_residues += pu_residues + ","
        pu_energy, pu_z_score = get_energy_and_z_score(
            bin_dir, f"{results_dir}/{pdb_id_chain}", pu_residues
        )
        energies[(i, j, start_pu, end_pu)] = [pu_energy, pu_z_score]
    dom_energy, dom_z_score = get_energy_and_z_score(
        bin_dir, f"{results_dir}/{pdb_id_chain}", dom_residues.rstrip(",")
    )
    energies[(i, j)] = [dom_energy, dom_z_score]


def write_peeling_results(disable_energies):
    """
    Parse Protein Peeling 3 results and calculate pseudo energy and AUL for
    all Protein Units.

    Args:
        - disable_energies (bool): Whether to calculate energies or not
    """
    peeling_num = os.path.join(
        RESULTS_DIR, "PDBs_Clean", pdb_id_chain, f"{pdb_id_chain}.num"
    )
    ori_resnums = []
    if os.path.exists(peeling_num):
        with open(peeling_num, "r") as f:
            ori_resnums = [int(resnum) for resnum in f.readline().split()]
    peeling_results = {}
    peeling_log = os.path.join(
        RESULTS_DIR, "PDBs_Clean", pdb_id_chain, "Peeling", "Peeling.log"
    )
    with open(peeling_log, "r") as f:
        next(f)
        nb_lvl = 1
        for line in f:
            if not line.startswith("#") and line.strip():
                line = line.split()
                peeling_results[nb_lvl] = {
                    "i/e": float(line[0]),
                    "i/i+e": float(line[1]),
                    "R2": float(line[2]),
                    "CI": float(line[3]),
                    "N": int(line[4]),
                    "PUs": sorted(
                        [
                            (
                                ori_resnums[int(line[5 + i]) - 1],
                                ori_resnums[int(line[5 + i + 1]) - 1],
                            )
                            for i in range(0, len(line[5:]) - 1, 2)
                        ],
                        key=lambda x: x[0],
                    ),
                }
                nb_lvl += 1

    logging.info("Write Peeling results")
    peeling = os.path.join(RESULTS_DIR, "PEELING_summary.txt")
    with open(peeling, "w") as f:
        for lvl, data in peeling_results.items():
            f.write(
                f"Peeling level {lvl}\n    Number of Protein Units: {data['N']}\n    Compaction Index: {round(data['CI'], 2)}\n"
            )
            for start_pu, end_pu in data["PUs"]:
                if not disable_energies:
                    pu_residues = ",".join(
                        [f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu + 1)]
                    )
                    pu_energy, pu_z_score = get_energy_and_z_score(
                        BIN_DIR, f"{RESULTS_DIR}/{pdb_id_chain}", pu_residues
                    )
                    pu_aul = (
                        int((1 - (1 / (pu_z_score) ** 2)) * 100)
                        if abs(pu_z_score) >= 1
                        else 0
                    )
                    f.write(
                        f"    {str(start_pu)+'-'+str(end_pu):>7}: AUL={pu_aul:3}% Z-score={round(pu_z_score, 1)}\n"
                    )
                else:
                    f.write(f"    {str(start_pu)+'-'+str(end_pu):>7}\n")


def generate_plots(i, part, mat, RESULTS_DIR, pus_colors):
    fig1, ax1 = plt.subplots(figsize=(6, 9), dpi=150)
    ax1.set_xlabel("Residues")
    ax1.set_ylabel("Residues")
    ax1.imshow(mat, cmap="RdPu")
    ax1.invert_yaxis()
    box1 = ax1.get_position()
    ax1.set_position(
        [box1.x0, box1.y0 + box1.height * 0.2, box1.width, box1.height * 0.9]
    )
    if i == 0:
        ax1.set_title(
            "Contact Probability Map of the\noptimal partition (all Protein Units)"
        )
    else:
        ax1.set_title(
            f"Contact Probability Map of the alternative\npartition n°{i} (all Protein Units)"
        )
    for j, domain in part["BOUNDARIES"].items():
        fig2, ax2 = plt.subplots(figsize=(6, 9), dpi=150)
        ax2.set_xlabel("Residues")
        ax2.set_ylabel("Residues")
        ax2.imshow(mat, cmap="RdPu")
        ax2.invert_yaxis()
        box2 = ax2.get_position()
        ax2.set_position(
            [box2.x0, box2.y0 + box2.height * 0.2, box2.width, box2.height * 0.9]
        )
        if i == 0:
            ax2.set_title(
                f"Contact Probability Map of the domain {j+1}\nof the optimal partition"
            )
        else:
            ax2.set_title(
                f"Contact Probability Map of the domain {j+1}\nof the alternative partition n°{i}"
            )
        for start_pu, end_pu in domain:
            fig3, ax3 = plt.subplots(figsize=(5, 6.5), dpi=150)
            ax3.set_xlabel("Residues")
            ax3.set_ylabel("Residues")
            ax3.imshow(mat, cmap="RdPu")
            ax3.invert_yaxis()
            if i == 0:
                ax3.set_title(
                    f"Contact Probability Map of PU {start_pu}-{end_pu} of the domain {j+1}\nof the optimal partition"
                )
            else:
                ax3.set_title(
                    f"Contact Probability Map of PU {start_pu}-{end_pu} of the domain {j+1}\nof the alternative partition n°{i}"
                )
            l = end_pu - start_pu
            rect = patches.Rectangle(
                (start_pu - 1, start_pu - 1),
                l,
                l,
                linewidth=1.5,
                edgecolor="#%02x%02x%02x" % pus_colors[(start_pu, end_pu)],
                facecolor="none",
            )
            rect.set_label(f"{start_pu}-{end_pu}")
            # Set labels on rectangles added to ax1 and ax2
            rect1 = patches.Rectangle(
                (start_pu - 1, start_pu - 1),
                l,
                l,
                linewidth=1.5,
                edgecolor="#%02x%02x%02x" % pus_colors[(start_pu, end_pu)],
                facecolor="none",
                label=f"{start_pu}-{end_pu}",
            )
            rect2 = patches.Rectangle(
                (start_pu - 1, start_pu - 1),
                l,
                l,
                linewidth=1.5,
                edgecolor="#%02x%02x%02x" % pus_colors[(start_pu, end_pu)],
                facecolor="none",
                label=f"{start_pu}-{end_pu}",
            )
            ax1.add_patch(rect1)
            ax2.add_patch(rect2)
            ax3.add_patch(rect)
            ax3.legend(
                title="Protein Unit",
                loc="upper center",
                bbox_to_anchor=(0.5, -0.15),
                fancybox=False,
                shadow=False,
                ncol=3,
                frameon=False
            )
            fig3.savefig(
                os.path.join(
                    RESULTS_DIR,
                    "Contact_Probability_Matrix",
                    f"contact_probability_matrix_alternative_{i}_domain_{j}_pu_{start_pu}_{end_pu}.png",
                ),
                bbox_inches='tight',
            )
            plt.close(fig3)
        ax2.legend(
            title="Protein Units",
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            fancybox=True,
            shadow=False,
            ncol=3,
            frameon=False
        )
        fig2.savefig(
            os.path.join(
                RESULTS_DIR,
                "Contact_Probability_Matrix",
                f"contact_probability_matrix_alternative_{i}_domain_{j}.png",
            ),
            bbox_inches='tight'
        )
        plt.close(fig2)
    ax1.legend(
        title="Protein Units",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        fancybox=True,
        shadow=False,
        ncol=3,
        frameon=False
    )
    fig1.savefig(
        os.path.join(
            RESULTS_DIR,
            "Contact_Probability_Matrix",
            f"contact_probability_matrix_alternative_{i}.png",
        ),
        bbox_inches='tight',
    )
    plt.close(fig1)

if __name__ == "__main__":
    
    start = time.time()

    def check_cpu(nb_cpu):
        """
        Check if the user input CPU nb is valid
        """
        try:
            nb_cpu = int(nb_cpu)
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                "Error option -c/--cpu: please input an integer"
            ) from e
        if 0 <= nb_cpu <= multiprocessing.cpu_count():
            return nb_cpu
        raise argparse.ArgumentTypeError(
            f"Error option -c/--cpu: nb_cpu should be 0 <= nb_cpu <= {multiprocessing.cpu_count()}"
        )

    def check_model(model):
        """
        Check if the user input model nb is valid
        """
        try:
            model = int(model)
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                "Error option -d/--model: please input an integer"
            ) from e
        if model >= 1:
            return model
        raise argparse.ArgumentTypeError(
            "Error option -d/--model: model should be >= 1"
        )

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            SWORD2: SWift and Optimized Recognition of protein Domains.
            The SWORD2 partitioning algorithm produces multiple alternative 
            domain assignments for a given protein structure. 
            This unique approach handles ambiguous protein structure partitioning, 
            admitting several solutions. The decomposition of the protein structure
            into domains is achieved through the hierarchical clustering of Protein Units, 
            evolutionarily preserved structural descriptors at the interface between 
            secondary structures and domains."""
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    optional = parser.add_argument_group("optional arguments")
    required = parser.add_argument_group("required arguments")
    parser.add_argument("--version", action="version", version="SWORD2 2.0.0")
    group.add_argument(
        "-u", "--uniprot-id", help="AlphaFold Uniprot Accession Id.", type=str
    )
    group.add_argument(
        "-m", "--mgnify-id", help="MGnify Id for the ESM Metagenomic Atlas.", type=str
    )
    group.add_argument(
        "-p", "--pdb-id", help="PDB id to download from the PDB database.", type=str
    )
    group.add_argument(
        "-i", "--input-file", help="Path to an input PDB or mmCIF file.", type=str
    )
    optional.add_argument(
        "-c", "--pdb-chain", help="PDB chain. If not specified, the first chain in the PDB file will be used.", type=str, default=None
    )
    optional.add_argument(
        "-d",
        "--model",
        help="Model to parse. Especially useful for NMR files which contain several models. Default is 1.",
        type=check_model,
        default=1,
    )
    optional.add_argument(
        "-x",
        "--cpu",
        help=f"Number of CPUs to use. Default all (0). Max on this computer is: {multiprocessing.cpu_count()}",
        default=0,
        type=check_cpu,
    )
    optional.add_argument(
        "-e",
        "--disable-energies",
        action="store_true",
        help="Disable the calculation of pseudo-energy of domains and PUs.",
    )
    optional.add_argument(
        "-l",
        "--disable-plots",
        action="store_true",
        help="Disable the generation of contact probability matrices plots.",
    )
    required.add_argument(
        "-o",
        "--output",
        help="Output directory. Results will be generated inside in a dedicated directory named after OUTPUT/PDBID_CHAIN/",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    uniprot_id = args.uniprot_id
    mgnify_id = args.mgnify_id
    pdb_id = args.pdb_id
    input_file = args.input_file
    pdb_chain = args.pdb_chain
    model = args.model
    output_dir = args.output
    nb_cpu = args.cpu if args.cpu != 0 else multiprocessing.cpu_count()
    disable_energies = args.disable_energies
    disable_plots = args.disable_plots

    

    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    BIN_DIR = os.path.join(BASE_DIR, "bin")
    SWORD_DIR = os.path.join(BIN_DIR, "SWORD/bin/SWORD")
    SWORD = os.path.join(SWORD_DIR, "SWORD")
    DISPLAY_SWORD2 = os.path.join(BIN_DIR, "display_SWORD2_output.pl")

    # Define a temporary RESULTS_DIR
    TEMP_RESULTS_DIR = output_dir  # Use output_dir for temporary storage

    # Parse and check the PDB, get the prot object and updated pdb_chain
    prot, pdb_chain = check_parsing_pdb(
        uniprot_id, mgnify_id, pdb_id, pdb_chain, model, input_file, TEMP_RESULTS_DIR
    )

    # Now, construct the pdb_id_chain string
    if input_file:
        pdb_id_chain = os.path.basename(os.path.splitext(input_file)[0])
    elif uniprot_id:
        pdb_id_chain = uniprot_id
    elif mgnify_id:
        pdb_id_chain = mgnify_id
    else:
        pdb_id_chain = pdb_id

    pdb_id_chain = pdb_id_chain + "_" + pdb_chain

    RESULTS_DIR = os.path.join(output_dir, pdb_id_chain)
    new_dir = False
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
    else:
        new_dir = True
        name_rep = time.strftime("_%d_%m_%Y_") + "".join(
            random.choice("0123456789") for _ in range(5)
        )
        RESULTS_DIR += name_rep
        os.makedirs(RESULTS_DIR)
    
    if new_dir:
        logging.warning(
            f"Results dir '{os.path.join(output_dir, pdb_id_chain)}' already exists --> New results dir '{RESULTS_DIR}'"
        )
        
    fh = log.FileHandler(os.path.join(RESULTS_DIR, "sword2.log"))
    logging.addHandler(fh)

    # If any files were downloaded to TEMP_RESULTS_DIR, move them to RESULTS_DIR
    if TEMP_RESULTS_DIR != RESULTS_DIR:
        for file_name in os.listdir(TEMP_RESULTS_DIR):
            full_file_name = os.path.join(TEMP_RESULTS_DIR, file_name)
            if os.path.isfile(full_file_name):
                shutil.move(full_file_name, RESULTS_DIR)

    # predict_time(prot) returns the time in seconds
    est_time_in_seconds = None
    if disable_energies and disable_plots:
        est_time_in_seconds = predict_time_no_energy_no_plots(prot)
    else:
        est_time_in_seconds = predict_time_full(prot)

    if est_time_in_seconds < 60:
        est_time_str = f"Estimated runtime: {est_time_in_seconds} seconds"
    else:
        minutes = est_time_in_seconds // 60
        seconds = est_time_in_seconds % 60
        est_time_str = f"Estimated runtime: {minutes} minutes and {seconds} seconds"

    logging.info("")
    logging.info(f">>>   {pdb_id_chain} ({len(set(prot.getResnums()))} aa)")
    logging.info(f">>>   {est_time_str}")
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
        prot = prot.select(
            "protein and not nonstdaa and not hetatm and not resnum " + res_to_remove
        )
    else:
        prot = prot.select("protein and not nonstdaa and not hetatm")
    hv = prot.getHierView()
    seq_resnums = [residue.getResnum() for residue in hv.iterResidues()]
    new_seq_resnums = list(range(1, len(seq_resnums) + 1))
    for i, residue in enumerate(hv.iterResidues()):
        residue.setResnum(new_seq_resnums[i])
    file_written = writePDB(pdb_chain_file, prot)
    if file_written is None:
        logging.warning(
            f"PDB file {pdb_chain_file} could not be written using function writePDB of ProDy"
        )
        sys.exit(1)

    # Remove the ".pdb" extension for SWORD
    os.rename(file_written, os.path.splitext(file_written)[0])

    # Now that the protein is clean, get some infos
    hv = prot.getHierView()
    ch = hv.getChain(pdb_chain)

    # Get the aa sequence
    seq = ch.getSequence()
    prot_len = len(seq)
    seq_resnums_norm = {num: i for i, num in enumerate(seq_resnums)}
    seq_resnums_norm_inv = {i: num for i, num in enumerate(seq_resnums)}

    logging.info("Launch SWORD")
    if not os.path.exists(SWORD_DIR + "/bin/Dssp/dsspcmbi"):
        subprocess.run(SWORD, capture_output=True)
    cmd_args = (
        f"{DISPLAY_SWORD2} '{SWORD} -i {pdb_id_chain} --dir {RESULTS_DIR} -max 9 -nbcpu {nb_cpu}'"
    )
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    with open(f"{RESULTS_DIR}/sword.err", "w") as f:
        f.write(output.stderr.decode("utf-8") + "\n")
    output = output.stdout.decode("utf-8").split("\n")
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

    if not disable_energies:
        logging.info("Calculate pseudo-energies of Domains")
        manager = multiprocessing.Manager()
        energies = manager.dict()
        for i, part in sword_results["DOMAINS"].items():
            with multiprocessing.Pool(processes=nb_cpu) as p:
                func = partial(
                    multiprocess_get_energy,
                    i,
                    pdb_chain,
                    pdb_id_chain,
                    RESULTS_DIR,
                    BIN_DIR,
                    energies,
                )
                p.imap_unordered(func, list(part["BOUNDARIES"].items()))
                p.close()
                p.join()
    else:
        energies = {}

    ############################
    # Write  SWORD partitionings
    ############################

    write_partitionings(sword_results, energies, disable_energies)
    write_partitionings_json(sword_results, energies, disable_energies)

    #######################
    # Write Peeling results
    #######################

    write_peeling_results(disable_energies)

    ############
    # Set colors
    ############

    pus_colors, dom_colors = define_colors(sword_results)

    ######################################################
    # Write Javascript formated histogram of SWORD Domains
    ######################################################

    write_domains_histogram(sword_results, dom_colors)

    #############
    # Contact map
    #############

    if not disable_plots:
        font = {"weight": "normal", "size": 11}
        plt.rc("font", **font)
        plt.rcParams["axes.linewidth"] = 0.5
        plt.rcParams["xtick.major.size"] = 1.5
        plt.rcParams["ytick.major.size"] = 1.5
        plt.rcParams["figure.max_open_warning"] = 0

        os.makedirs(os.path.join(RESULTS_DIR, "Contact_Probability_Matrix"), exist_ok=True)
        proba_mat_file = os.path.join(
            RESULTS_DIR, "PDBs_Clean", pdb_id_chain, "file_proba_contact.mat"
        )
        mat = np.loadtxt(proba_mat_file)

        # Use multiprocessing to parallelize plot generation
        logging.info("Generate contact probability matrices")
        with multiprocessing.Pool(processes=nb_cpu) as pool:
            func = partial(generate_plots, mat=mat, RESULTS_DIR=RESULTS_DIR, pus_colors=pus_colors)
            pool.starmap(func, sword_results["DOMAINS"].items())
    
    #########################
    # Junctions consistencies
    #########################

    logging.info("Calculate junctions consistencies")

    cmd_args = f"{BIN_DIR}/stat_pu_domains_from_SWORD.pl {RESULTS_DIR}/sword.txt"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8").split("\n")
    junctions = {}
    with open(f"{RESULTS_DIR}/junctions_consistencies.txt", "w") as f:
        for line in output:
            jctn_found = re.search(r"^(\d+)\s+(\d+)\s+(\d\.\d+)\s+(\d+\.\d+)$", line)
            if jctn_found:
                jct = int(jctn_found.group(1))
                cnt = int(jctn_found.group(2))
                raw = float(jctn_found.group(3))
                weight = float(jctn_found.group(4))
                junctions[jct] = {"cnt": cnt, "raw": raw, "weight": weight}
            f.write(line + "\n")
        del junctions[min(junctions.keys())]
        del junctions[max(junctions.keys())]

    # Mapping of authors PDB residues numbers with the new one presented on the server
    with open(f"{RESULTS_DIR}/mapping_auth_resnums.txt", "w") as f1:
        f1.write(
            "# Mapping of authors PDB residues numbers \n# with the new one presented on the server\nORIGINAL RENUM\n"
        )
        for i, j in enumerate(seq_resnums, start=1):
            f1.write(f"{j} {i}\n")

    logging.info("Clean and prepare results")
    shutil.rmtree(
        os.path.join(RESULTS_DIR, "PDBs_Stand"), ignore_errors=False, onerror=None
    )
    shutil.move(
        os.path.join(RESULTS_DIR, "PDBs_Clean"), os.path.join(RESULTS_DIR, "SWORD")
    )
    os.makedirs(os.path.join(RESULTS_DIR, "Junctions"), exist_ok=True)
    shutil.move(
        os.path.join(RESULTS_DIR, "junctions_consistencies.txt"),
        os.path.join(RESULTS_DIR, "Junctions"),
    )
    os.system(f"mv {RESULTS_DIR}/SWORD/*/Peeling {RESULTS_DIR}/Protein_Units")
    logging.info(f"Results can be found here: {RESULTS_DIR}")

    logging.info(f"Total runtime: {int(time.time()-start)} seconds")
