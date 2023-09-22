import os
from ops.os_operation import mkdir,functToDeleteItems
from ops.io_utils import download_file
from ops.pdb_utils import count_atom_line,filter_chain_cif,cif2pdb,filter_chain_pdb,count_residues
from ops.fasta_utils import read_fasta,write_all_fasta
from ops.io_utils import write_pickle,load_pickle
from collections import defaultdict
def parse_blast_output(blast_file):
    match_dict=defaultdict(list)
    read_flag=False
    with open(blast_file,'r') as rfile:
        all_lines = rfile.readlines()
        for line_index,line in enumerate(all_lines):
            if line.startswith("Query="):
                line =line.strip("\n")
                line = line.replace("Query=","")
                current_id = line.replace(" ","")
            elif line.startswith("Sequences producing significant"):
                read_flag=True
                continue
            if read_flag:
                line=line.strip("\n")
                if len(line)>0:
                    line = line.split()
                    match_id = line[0]
                    score = float(line[-1])
                    match_dict[current_id].append([match_id,score])
                    for k in range(1,10):
                        new_line=all_lines[k+line_index]
                        new_line = new_line.strip("\n").split()
                        match_id = new_line[0]
                        score = float(new_line[-1])
                        match_dict[current_id].append([match_id,score])
                    read_flag=False
    return match_dict

def download_pdb(pdb_id,current_chain_dir,final_pdb_path):
    pdb = pdb_id.split("_")[0]
    chain_id = pdb_id.split("_")[1]
    download_link = "https://files.rcsb.org/download/%s.cif"%pdb
    cif_file = os.path.join(current_chain_dir,"input.cif")
    download_file(download_link,cif_file)
    if os.path.exists(cif_file) and count_atom_line(cif_file)>=50:
        #segment the specific chains
        chain_cif = os.path.join(current_chain_dir,"input_%s.cif"%chain_id)
        filter_chain_cif(cif_file,chain_id,chain_cif)

        #then convert cif file format to pdb
        cif2pdb(chain_cif,final_pdb_path)
    else:
        #download the pdb if the old one did not exist
        download_link = "https://files.rcsb.org/download/%s.pdb"%pdb
        pdb_file = os.path.join(current_chain_dir,"input.pdb")
        download_file(download_link,pdb_file)
        filter_chain_pdb(pdb_file,chain_id,final_pdb_path)
    return final_pdb_path
def fasta_searchdb(params,save_path):
    single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
    mkdir(single_chain_pdb_dir)
    fitting_pickle_path = os.path.join(save_path,"fitting_dict.pkl")
    if os.path.exists(fitting_pickle_path) and os.path.getsize(fitting_pickle_path)>10:
        fitting_dict = load_pickle(fitting_pickle_path)
        return fitting_dict
    #reorganize data
    fasta_path = os.path.abspath(params['P'])
    chain_dict = read_fasta(fasta_path)
    fasta_path = os.path.join(single_chain_pdb_dir,"input.fasta")
    write_all_fasta(chain_dict,fasta_path)

    #first do experimental database search
    output_path = os.path.join(single_chain_pdb_dir,"exp_search.out")
    search_command = "blastp -query %s -db %s -out %s -num_threads %d"%(fasta_path,params['db_exp_path'],
                                                                        output_path, params['search_thread'])
    os.system(search_command)

    matched_dict = {}#[key]: chain id list, [value]: the structure id
    fitting_dict={}
    #parse the information
    exp_match_dict = parse_blast_output(output_path)
    #merge only identical chains, otherwise, rely on combine search
    for key in exp_match_dict:
        current_match_list = exp_match_dict[key]
        for k in range(len(current_match_list)):
            match_id, evalue = current_match_list[k]
            if evalue==0:
                #fetch the pdb to see if the protein size really reasonable
                expected_seq_length = len(chain_dict[key])*params['search']['length_ratio']
                chain_name_list = key.replace(",","-")
                current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list))
                mkdir(current_chain_dir)
                final_pdb_path = os.path.join(single_chain_pdb_dir,chain_name_list+".pdb")
                download_pdb(match_id,current_chain_dir,final_pdb_path)
                actual_structure_length = count_residues(final_pdb_path)
                if actual_structure_length>=expected_seq_length and actual_structure_length<=len(chain_dict[key]):
                    matched_dict[key]="PDB:"+match_id
                    final_chain_list = chain_name_list.split("-")
                    fitting_dict[final_pdb_path]=final_chain_list
                    break
                else:
                    os.remove(final_pdb_path)
                    functToDeleteItems(current_chain_dir)
            if evalue>0:
                break
    matched_keys = matched_dict.keys()
    #write a new fasta to search

    remain_chain_dict = {k:v for k,v in chain_dict.items() if k not in matched_keys}
    print("experimental db search finished, get %s/%s matched single chain structure"%(len(matched_keys),len(chain_dict)))
    if len(remain_chain_dict)>0:
        print("continue PDB+AFDB search for remained %d chains"%len(remain_chain_dict))
        remain_fasta_path = os.path.join(single_chain_pdb_dir,"remain_search.fasta")
        write_all_fasta(remain_chain_dict,remain_fasta_path)
        output_path = os.path.join(single_chain_pdb_dir,"expaf_search.out")
        search_command = "blastp -query %s -db %s -out %s -num_threads %d"%(remain_fasta_path,params['db_path'],
                                                                            output_path, params['search_thread'])
        os.system(search_command)
        expaf_match_dict = parse_blast_output(output_path)
        if len(expaf_match_dict)!=len(remain_chain_dict):
            for tmp_key in remain_chain_dict:
                if tmp_key not in expaf_match_dict:
                    print("warning, query %s can not find any templates in database"%tmp_key)

        for key in expaf_match_dict:
            current_match_list = exp_match_dict[key]
            for k in range(len(current_match_list)):
                match_id, evalue = current_match_list[k]
                if "AFDB" in match_id:
                    matched_dict[key]=match_id
                    break
                else:
                    expected_seq_length = len(chain_dict[key])*params['search']['length_ratio']
                    chain_name_list = key.replace(",","-")
                    current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list))
                    mkdir(current_chain_dir)
                    final_pdb_path = os.path.join(single_chain_pdb_dir,chain_name_list+".pdb")
                    download_pdb(match_id,current_chain_dir,final_pdb_path)
                    actual_structure_length = count_residues(final_pdb_path)
                    if actual_structure_length>=expected_seq_length  and actual_structure_length<=len(chain_dict[key]):
                        matched_dict[key]="PDB:"+match_id
                        final_chain_list = chain_name_list.split("-")
                        fitting_dict[final_pdb_path]=final_chain_list
                        break
                    else:
                        os.remove(final_pdb_path)
                        functToDeleteItems(current_chain_dir)

    print("DB search finished! Match relationship ",matched_dict)
    #get the matched dict


    for chain_name_list in chain_dict:
        if chain_name_list not in matched_dict:
            print("*"*100)
            print("WARNING! query %s can not find any templates in database"%chain_name_list)
            print("*"*100)
            continue
        matched_id = matched_dict[chain_name_list]
        chain_name_list = chain_name_list.replace(",","-")
        current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list))
        mkdir(current_chain_dir)
        final_pdb_path = os.path.join(single_chain_pdb_dir,chain_name_list+".pdb")
        if final_pdb_path in fitting_dict:
            continue
        if os.path.exists(final_pdb_path) and count_atom_line(final_pdb_path)>=50:
            final_chain_list = chain_name_list.split("-")
            fitting_dict[final_pdb_path]=final_chain_list
            continue
        split_info = matched_id.split(":")
        database = split_info[0]
        pdb_id = split_info[1]
        if database=="AFDB":
            #alphafold db
            download_link = "https://alphafold.ebi.ac.uk/files/%s-model_v4.pdb"%pdb_id
            download_file(download_link,final_pdb_path)
        else:
            download_pdb(pdb_id,current_chain_dir,final_pdb_path)
        final_chain_list = chain_name_list.split("-")
        fitting_dict[final_pdb_path]=final_chain_list
    print("collecting finish: fitting dict: ",fitting_dict)
    write_pickle(fitting_dict,fitting_pickle_path)
    return fitting_dict
