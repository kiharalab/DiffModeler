import os
import shutil
from ops.io_utils import load_pickle
from modeling.pdb_utils import extract_residue_locations

def sort_dict_by_value_desc(d: dict) -> dict:
    """
    Given a dictionary `d`, returns a new dictionary with the same keys as `d`,
    but with the values sorted in descending order.
    """
    sorted_items = sorted(d.items(), key=lambda x: x[1], reverse=True)
    return {k: v for k, v in sorted_items}
# def evaluate_chain_trace(diffmap_ldp_path,current_pdb_path):
#     #recall score is integrated into VESPER
#
#
# def calculate_candidates_score(current_output_dir,diffmap_ldp_path):
#
#     listfiles = [x for x in os.listdir(current_output_dir) if "vesper" in x and ".pdb" in x
#                  and "output" not in x and ".txt" not in x]
#     assert  len(listfiles)>0#avoid running those cases structure have not finished run
#     score_dict={}
#     check_flag=False
#     for pdb_entry in listfiles:
#         current_pdb_path = os.path.join(current_output_dir,pdb_entry)
#         try:
#             bb_recall = evaluate_chain_trace(diffmap_ldp_path,current_pdb_path)
#             check_flag =True
#         except:
#             bb_recall = -1
#         print("eval top: %s %.2f"%(pdb_entry,bb_recall))
#         score_dict[pdb_entry]=bb_recall
#
#     if check_flag is False:
#         exit()#no reasonable fitting exists, not reasonable outcome
#     return score_dict

def read_score(new_score_dict,pdb_dir,output_path):
    pdb_dir = os.path.abspath(pdb_dir)
    listoldpdb = [x for x in os.listdir(pdb_dir) if ".pdb" in x]
    listoldpdb.sort()
    for item in listoldpdb:
        old_pdb_path = os.path.join(pdb_dir,item)
        key_id = item.split("_")[0]
        new_pdb_path = os.path.join(pdb_dir,key_id+".pdb")
        shutil.move(old_pdb_path,new_pdb_path)

    list_score=[]
    check_flag=False
    with open(output_path,'r') as file:
        for line in file:
            if line.startswith("#0") and not check_flag:
                check_flag=True
            if check_flag and line.startswith("LDP Recall Score"):
                score=float(line.strip("\n").replace("LDP Recall Score:",""))
                list_score.append(score)

    assert len(list_score)==len(listoldpdb)
    for kk in range(len(list_score)):
        pdb_path = os.path.join(pdb_dir,"#%d.pdb"%kk)
        score= list_score[kk]
        new_score_dict[pdb_path]=score
    new_score_dict=sort_dict_by_value_desc(new_score_dict)
    return new_score_dict


def build_score_pool(fitting_dir,fitting_dict):
    overall_score_dict = {}
    for pdb_path in fitting_dict:
        chain_list = fitting_dict[pdb_path]
        for current_chain in chain_list:
            cur_fit_dir = os.path.join(fitting_dir,current_chain)
            score_path = os.path.join(cur_fit_dir,"score.pkl")
            if not os.path.exists(score_path):
                print("%s score is not calculated"%cur_fit_dir)
                exit()
            current_score_info = load_pickle(score_path)
            for path_id in current_score_info:
                overall_score_dict[current_chain+","+path_id]=current_score_info[path_id]
    return overall_score_dict


def filter_score_dict(overall_score_dict,score_cutoff):
    new_score={}
    for key in overall_score_dict:
        cur_score = overall_score_dict[key]
        if cur_score>=score_cutoff:
            new_score[key]=overall_score_dict[key]
    return new_score


def add_structure_size_score(fitting_dict,fitting_dir,overall_score_dict):
    #get the chain length dict and then get a relative score, with range of 5,minimal 0, max 5
    chain_length_dict={}
    for pdb_path in fitting_dict:
        chain_list = fitting_dict[pdb_path]
        for current_chain in chain_list:
            cur_fit_dir = os.path.join(fitting_dir,current_chain)

            current_path = os.path.join(cur_fit_dir,"top1.pdb")
            res_locations=extract_residue_locations(current_path)
            chain_length_dict[current_chain]=len(res_locations)
    chain_length_score ={}
    chain_length_list = list(chain_length_dict.values())
    min_chain_length=min(chain_length_list)
    max_chain_length=max(chain_length_list)
    for pdb_path in fitting_dict:
        chain_list = fitting_dict[pdb_path]
        for chain in chain_list:
            if max_chain_length==min_chain_length:
                chain_length_score[chain]=0
            else:
                chain_length_score[chain]=5*(chain_length_dict[chain]-min_chain_length)/(max_chain_length-min_chain_length)

    #recalculate overall_score_dict
    for key in overall_score_dict:
        current_chain = key.split(",")[0]
        overall_score_dict[key]+=chain_length_score[current_chain]
    return overall_score_dict,chain_length_score


def clean_score_dict(score_dict,chain_name):
    """
    remove records of specified chain in score_dict
    :param score_dict:
    :param chain_name:
    :return:
    """
    new_score_dict={}
    for tmp_key in score_dict:
        tmp_chain= tmp_key.split(",")[0]
        if tmp_chain==chain_name:
            continue
        new_score_dict[tmp_key]=score_dict[tmp_key]
    score_dict=new_score_dict
    return score_dict

def find_biggest_unvisited_chain(chain_length_score,chain_visit_dict):
    """
    find the chains that are biggest that have not been fitted.
    :param chain_length_score:
    :param chain_visit_dict:
    :return:
    """
    chain_length_score=sort_dict_by_value_desc(chain_length_score)
    for key in chain_length_score:#first find biggest that has not visited
        if chain_visit_dict[key]==0:
            select_chain=key
            break
    return select_chain
