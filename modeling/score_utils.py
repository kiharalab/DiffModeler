import os

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
    listpdb = [x for x in os.listdir(pdb_dir) if ".pdb" in x]
    listpdb.sort()
    list_score=[]
    with open(output_path,'r') as file:
        for line in file:
            if line.startswith("LDP Recall Score"):
                score=float(line.strip("\n").replace("LDP Recall Score:",""))
                list_score.append(score)

    assert len(list_score)==len(listpdb)
    for kk in range(len(listpdb)):
        pdb_path = os.path.join(pdb_dir,listpdb[kk])
        score= list_score[kk]
        new_score_dict[pdb_path]=score
    new_score_dict=sort_dict_by_value_desc(new_score_dict)
    return new_score_dict
