import os
from ops.os_operation import mkdir
from modeling.score_utils import sort_dict_by_value_desc,clean_score_dict,find_biggest_unvisited_chain
from modeling.pdb_utils import find_identical_chain,find_chain_pdb_path,remove_overlap_pdb,rename_chains_cif,calculate_overlap_score
from modeling.fit_structure_chain import fit_single_chain
from modeling.map_utils import mask_map_by_pdb,segment_map,format_map
from modeling.pdb_utils import extract_residue_locations
import shutil
from ops.io_utils import write_pickle
def iterative_fitting(diff_trace_map,diff_ldpmap_path,
                  modeling_dir,score_dict,fitting_dict,
                  map_ldp_pdb_path,chain_visit_dict,chain_length_score,params):
    print("iterative fitting, still remains:")
    print(chain_visit_dict)
    if len(score_dict)>0:
        score_dict=sort_dict_by_value_desc(score_dict)
        search_key = list(score_dict.keys())[0]
        current_key= search_key
        current_chain = current_key.split(",")[0]
        current_score = score_dict[search_key]
        if chain_visit_dict[current_chain]==1:
            #if the top scored chain have been fitted, check if it has identical chains that has not been fitted
            current_chain_list = find_identical_chain(fitting_dict,current_chain,chain_visit_dict)
            if current_chain_list is None or len(current_chain_list)==0:
                #no need to consider this fit since it has no identical chains and it has been fitted
                score_dict = clean_score_dict(score_dict,current_chain)
                iterative_fitting(diff_trace_map,diff_ldpmap_path,
                  modeling_dir,score_dict,fitting_dict,
                  map_ldp_pdb_path,chain_visit_dict,chain_length_score,params)
                return
            current_assign_chain = current_chain_list[0]#use its identical chains to assign chain name
        else:
            current_assign_chain = current_chain
        current_fitpdb_path =  current_key.split(",")[1]
    else:
        #if score pool has become empty, consider remained chains to fit the largest for global fitting in remained mask map
        #pick the largetest
        select_chain = find_biggest_unvisited_chain(chain_length_score,chain_visit_dict)
        current_score = chain_length_score[select_chain]
        current_fitpdb_path = find_chain_pdb_path(fitting_dict,select_chain)
        current_chain = select_chain
        current_assign_chain = current_chain

    global_refit_flag=False
    if current_score<params["assembling"]['score_min_cutoff']:
        #do global refitting
        # if score dict is empty will for sure global fitting since sequence score max score is 5
        #do global refitting
        current_output_dir = os.path.join(modeling_dir,"globalsearch_%s"%current_assign_chain)
        mkdir(current_output_dir)
        print("score %.2f: 1st global refinment!"%(current_score))
        new_score_dict=fit_single_chain(diff_ldpmap_path,current_fitpdb_path,
                                        current_output_dir,map_ldp_pdb_path,params,global_mode=1)
        #clean the global refitting to remove any structure
        unclean_global_score_dict = new_score_dict
        #check all fitting  structure and clean overlap
        #in principle should no such overlap cases, but may happen in some severe conditions
        clash_distance= params['assembling']['clash_distance']
        ratio_cutoff = params['assembling']['overlap_ratio_limit']
        for key in chain_visit_dict:
            if chain_visit_dict[key]==1:
                previous_fit_dir = os.path.join(modeling_dir,"iterative_%s"%key)
                previous_fit_pdb = os.path.join(previous_fit_dir,"final.pdb")
                new_score_dict = remove_overlap_pdb(new_score_dict,previous_fit_pdb,clash_distance,ratio_cutoff,split_key=False)
                if len(new_score_dict)<=0:
                    break
        print("global fitting cleaning remains %d/%d"%(len(new_score_dict),len(unclean_global_score_dict)))
        if len(new_score_dict)>0:
            new_score_dict=sort_dict_by_value_desc(new_score_dict)#if still remain satisfied, pick the remained top1
            current_file_name=list(new_score_dict.keys())[0]
        else:
            new_score_dict=unclean_global_score_dict
            clash_dict={}
            #calculate all the clashes and sum up them, pick the one with lowest sum clash score
            for key in chain_visit_dict:
                if chain_visit_dict[key]==1:
                    previous_fit_dir = os.path.join(modeling_dir,"iterative_%s"%key)
                    previous_fit_pdb = os.path.join(previous_fit_dir,"final.pdb")
                    clash_dict=calculate_overlap_score(new_score_dict,clash_dict,previous_fit_pdb,clash_distance,split_key=False)
            clash_dict=sort_dict_by_value_desc(clash_dict)
            current_file_name=list(clash_dict.keys())[-1]
            print("repicking, finding the mininum clash one with %f"%clash_dict[current_file_name])
        current_key = current_file_name#it is actually a path
        current_fitpdb_path = os.path.join(current_output_dir,"top1.pdb")
        shutil.copy(current_file_name,current_fitpdb_path)
        new_score = new_score_dict[current_file_name]+chain_length_score[current_assign_chain]
        print("after global fitting score %.2f compared to previous %.2f"%(new_score,current_score))
        current_score=new_score
        global_refit_flag=True
    #local fitting
    current_output_dir = os.path.join(modeling_dir,"iterative_%s"%current_assign_chain)
    mkdir(current_output_dir)
    print("score %.2f: local refinment!"%(current_score))
    segment_map_path0 = os.path.join(current_output_dir,"segment_map_tmp.mrc")
    mask_map_by_pdb(diff_trace_map,segment_map_path0,current_fitpdb_path,keep_label=True)
    segment_map_path0 = format_map(segment_map_path0,os.path.join(current_output_dir,"segment_map.mrc"))
    segment_map_path = os.path.join(current_output_dir,"segment_map_small.mrc")
    segment_map(segment_map_path0,segment_map_path)

    # we skipped chimera local release in open code to remove the dependency to chimera
    final_pdb_output = os.path.join(current_output_dir,"final.pdb")
    if not global_refit_flag:
        #for local fitting using map
        new_score_dict = fit_single_chain(segment_map_path,current_fitpdb_path,
                                          current_output_dir,map_ldp_pdb_path,params,global_mode=0)
        current_file_name=list(new_score_dict.keys())[0]
        new_score = new_score_dict[current_file_name]
        new_score = new_score+chain_length_score[current_assign_chain]
        print("finally selected %s score:%.2f"%(current_file_name,new_score_dict[current_file_name]))
        print("before vesper local fitting score: %.2f"%current_score)
        print("after vesper local fitting score: %.2f"%new_score)
        if new_score>=current_score:
            shutil.copy(current_file_name,final_pdb_output)
        else:
            print("still used the begining top structure %s"%current_fitpdb_path)
            print("Score: %.2f"%current_score)
            shutil.copy(current_fitpdb_path,final_pdb_output)
    else:
        #use global refit one
        shutil.copy(current_fitpdb_path,final_pdb_output)

    #mask map and map ldp according to fitted structure.
    diff_traced_map_new = os.path.join(current_output_dir,"iterative_%s.mrc"%current_assign_chain)
    diff_traced_map_tmp = os.path.join(current_output_dir,"iterative_%s_tmp.mrc"%current_assign_chain)
    mask_map_by_pdb(diff_trace_map,diff_traced_map_tmp,final_pdb_output,keep_label=False)
    format_map(diff_traced_map_tmp,diff_traced_map_new)

    diff_ldpmap_path_new = os.path.join(current_output_dir,"iterativeLDP_%s.mrc"%current_assign_chain)
    diff_ldpmap_path_tmp = os.path.join(current_output_dir,"iterativeLDP_%s.mrc"%current_assign_chain)
    mask_map_by_pdb(diff_ldpmap_path,diff_ldpmap_path_tmp,final_pdb_output,keep_label=False)
    format_map(diff_ldpmap_path_tmp,diff_ldpmap_path_new)

    clash_distance= params['assembling']['clash_distance']
    ratio_cutoff = params['assembling']['overlap_ratio_limit']
    score_dict = remove_overlap_pdb(score_dict,final_pdb_output,clash_distance,ratio_cutoff)

    score_path_new = os.path.join(current_output_dir,"remain_score.pkl"%score_dict)
    write_pickle(score_dict,score_path_new)
    print("initial candidates still remained %d"%(len(score_dict)))
    #adjust chain_visit_dict
    chain_visit_dict[current_assign_chain]=1
    chain_remains =[]
    for key in chain_visit_dict:
        if chain_visit_dict[key]==0:
            chain_remains.append(key)
    print("still remains %d chains"%len(chain_remains))

    finalnew_pdb_path = os.path.join(current_output_dir,"final_rechain.cif")
    rename_chains_cif(final_pdb_output,current_assign_chain,finalnew_pdb_path)



    if len(chain_remains)>0:
        iterative_fitting(diff_traced_map_new,diff_ldpmap_path_new,
                  modeling_dir,score_dict,fitting_dict,
                  map_ldp_pdb_path,chain_visit_dict,chain_length_score,params)




