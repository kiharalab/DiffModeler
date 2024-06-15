
import os
from ops.io_utils import write_pickle,load_pickle
from modeling.score_utils import build_score_pool,filter_score_dict,add_structure_size_score
from modeling.pdb_utils import collect_final_pdb
from modeling.iterative_fitting import iterative_fitting
from ops.os_operation import mkdir
def  assemble_structure(diff_trace_map,fitting_dict,fitting_dir,modeling_dir,params):
    mkdir(modeling_dir)
    #first obtain the ldp path
    bandwidth = 2
    output_pdb_dir = os.path.join(fitting_dir,"LDP_%s"%(bandwidth))
    map_ldp_pdb_path = os.path.join(output_pdb_dir,"chain_LDP.pdb")
    diff_ldpmap_path = os.path.join(output_pdb_dir,"chain_LDP.mrc")
    if not os.path.exists(map_ldp_pdb_path):
        print("no LDP has been generated in fitting stage!Failed!")
        return

    overall_score_pickle = os.path.join(modeling_dir,"Initial_Score.pkl")
    if not os.path.exists(overall_score_pickle):
        overall_score_dict = build_score_pool(fitting_dir,fitting_dict)
        write_pickle(overall_score_dict,overall_score_pickle)
    else:
        overall_score_dict = load_pickle(overall_score_pickle)
    print("Initial Score Pool Constructed")

    overall_score_list = [overall_score_dict[k] for k in overall_score_dict]
    overall_score_list = sorted(overall_score_list,reverse=True)
    score_cutoff= overall_score_list[int(len(overall_score_list)*params['assembling']['score_cutoff'])]
    score_cutoff= min(score_cutoff,params['assembling']['score_max_cutoff'])
    score_cutoff = max(score_cutoff,params['assembling']['score_min_cutoff'])
    print("current acceptable score is %.2f"%score_cutoff)
    overall_score_dict = filter_score_dict(overall_score_dict,score_cutoff)
    print("%d/%d remained to assemble structures"%(len(overall_score_dict),len(overall_score_list)))

    overall_score_dict,chain_length_score =  add_structure_size_score(fitting_dict,fitting_dir,overall_score_dict)

    chain_visit_dict= {}
    for pdb_path in fitting_dict:
        chain_list = fitting_dict[pdb_path]
        for current_chain in chain_list:
            chain_visit_dict[current_chain]=0
    iterative_fitting(diff_trace_map,diff_ldpmap_path,
                  modeling_dir,overall_score_dict,fitting_dict,
                  map_ldp_pdb_path,chain_visit_dict,chain_length_score,params)

    final_path=collect_final_pdb(modeling_dir,chain_visit_dict)
    if params['domain']:

        #also generate a cif with original chain id instead of domain chain ids [chainid_domainid] format
        from modeling.pdb_utils import collect_domain_pdb
        final_new_path = collect_domain_pdb(modeling_dir,chain_visit_dict)
        print("The domain based chain name cif is saved at %s"%final_path)
        final_path = final_new_path 
    return final_path
