from ops.os_operation import mkdir,functToDeleteItems
from modeling.pdb_utils import filter_backbone,generate_move_structure
import os
from graph.build_ldp import build_ldp
from ops.map_utils import SimuMap1A
import shutil
from ops.io_utils import load_pickle,write_pickle
from modeling.score_utils import read_score
from modeling.map_utils import mask_map_by_pdb


def fit_structure_chain(input_map_path,fitting_dict,fitting_dir,params):
    mkdir(fitting_dir)
    bandwidth_list=[1,2]
    fit_target_map_list=[]
    for bandwidth in bandwidth_list:
        output_pdb_dir = os.path.join(fitting_dir,"LDP_%s"%(bandwidth))
        params['LDP']['g']=bandwidth
        map_ldp_pdb_path = os.path.join(output_pdb_dir,"chain_LDP.pdb")
        if not os.path.exists(map_ldp_pdb_path):
            build_ldp(input_map_path,output_pdb_dir,params['LDP'])
        fit_target_map_list.append(os.path.join(output_pdb_dir,"chain_LDP.mrc"))
        verify_pdb_path = map_ldp_pdb_path
    root_fit_target_list = fit_target_map_list

    for pdb_path in fitting_dict:
        chain_list = fitting_dict[pdb_path]
        assert len(chain_list)>1
        tmp_chain = chain_list[0]
        backbone_pdb = os.path.join(fitting_dir,"backbone_%s.pdb"%tmp_chain)
        filter_backbone(pdb_path,backbone_pdb)
        fit_target_map_list = root_fit_target_list
        #generate simulated map from the pdb
        gen_reso_map_path = os.path.join(fitting_dir,"simumap_backbone_1A_%s.mrc"%tmp_chain)
        SimuMap1A(backbone_pdb,gen_reso_map_path)
        vesper_script = os.path.join(os.getcwd(),"VESPER_CUDA")
        vesper_script = os.path.join(vesper_script,"main.py")
        if not os.path.exists(vesper_script):
            print("missing vesper running script %s"%vesper_script)
            exit()
        for current_chain in chain_list:
            cur_fit_dir = os.path.join(fitting_dir,current_chain)
            mkdir(cur_fit_dir)

            score_path = os.path.join(cur_fit_dir,"score.pkl")
            # listfiles = [x for x in os.listdir(cur_fit_dir) if "vesper" in x and ".pdb" in x
            #      and "output" not in x and ".txt" not in x]

            if not os.path.exists(score_path):
                functToDeleteItems(cur_fit_dir)
                #count_model=0
                new_score_dict={}
                for k,fit_target_map in enumerate(fit_target_map_list):
                    #bandwidth = bandwidth_list[k]
                    cur_fit_experiment_path = os.path.join(cur_fit_dir,"fit_experiment_%d"%k)
                    mkdir(cur_fit_experiment_path)
                    output_path = os.path.join(cur_fit_dir,"vesper_simu_output_%d.out"%k)
                    command_line="python3 %s orig -a %s -t %f -b %s -T %f -g %f -s %f " \
                                 "-A %f -N %d -M %s -gpu 0 -o %s -pdbin %s -ca %s -ldp %s -c %d >%s"\
                                 %(vesper_script,fit_target_map,params['vesper']['ldp_cutoff'],gen_reso_map_path,
                                 params['vesper']['simu_cutoff'],params['vesper']['kernel_size'],
                                   params['vesper']['voxel_spacing'],params['vesper']['angle_spacing'],
                                   params['vesper']['num_models'],params['vesper']['rank_mode'],cur_fit_experiment_path,
                                   pdb_path,backbone_pdb,verify_pdb_path,params['vesper']['thread'],output_path
                                   )
                    os.system(command_line)
                    pdb_dir=os.path.join(cur_fit_experiment_path,"PDB")
                    read_score(new_score_dict,pdb_dir,output_path)
                    #count_model = generate_move_structure(pdb_path,cur_fit_dir,output_path,count_model)
                #get all generated vesper model to calculate LDP correlation
                # new_score_dict=calculate_candidates_score(cur_fit_dir,verify_pdb_path)
                write_pickle(new_score_dict,score_path)
                #use the top model to mask the corresponding regions to continue fitting
                current_file_name=list(new_score_dict.keys())[0]
                current_fitpdb_path = os.path.join(cur_fit_dir,current_file_name)
                new_score = new_score_dict[current_file_name]
                print("top global fitting score %.2f "%(new_score))
                top_fitpdb_path =  os.path.join(cur_fit_dir,"top1.pdb")
                shutil.copy(current_fitpdb_path,top_fitpdb_path)
            else:
                new_score_dict = load_pickle(score_path)
                current_file_name=list(new_score_dict.keys())[0]
                current_fitpdb_path = os.path.join(cur_fit_dir,current_file_name)
                new_score = new_score_dict[current_file_name]
                print("top global fitting score %.2f "%(new_score))
                top_fitpdb_path =  os.path.join(cur_fit_dir,"top1.pdb")
                shutil.copy(current_fitpdb_path,top_fitpdb_path)

            #apply mask
            new_fit_target_map_list=[]
            for k,fit_target_map in enumerate(fit_target_map_list):
                diff_map_path_new = os.path.join(cur_fit_dir,"iterative_%d.mrc"%k)
                if not os.path.exists(diff_map_path_new):
                    mask_map_by_pdb(fit_target_map,diff_map_path_new,top_fitpdb_path,keep_label=False)
                new_fit_target_map_list.append(diff_map_path_new)
            fit_target_map_list = new_fit_target_map_list

def fit_single_chain(input_map_path,input_pdb_path,output_dir,ldp_pdb_path,params,global_mode=0):
    mkdir(output_dir)
    #functToDeleteItems(output_dir)
    #generate backbone from pdb
    backbone_pdb = os.path.join(output_dir,"backbone.pdb")
    filter_backbone(input_pdb_path,backbone_pdb)

    gen_reso_map_path = os.path.join(output_dir,"simumap_backbone_1A.mrc")
    SimuMap1A(backbone_pdb,gen_reso_map_path)

    kernel_size=2
    voxel_spacing =  2

    if global_mode==1:
        output_pdb_path = os.path.join(output_dir,"vesper_globalfit.out")
        angle_spacing = 10
    else:
        output_pdb_path = os.path.join(output_dir,"vesper_localfit.out")
        angle_spacing = 5
    score_path = os.path.join(output_dir,"score.pkl")
    vesper_script = os.path.join(os.getcwd(),"VESPER_CUDA")
    vesper_script = os.path.join(vesper_script,"main.py")
    if os.path.exists(os.path.join(output_dir,"PDB")):
        functToDeleteItems(os.path.join(output_dir,"PDB"))
    command_line="python3 %s orig -a %s -t %f -b %s -T %f -g %f -s %f " \
                                 "-A %f -N %d -M %s -gpu 0 -o %s -pdbin %s -ca %s -ldp %s -c %d >%s"\
                                 %(vesper_script,input_map_path,params['vesper']['ldp_cutoff'],gen_reso_map_path,
                                 params['vesper']['simu_cutoff'],kernel_size,
                                   voxel_spacing,angle_spacing,
                                   params['vesper']['num_models'],params['vesper']['rank_mode'],output_dir,
                                   input_pdb_path,backbone_pdb,ldp_pdb_path,params['vesper']['thread'],output_pdb_path
                                   )
    os.system(command_line)
    new_score_dict={}
    pdb_dir=os.path.join(output_dir,"PDB")
    read_score(new_score_dict,pdb_dir,output_pdb_path)#score sorted inside the func
    write_pickle(new_score_dict,score_path)
    return new_score_dict
