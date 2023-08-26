
import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
import shutil
def init_save_path(origin_map_path):
    save_path = os.path.join(os.getcwd(), 'Predict_Result')
    mkdir(save_path)
    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    map_name = map_name.replace(".map", "")
    map_name = map_name.replace("(","").replace(")","")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    return save_path,map_name

def set_up_envrionment(params):
    if params['resolution']>20:
        print("maps with %.2f resolution is not supported! We only support maps with resolution 0-20A!"%params['resolution'])
        exit()
    gpu_id = params['gpu']
    if gpu_id is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
    cur_map_path = os.path.abspath(params['F'])
    if cur_map_path.endswith(".gz"):
        from ops.os_operation import unzip_gz
        cur_map_path = unzip_gz(cur_map_path)

    if params['output'] is None:
        save_path,map_name = init_save_path(cur_map_path)
    else:
        save_path=params['output']
        map_name="input"
        mkdir(save_path)
    from data_processing.Unify_Map import Unify_Map
    cur_map_path = Unify_Map(cur_map_path,os.path.join(save_path,map_name+"_unified.mrc"))
    from data_processing.Resize_Map import Resize_Map
    cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
    return save_path,cur_map_path

def diffusion_trace_map(save_path,cur_map_path,params):
    if params['resolution']>=2:
        from predict.infer_diffusion import infer_diffem
        diffusion_dir = os.path.join(save_path,"infer_diffusion")
        diff_trace_map=infer_diffem(cur_map_path,diffusion_dir,params)
    else:
        print("skip diffusion with very high resolution map %f"%params['resolution'])
        diff_trace_map = cur_map_path
    print("Diffusion process finished! Traced map saved here %s"%diff_trace_map)
    # segment this difftrace map to save time
    from modeling.map_utils import segment_map
    diff_new_trace_map = os.path.join(save_path,"diffusion.mrc")
    segment_map(diff_trace_map,diff_new_trace_map,contour=0)
    return diff_new_trace_map

def construct_single_chain_candidate(params,save_path):
    #first build a dict from the input text configure file
    single_chain_pdb_input = os.path.abspath(params['P'])

    if not os.path.isdir(single_chain_pdb_input):
        single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
        from ops.os_operation import extract_compressed_file
        extract_compressed_file(single_chain_pdb_input,single_chain_pdb_dir)

    else:
        from ops.os_operation import copy_directory
        single_chain_pdb_dir = copy_directory(single_chain_pdb_input,os.path.join(save_path,"single_chain_pdb"))
    from ops.io_utils import read_structure_txt
    fitting_dict = read_structure_txt(single_chain_pdb_dir,os.path.abspath(params['M']))
    return fitting_dict



if __name__ == "__main__":
    params = argparser()

    if params['mode']==0:

        #diffusion inference
        save_path,cur_map_path = set_up_envrionment(params)
        diff_trace_map = diffusion_trace_map(save_path,cur_map_path,params)
        #first build a dict from the input text configure file
        fitting_dict = construct_single_chain_candidate(params,save_path)

        #VESPER singl-chain fitting process
        fitting_dir = os.path.join(save_path,"structure_modeling")
        from modeling.fit_structure_chain import fit_structure_chain
        fit_structure_chain(diff_trace_map,fitting_dict,fitting_dir,params)

        #VESPER assembling
        modeling_dir = os.path.join(save_path,"structure_assembling")
        from modeling.assemble_structure import assemble_structure
        source_cif = assemble_structure(diff_trace_map,fitting_dict,fitting_dir,modeling_dir,params)
        output_cif = os.path.join(save_path,"DiffModeler.cif")
        shutil.copy(source_cif,output_cif)
    elif params['mode']==1:
        #diffusion inference
        save_path,cur_map_path = set_up_envrionment(params)
        diff_trace_map = diffusion_trace_map(save_path,cur_map_path,params)
        from ops.fasta2pool import fasta2pool
        fitting_dict = fasta2pool(params,save_path)

        #VESPER singl-chain fitting process
        fitting_dir = os.path.join(save_path,"structure_modeling")
        from modeling.fit_structure_chain import fit_structure_chain
        fit_structure_chain(diff_trace_map,fitting_dict,fitting_dir,params)

        #VESPER assembling
        modeling_dir = os.path.join(save_path,"structure_assembling")
        from modeling.assemble_structure import assemble_structure
        source_cif = assemble_structure(diff_trace_map,fitting_dict,fitting_dir,modeling_dir,params)
        output_cif = os.path.join(save_path,"DiffModeler.cif")
        shutil.copy(source_cif,output_cif)
