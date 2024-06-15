
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
        map_name="input_diffmodeler" #to avoid server same name bugs
        mkdir(save_path)
    try:
        print("pre-compile VESPER to accelerate!")
        running_dir = os.path.dirname(os.path.abspath(__file__))
        os.system(f"cd {running_dir}; python -O -m compileall VESPER_CUDA")
    except:
        print("pre-compile VESPER failed! No impact to main scripts!")
    save_path = os.path.abspath(save_path)
    from data_processing.Unify_Map import Unify_Map
    cur_map_path = Unify_Map(cur_map_path,os.path.join(save_path,map_name+"_unified.mrc"))
    from data_processing.Resize_Map import Resize_Map
    cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
    if params['contour']<0:
        #change contour level to 0 and increase all the density
        from ops.map_utils import increase_map_density
        cur_map_path = increase_map_density(cur_map_path,os.path.join(save_path,map_name+"_increase.mrc"),params['contour'])
        params['contour']=0
    from modeling.map_utils import segment_map
    new_map_path = os.path.join(save_path,map_name+"_segment.mrc")
    segment_map(cur_map_path,new_map_path,contour=0)
    return save_path,new_map_path

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
from ops.io_utils import delete_dir
def construct_single_chain_candidate(params,save_path):
    #first build a dict from the input text configure file
    single_chain_pdb_input = os.path.abspath(params['P'])
    single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
    
    #if os.path.exists(single_chain_pdb_dir):
        #shutil.rmtree(single_chain_pdb_dir)
    #delete_dir(single_chain_pdb_dir)
    if not os.path.isdir(single_chain_pdb_input):
        os.makedirs(single_chain_pdb_dir,exist_ok=True)
        from ops.os_operation import extract_compressed_file
        single_chain_pdb_dir=extract_compressed_file(single_chain_pdb_input,single_chain_pdb_dir)
    else:
        if os.path.exists(single_chain_pdb_dir):
            shutil.rmtree(single_chain_pdb_dir)
        from ops.os_operation import copy_directory
        single_chain_pdb_dir = copy_directory(single_chain_pdb_input,single_chain_pdb_dir)
    from ops.io_utils import read_structure_txt
    fitting_dict = read_structure_txt(single_chain_pdb_dir,os.path.abspath(params['M']))
    return fitting_dict



if __name__ == "__main__":
    params = argparser()
    save_path,cur_map_path = set_up_envrionment(params)
    running_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(running_dir)
    if not os.path.isabs(params['model']['path']):
        params['model']['path'] =os.path.join(running_dir,params['model']['path'])
    if not os.path.isabs(params['db_exp_path']):
        params['db_exp_path'] = os.path.join(running_dir,params['db_exp_path'])
    if not os.path.isabs(params['db_path']):
        params['db_path'] = os.path.join(running_dir,params['db_path'])
    if params['mode']==0:
        #first build a dict from the input text configure file
        fitting_dict = construct_single_chain_candidate(params,save_path)
    elif params['mode']==1:
        from ops.fasta2pool import fasta2pool
        fitting_dict = fasta2pool(params,save_path)
    elif params['mode']==2:
        from ops.fasta_searchdb import fasta_searchdb
        fitting_dict = fasta_searchdb(params,save_path)
        if params['seq_search']:
            print("sequence search finished!")
            exit()
    elif params['mode']==3:
        #support fasta+template mixed mode. Will first check if template is provided, will use the provided one first. Then for remained results, do db search
        fitting_dict = construct_single_chain_candidate(params,save_path)
        params['fasta_path'] = os.path.abspath(params['fasta_path'])
        from ops.fasta_utils import read_fasta
        chain_dict = read_fasta(params['fasta_path'])
        refined_fasta_path = os.path.join(save_path,"refined_input.fasta")
        from ops.fasta_utils import refine_fasta_input
        fitting_dict = refine_fasta_input(chain_dict,fitting_dict,refined_fasta_path)
        print("updated fitting_dict:",fitting_dict)
        params['P']=refined_fasta_path
        from ops.fasta_searchdb import fasta_searchdb
        additional_fitting_dict = fasta_searchdb(params,save_path)
        #merge two fitting dict
        for key in additional_fitting_dict:
            fitting_dict[key]=additional_fitting_dict[key]
        print("final fitting dict:",fitting_dict)
    else:
        print("mode %d is not supported!"%params['mode'])
        exit()

    if len(fitting_dict)==0:
        print("Empty Template candiate, DiffModeler can not run!!!")
        exit()
    #clean fitting dict to avoid strange pdb cause entire program fail
    final_template_dir = os.path.join(save_path,"final_template_input")
    from ops.pdb_utils import clean_pdb_template
    final_fitting_dict=clean_pdb_template(fitting_dict,final_template_dir)

    if params['domain']:
        from ops.domain_utils import prepare_domain_input
        domain_template_dir = os.path.join(save_path,"domain_template_input")
        final_fitting_dict = prepare_domain_input(final_fitting_dict,
                                            domain_template_dir,
                                            num_cpu=params['search_thread'])
        print("Domain split finished!",final_fitting_dict)
    #diffusion inference
    diff_trace_map = diffusion_trace_map(save_path,cur_map_path,params)

    #VESPER singl-chain fitting process
    fitting_dir = os.path.join(save_path,"structure_modeling")
    from modeling.fit_structure_chain import fit_structure_chain
    fit_structure_chain(diff_trace_map,final_fitting_dict,fitting_dir,params)

    #VESPER assembling
    modeling_dir = os.path.join(save_path,"structure_assembling")
    from modeling.assemble_structure import assemble_structure
    source_cif = assemble_structure(diff_trace_map,final_fitting_dict,fitting_dir,modeling_dir,params)
    output_cif = os.path.join(save_path,"DiffModeler.cif")
    shutil.copy(source_cif,output_cif)
    print(f"Please check DiffModeler's output structure in {output_cif}")




