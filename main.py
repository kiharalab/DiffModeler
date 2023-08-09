
import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
def init_save_path(origin_map_path):
    save_path = os.path.join(os.getcwd(), 'Predict_Result')
    mkdir(save_path)
    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    map_name = map_name.replace(".map", "")
    map_name = map_name.replace("(","").replace(")","")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    return save_path,map_name


if __name__ == "__main__":
    params = argparser()

    if params['mode']==0:
        gpu_id = params['gpu']
        if gpu_id is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        cur_map_path = os.path.abspath(params['F'])
        if ".gz"==cur_map_path[-3:]:
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
        #diffusion inference
        from predict.infer_diffusion import infer_diffem
        diffusion_dir = os.path.join(save_path,"infer_diffusion")
        diff_trace_map=infer_diffem(cur_map_path,diffusion_dir,params)
        print("Diffusion process finished! Traced map saved here %s"%diff_trace_map)
        #first build a dict from the input text configure file
        from ops.io_utils import read_structure_txt
        fitting_dict = read_structure_txt(os.path.abspath(params['M']))

        #VESPER singl-chain fitting process
        fitting_dir = os.path.join(save_path,"structure_modeling")
        from modeling.fit_structure_chain import fit_structure_chain
        fit_structure_chain(diff_trace_map,fitting_dict,fitting_dir,params)


