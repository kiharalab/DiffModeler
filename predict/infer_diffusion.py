
import os
import shutil

from model.DDIM import DDIM
from ops.os_operation import mkdir
from data_processing.generate_input_data import generate_infer_data
from data_processing.Test_Dataset import Single_Dataset
import torch
from predict.output_utils import merge_diffusion_map
from predict.infer_batch import infer_batch
from ops.io_utils import delete_dir
def infer_diffem(input_map_path,save_dir,params):
    mkdir(save_dir)
    final_diffusion_path = os.path.join(save_dir,"trace_backbone.mrc")
    if os.path.exists(final_diffusion_path):
        return final_diffusion_path
    if not os.path.exists(params['model']['path']):
        print("please configure the diffusion model path in config json file.")
        exit()
    params['model']['path'] = os.path.abspath(params['model']['path'])
    #init model
    ddim_runner = DDIM(params)
    model = ddim_runner.netG
    save_input_dir = os.path.join(save_dir,"Input")

    #configure input box
    Coord_Box,map_data,contour = generate_infer_data(input_map_path,save_input_dir,params['contour'],params)
    test_dataset = Single_Dataset(save_input_dir,"input_")

    test_loader = torch.utils.data.DataLoader(
        test_dataset,
        pin_memory=True,
        batch_size=params['model']['batch_size'],
        num_workers=params['model']['num_workers'],
        drop_last=False)

    print("dataset loading finished with %d boxes!"%len(test_dataset))
    #load model
    resume_model_path=os.path.abspath(params['model']['path'])
    state_dict = torch.load(resume_model_path)
    msg=model.load_state_dict(state_dict['state_dict'])
    print("loading model message: ",msg)

    test_log_path = os.path.join(save_dir,"test.log")
    save_diffusion_path = os.path.join(save_dir,"diffusion_box")
    mkdir(save_diffusion_path)
    infer_batch(test_loader, ddim_runner,model, test_log_path, save_diffusion_path, params)


    #for each stage, collecting the overall changed output
    listfiles = os.listdir(save_diffusion_path)
    listfiles.sort()
    final_save_path = os.path.join(save_dir,"diffusion_map")
    mkdir(final_save_path)
    box_shape = ( params['model']['diffusion']['channels'], map_data.shape[0],map_data.shape[1],map_data.shape[2])
    for item in listfiles:
        cur_box_path = os.path.join(save_diffusion_path,item)
        cur_output_path = os.path.join(final_save_path,item)
        mkdir(cur_output_path)
        merge_diffusion_map(input_map_path,map_data,cur_box_path,cur_output_path,Coord_Box,box_shape,params,item,contour)

    id_sort= [int(x.replace("sample_","")) for x in listfiles]
    max_id = int(max(id_sort))
    final_map_path = os.path.join(final_save_path,"sample_%d"%max_id)
    final_map_path = os.path.join(final_map_path,"sample_%d_c0_mask.mrc"%max_id)
    shutil.copy(final_map_path,final_diffusion_path)

    #clean input and box dir, only keep the final map output to save disk space
    #shutil.rmtree(save_input_dir)
    #shutil.rmtree(save_diffusion_path)
    delete_dir(save_input_dir)
    delete_dir(save_diffusion_path)
    return final_diffusion_path
