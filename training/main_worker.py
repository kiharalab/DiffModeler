import os
import datetime
import time 
import shutil
import torch.backends.cudnn as cudnn
import random
import torch
from torch.optim import lr_scheduler
from data_processing.Train_Dataset import Train_Dataset
from training.train_utils import save_checkpoint
from training.train_diffmodeler import train_diffmodeler
from training.valid_diffmodeler import valid_diffmodeler
def init_save_path(params):
    """
    configure the save path for the training
    """
    output_path = os.path.abspath(params['output'])
    today = datetime.date.today()
    formatted_today = today.strftime('%y%m%d')
    now = time.strftime("%H:%M:%S")
    log_path = os.path.join(output_path,formatted_today+now)
    model_path = os.path.join(log_path,'model')
    os.makedirs(model_path,exist_ok=True)
    return log_path,model_path
def read_mapid_info(train_info_path):
    train_list = []
    with open(train_info_path,'r') as file:
        line = file.readline()
        while line:
            line = line.strip("\n")
            split_result = line.split()
            train_list.append(split_result[0])
            line =file.readline()
    return train_list

def configure_dataset(params):
    """
    configure the dataset for the training
    """
    config_info_txt = os.path.abspath(params['info_txt'])
    train_list = read_mapid_info(config_info_txt)
    random.seed(params['train']['rand_seed'])
    random.shuffle(train_list)
    train_list_len = len(train_list)
    train_list_len = int(train_list_len*params['train']['train_ratio'])
    val_list = train_list[train_list_len:]
    train_list = train_list[:train_list_len]
    input_dataset_path = os.path.abspath(params['F'])
    
    train_dataset = Train_Dataset(input_dataset_path,train_list,
                                  params['model']['diffusion']['box_size'],
                                  params['data']['input_channel'],
                                  params['data']['output_channel'],is_train=True)
    val_dataset = Train_Dataset(input_dataset_path,val_list,
                                params['model']['diffusion']['box_size'],
                                params['data']['input_channel'],
                                params['data']['output_channel'],is_train=False)
    
    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=params['train']['batch_size'],
        shuffle=True,
        num_workers=params['train']['num_workers'],
        pin_memory=True,
        drop_last=True)
    val_loader = torch.utils.data.DataLoader(
        val_dataset,
        batch_size=params['train']['batch_size'],
        shuffle=False,
        num_workers=params['train']['num_workers'],
        pin_memory=True,
        drop_last=False)

    return train_loader,val_loader

def main_worker(params):
    log_path,model_path = init_save_path(params)
    #copy the config path to log dir to record the configuration
    config_path = os.path.join(log_path,'config.json')
    shutil.copy(params['config'],config_path)

    from model.DDIM import DDIM
    ddim_runner = DDIM(params)
    model = ddim_runner.netG
    optimizer = ddim_runner.optG
    print("model and optimizer initialized.")
    cudnn.benchmark = True
    scheduler = lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=params['train']['epoch'], eta_min=params['train']["optimizer"]["min_lr"])

    #configure the dataset
    train_loader,val_loader = configure_dataset(params)
    print("dataset configured.")


    start_epoch=0
    best_loss =10000
    if params['resume']['flag']:
        resume_model_path=os.path.abspath(params['resume']['path'])
        state_dict = torch.load(resume_model_path)
        model.load_state_dict(state_dict['state_dict'])
        optimizer.load_state_dict(state_dict['optimizer'])
        scheduler.load_state_dict(state_dict['scheduler'])
        best_loss = state_dict['best_loss']
        start_epoch = state_dict['epoch']

    train_log_path = os.path.join(log_path, "train.log")
    val_log_path = os.path.join(log_path, "val.log")

    for epoch in range(start_epoch,params['train']['epoch']):
        print('Epoch [%d/%d]' % (epoch, params['train']['epoch']))
        train_loss = train_diffmodeler(train_loader,ddim_runner,model,
                                       epoch,train_log_path,params)
        print("*Training loss %.5f:"%train_loss)

        val_loss = valid_diffmodeler(val_loader,ddim_runner,model,
                                     epoch,val_log_path,params) 

        print("*Validation loss %.5f:"%val_loss)
        scheduler.step()
        is_best = best_loss > val_loss
        best_loss = min(best_loss, val_loss)
        save_dict = {
            'epoch': epoch + 1,
            'best_loss': best_loss,
            'val_loss': val_loss,
            'state_dict': model.state_dict(),
            'optimizer': optimizer.state_dict(),
            'scheduler': scheduler.state_dict(),
        }
        save_model_path = os.path.join(model_path, "checkpoint.pth.tar")
        save_checkpoint(save_dict, is_best, save_model_path)
        if epoch % params['train']['save_checkpoint_freq'] == 0 and epoch > 0:
            save_model_path = os.path.join(model_path, "model%d.pth.tar" % epoch)
            save_checkpoint(save_dict, False, save_model_path)


