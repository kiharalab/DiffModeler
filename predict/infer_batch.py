

from predict.Logger import AverageMeter,ProgressMeter
import time
import torch
import torch.nn as nn
import os
from ops.os_operation import mkdir
import numpy as np
def infer_batch(val_loader, ddim_runner,model,
                 val_log_path,save_output_dir, params,epoch=0):
    model.eval()
    avg_meters = {'data_time':AverageMeter('data_time'),
                  'train_time':AverageMeter('train_time'),
                  'loss': AverageMeter('loss'),}

    progress = ProgressMeter(
        len(val_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         avg_meters['loss']]
        ,prefix="Epoch: [{}]".format(epoch))

    end_time=time.time()

    with torch.no_grad():
        for batch_idx,data in enumerate(val_loader):
            ddim_runner.feed_data(data)
            current_ids = data['ID']
            batch_size = len(current_ids)
            avg_meters['data_time'].update(time.time()-end_time,batch_size )
            ddim_runner.test(continous=True)
            predict_backbone = ddim_runner.SR.detach().cpu().numpy() #Kept N results during iterations
            # Write all the box segmented numpy into a file
            for kk in range(len(predict_backbone)//batch_size):
                cur_save_dir = os.path.join(save_output_dir,"sample_"+str(kk))
                mkdir(cur_save_dir)
                cur_predict_backbones = predict_backbone[kk*batch_size:(kk+1)*batch_size]
                for jj in range(len(current_ids)):
                    cur_id = int(current_ids[jj])
                    cur_backbone = cur_predict_backbones[jj]
                    cur_save_path = os.path.join(cur_save_dir,"Predict_%d.npy"%cur_id)
                    np.save(cur_save_path,cur_backbone)

            avg_meters['train_time'].update(time.time()-end_time,batch_size )
            end_time = time.time()
            progress.display(batch_idx)
            progress.write_record(batch_idx,val_log_path)



