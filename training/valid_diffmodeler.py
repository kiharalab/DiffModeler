
from training.train_utils import AverageMeter,ProgressMeter
import time
import torch 
def valid_diffmodeler(valid_loader,ddim_runner,model,
                                       epoch,valid_log_path,params):
    model.eval()
    avg_meters = {'data_time':AverageMeter('data_time'),
                  'train_time':AverageMeter('train_time'),
                  'loss': AverageMeter('loss'),}
    progress = ProgressMeter(
        len(valid_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         avg_meters['loss']]
        ,prefix="Epoch: [{}]".format(epoch))
    end_time=time.time()
    batch_size = params['train']['batch_size']
    with torch.no_grad():
        for batch_idx,data in enumerate(valid_loader):
            ddim_runner.feed_data(data)
            avg_meters['data_time'].update(time.time()-end_time,batch_size )
            loss = ddim_runner.calculate_loss()
            avg_meters['loss'].update(loss['loss'].item(), batch_size)
            avg_meters['train_time'].update(time.time()-end_time,batch_size )
            end_time = time.time()
            progress.display(batch_idx)
            progress.write_record(batch_idx,valid_log_path)
    return avg_meters['loss'].avg