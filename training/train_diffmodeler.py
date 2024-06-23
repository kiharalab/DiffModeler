
from training.train_utils import AverageMeter,ProgressMeter
import time

def train_diffmodeler(train_loader,ddim_runner,model,
                                       epoch,train_log_path,params):
    model.train()
    avg_meters = {'data_time':AverageMeter('data_time'),
                  'train_time':AverageMeter('train_time'),
                  'iou':AverageMeter('iou'),
                  'loss': AverageMeter('loss'),}
    
    progress = ProgressMeter(
        len(train_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         avg_meters['iou'],
         avg_meters['loss']]
        ,prefix="Epoch: [{}]".format(epoch))
    end_time=time.time()
    batch_size = params['train']['batch_size']

    for batch_idx,data in enumerate(train_loader):
        ddim_runner.feed_data(data)
        avg_meters['data_time'].update(time.time()-end_time,batch_size)
        loss_dict = ddim_runner.optimize_parameters()
        avg_meters['loss'].update(loss_dict['loss'], batch_size)
        avg_meters['iou'].update(loss_dict['iou'], batch_size)
        avg_meters['train_time'].update(time.time()-end_time,batch_size )
        end_time = time.time()
        progress.display(batch_idx)
        progress.write_record(batch_idx,train_log_path)
    return avg_meters['loss'].avg
