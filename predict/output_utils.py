import numpy as np
import os
from ops.map_utils import save_label_map

def merge_diffusion_map(input_map_path,input_map_data, cur_box_path,output_dir,Coord_Box,
                        output_box_shape,params,save_key,contour):
    output_box = np.zeros(output_box_shape)-1#-1000
    num_channels,scan_x,scan_y,scan_z = output_box_shape

    box_size = params['model']['diffusion']['box_size']

    for k in range(len(Coord_Box)):

        x_start, y_start, z_start = Coord_Box[k]
        x_end = min(x_start + box_size, scan_x)
        y_end = min(y_start+box_size,scan_y)
        z_end = min(z_start+box_size,scan_z)

        if x_end < scan_x:
            x_start = x_start
        else:
            x_start = x_end - box_size
            if x_start<0:
                x_start=0
        if y_end < scan_y:
            y_start = y_start
        else:
            y_start = y_end - box_size
            if y_start<0:
                y_start=0
        if z_end < scan_z:
            z_start = z_start
        else:
            z_start = z_end - box_size
            if z_start<0:
                z_start=0
        cur_pred_path = os.path.join(cur_box_path,"Predict_%d.npy"%k)
        cur_prediction = np.load(cur_pred_path)
        output_box[:,x_start:x_end,y_start:y_end,z_start:z_end]=\
            np.maximum(output_box[:,x_start:x_end,y_start:y_end,z_start:z_end],
                       cur_prediction[:,:x_end-x_start,:y_end-y_start,:z_end-z_start])
    for k in range(len(output_box)):
        cur_save_path=os.path.join(output_dir,save_key+"_c%d.mrc"%k)
        save_label_map(input_map_path,cur_save_path,output_box[k])
        cur_save_path=os.path.join(output_dir,save_key+"_c%d_mask.mrc"%k)
        tmp_box = output_box[k]
        mask = input_map_data<=contour
        tmp_box[mask]=0
        save_label_map(input_map_path,cur_save_path,tmp_box)





