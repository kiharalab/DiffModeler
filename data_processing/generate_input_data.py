import mrcfile
import numpy as np
import os
from ops.map_utils import find_top_density,permute_ns_coord_to_pdb
from progress.bar import Bar
from ops.os_operation import mkdir
def gen_input_box(map_data,box_size,stride,contour,train_save_path):
    scan_x, scan_y, scan_z = map_data.shape
    count_voxel = 0
    count_iter=0
    Coord_Voxel = []

    bar = Bar('Preparing Input: ', max=int(np.ceil(scan_x/stride)*np.ceil(scan_y/stride)*np.ceil(scan_z/stride)))


    for x in range(0, scan_x, stride):
        x_end = min(x + box_size, scan_x)
        for y in range(0, scan_y, stride):
            y_end = min(y + box_size, scan_y)
            for z in range(0, scan_z, stride):
                count_iter+=1
                bar.next()
                #print("1st stage: %.4f percent scanning finished"%(count_iter*100/(scan_x*scan_y*scan_z/(stride**3))),"location %d %d %d"%(x,y,z))
                z_end = min(z + box_size, scan_z)
                if x_end < scan_x:
                    x_start = x
                else:
                    x_start = x_end - box_size

                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y
                else:
                    y_start = y_end - box_size

                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z
                else:
                    z_start = z_end - box_size

                    if z_start<0:
                        z_start=0
                #already normalized
                segment_map_voxel = np.zeros([box_size,box_size,box_size])
                segment_map_voxel[:x_end-x_start,:y_end-y_start,:z_end-z_start]=map_data[x_start:x_end, y_start:y_end, z_start:z_end]
                if contour<=0:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel>0))
                    meaningful_density_ratio = meaningful_density_count/float(box_size**3)
                    if meaningful_density_ratio<=0.001:
                        #print("no meaningful density ratio %f in current scanned box, skip it!"%meaningful_density_ratio)
                        continue
                else:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel > contour))
                    meaningful_density_ratio = meaningful_density_count / float(box_size ** 3)
                    if meaningful_density_ratio <= 0.001:
                       # print("no meaningful density ratio in current scanned box, skip it!")
                        continue
                cur_path = os.path.join(train_save_path,"input_"+str(count_voxel)+".npy")
                np.save(cur_path,segment_map_voxel)
                Coord_Voxel.append([x_start,y_start,z_start])
                count_voxel+=1
    bar.finish()
    Coord_Voxel = np.array(Coord_Voxel)
    coord_path = os.path.join(train_save_path,"Coord.npy")
    np.save(coord_path,Coord_Voxel)
    print("In total we prepared %d boxes as input"%(len(Coord_Voxel)))
    return Coord_Voxel

def generate_infer_data(input_map_path,save_input_dir,contour,params):
    box_size= params["model"]["diffusion"]['box_size']
    stride = params["model"]["diffusion"]['stride']

    mkdir(save_input_dir)
    with mrcfile.open(input_map_path, permissive=True) as map_mrc:
         #normalize data
        map_data = np.array(map_mrc.data)
        # get the value serve as 1 in normalization
        map_data[map_data < 0] = 0
        print("map density range: %f %f"%(0,np.max(map_data)))
        percentile_98 = find_top_density(map_data,0.98)

        print("map hist log percentage 98: ",percentile_98)
        map_data[map_data > percentile_98] = percentile_98
        min_value = np.min(map_data)
        max_value = np.max(map_data)
        map_data = (map_data-min_value)/(max_value-min_value)
        nxstart, nystart, nzstart = map_mrc.header.nxstart, \
                                    map_mrc.header.nystart, \
                                    map_mrc.header.nzstart
        orig = map_mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(", "")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        nstart = [nxstart, nystart, nzstart]
        mapc = map_mrc.header.mapc
        mapr = map_mrc.header.mapr
        maps = map_mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nstart = permute_ns_coord_to_pdb(nstart, mapc, mapr, maps)
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]) + float(nstart[k]))

        print("Origin:", new_origin)

        #adjust the contour level by the maximum value
        print("given contour %f"%contour)
        contour = contour/percentile_98
        print("revised contour %f"%contour)
        coord_path = os.path.join(save_input_dir, "Coord.npy")
        if os.path.exists(coord_path):
            Coord_Voxel = np.load(coord_path)
        else:
            Coord_Voxel = gen_input_box(map_data,box_size,stride,contour,save_input_dir)
        return Coord_Voxel,map_data,contour
