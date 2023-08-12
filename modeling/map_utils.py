import os
from ops.pdb_utils import read_proteinpdb_data
from ops.io_utils import write_pickle
from ops.map_utils import process_map_data,save_label_map
import numpy as np
from ops.map_coord_utils import permute_ns_coord_to_pdb,permute_pdb_coord_to_map

def assign_label_special(map_data, mapc, mapr, maps, origin ,nxstart,nystart,nzstart,cif_info_dict):
    #here nstart is useless, must removed because pdb2vol's generated map
    """
    :param map_data: map array
    :param mapc: map column indicator
    :param mapr: map row indicator
    :param maps: map section indicator
    :param origin: starting point of map, x, y,z
    :param cif_info_dict:
    :return:
    """
    output_atom_label = np.zeros(map_data.shape)#-1 denotes background
    distance_array = np.ones(map_data.shape) * 100 #record the closest distance
    count_check = 0
    for tmp_id in cif_info_dict:
        tmp_dict = cif_info_dict[tmp_id]
        current_coord = tmp_dict['coord']
        current_atom = tmp_dict['atom']
        current_radius =tmp_dict['radius']
        current_nuc = tmp_dict['nuc']
        nstart = [nxstart,nystart,nzstart]
        nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
        output_coord = []
        for k in range(3):
            output_coord.append(current_coord[k]-origin[k])
        new_x, new_y, new_z = permute_pdb_coord_to_map(output_coord,mapc,mapr,maps)
        #print("input x %.4f, y %.4f, z %.4f"%(output_coord[0],output_coord[1],output_coord[2]))
        #print("output x %.4f, y %.4f, z %.4f"%(new_x, new_y, new_z))
        check_radius = int(np.ceil(current_radius))
        int_x, int_y,int_z = int(new_x),int(new_y),int(new_z)
        square_limit = current_radius ** 2
        for ii in range(int_x-check_radius,int_x+check_radius+1):
            if ii >= distance_array.shape[0]:
                continue
            for jj in range(int_y-check_radius,int_y+check_radius+1):
                if jj >=distance_array.shape[1]:
                    continue
                for kk in range(int_z-check_radius,int_z+check_radius+1):
                    if kk>=distance_array.shape[2]:
                        continue
                    cur_distance = (ii-new_x)**2+(jj-new_y)**2+(kk-new_z)**2
                    if cur_distance<= square_limit:
                        record_distance = distance_array[ii,jj,kk]
                        if record_distance>=cur_distance:
                            #update record
                            distance_array[ii,jj,kk] = cur_distance
                            output_atom_label[ii,jj,kk] = current_atom
        count_check+=1
        if count_check%100==0:
            print("Calcu finished %d/%d"%(count_check,len(cif_info_dict)))
    return distance_array,output_atom_label

def mask_map_by_pdb(input_map_path,output_map_path,final_pdb_output,cutoff=2):
    output_dir = os.path.split(output_map_path)[0]
    pdb_info_dict = read_proteinpdb_data(final_pdb_output,run_type=1,atom_cutoff=cutoff)
    pdb_info_path = os.path.join(output_dir,"maskprotein_info_%d.pkl"%cutoff)
    write_pickle(pdb_info_dict,pdb_info_path)
    #pickle.dump(pdb_info_dict, open(pdb_info_path, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    map_data, mapc, mapr, maps, origin,nxstart,nystart,nzstart = process_map_data(input_map_path)
    distance_array, output_atom_label = \
                assign_label_special(map_data, mapc, mapr, maps, origin,nxstart,nystart,nzstart,  pdb_info_dict)
    map_data = np.array(map_data)
    count_meaningful = len(map_data[map_data!=0])
    count_masked = len(output_atom_label[output_atom_label!=0])
    print("mask %d grid points out of %d meaningful grid points. "
          "Ratio:%.2f"%(count_masked,count_meaningful,count_masked/count_meaningful))
    map_data[output_atom_label!=0]=0
    save_label_map(input_map_path,output_map_path,map_data)

