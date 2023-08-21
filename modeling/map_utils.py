import os
from ops.pdb_utils import read_proteinpdb_data
from ops.io_utils import write_pickle
from ops.map_utils import process_map_data,save_label_map
import numpy as np
from ops.map_coord_utils import permute_ns_coord_to_pdb,permute_pdb_coord_to_map
import mrcfile


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

def mask_map_by_pdb_slow(input_map_path,output_map_path,final_pdb_output,cutoff=2,keep_label=False):
    """

    :param input_map_path:
    :param output_map_path:
    :param final_pdb_output:
    :param cutoff:
    :param keep_label: keep labeled region, or keep unlabeled region.
    :return:
    """
    output_dir = os.path.split(output_map_path)[0]
    pdb_info_dict = read_proteinpdb_data(final_pdb_output,run_type=1,atom_cutoff=cutoff)
    pdb_info_path = os.path.join(output_dir,"maskprotein_info_%d.pkl"%cutoff)
    write_pickle(pdb_info_dict,pdb_info_path)
    #pickle.dump(pdb_info_dict, open(pdb_info_path, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    map_data, mapc, mapr, maps, origin,nxstart,nystart,nzstart = process_map_data(input_map_path)
    distance_array, output_atom_label = \
                assign_label_special(map_data, mapc, mapr, maps, origin,nxstart,nystart,nzstart,  pdb_info_dict)
    map_data = np.array(map_data)
    # count_meaningful = len(map_data[map_data!=0])
    # count_masked = len(output_atom_label[output_atom_label!=0])
    # print("mask %d grid points out of %d meaningful grid points. "
    #       "Ratio:%.2f"%(count_masked,count_meaningful,count_masked/count_meaningful))
    if keep_label:
        #unassigned region all set to 0, only focused on the labeled region
        map_data[output_atom_label==0]=0
    else:
        map_data[output_atom_label!=0]=0
    save_label_map(input_map_path,output_map_path,map_data)
from TEMPy import *
from TEMPy.protein.structure_parser import *
from TEMPy.protein.structure_blurrer import *
from TEMPy.protein.scoring_functions import *
from TEMPy.maps.map_parser import *
def mask_map_by_pdb(input_map_path,output_map_path,final_pdb_output,keep_label=False):
    """

    :param input_map_path:
    :param output_map_path:
    :param final_pdb_output:
    :param cutoff:
    :param keep_label: keep labeled region, or keep unlabeled region.
    :return:
    """
    #first use pdb to blue simulated map
    if keep_label:
        resolution =10#includes bigger map as possible
    else:
        resolution = 2
    sb = StructureBlurrer()
    pdb1 = PDBParser.read_PDB_file('PDB1',final_pdb_output)
    map1 = MapParser.readMRC(input_map_path)
    sim_map = sb.gaussian_blur_real_space(prot=pdb1,densMap=map1,resolution=resolution)
    bin_map1 = map1.fullMap
    bin_map2 = sim_map.fullMap
    if keep_label:
        mask = bin_map2>=0.01
    else:
        mask = bin_map2<0.001
    bin_map1 = bin_map1*mask
    map1.fullMap = bin_map1
    map1.write_to_MRC_file(output_map_path)

def segment_map(input_map,output_map,contour=0):
    """
    segment meaningful region of a map
    :param input_map:
    :param output_map:
    :return:
    generate a new small size map
    """
    with mrcfile.open(input_map, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        #check the useful density in the input
        input_data = mrc.data
        useful_index = np.argwhere(input_data>contour)
        min_x = int(np.min(useful_index[:,0]))
        max_x = int(np.max(useful_index[:,0]))
        min_y = int(np.min(useful_index[:,1]))
        max_y = int(np.max(useful_index[:,1]))
        min_z = int(np.min(useful_index[:,2]))
        max_z = int(np.max(useful_index[:,2]))
        mapc = mrc.header.mapc
        mapr = mrc.header.mapr
        maps = mrc.header.maps
        new_data = input_data[min_x:max_x,min_y:max_y,min_z:max_z]
        shift_start = permute_ns_coord_to_pdb([min_x,min_y,min_z],mapc,mapr,maps)
        origin = np.array(mrc.header.origin.tolist(), dtype=np.float32)
        origin = np.array(origin)+np.array(shift_start)
        mrc_new = mrcfile.new(output_map, data=new_data, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nx = int(max_x-min_x)
        mrc_new.header.ny = int(max_y-min_y)
        mrc_new.header.nz = int(max_z-min_z)
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mx = int(max_x-min_x)
        mrc_new.header.my = int(max_y-min_y)
        mrc_new.header.mz = int(max_z-min_z)
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.cella.x = int(max_x-min_x)
        mrc_new.header.cella.y = int(max_y-min_y)
        mrc_new.header.cella.z = int(max_z-min_z)
        #mrc_new.header.origin = origin
        (mrc_new.header.origin.x, mrc_new.header.origin.y, mrc_new.header.origin.z) = origin
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()

