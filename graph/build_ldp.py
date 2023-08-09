from structure.MRC import MRC
import numpy as np
import mrcfile
from ops.map_utils import process_map_data
from structure.Points import Points
from graph.io_utils import save_LDP_map
from graph.visualize_utils import Show_Graph_Connect,Show_Bfactor_cif
from ops.os_operation import mkdir
import os
from graph.LDP_ops import Extract_LDP_coord,calculate_merge_point_density

def build_ldp(input_map_path,output_pdb_path,params):
    with mrcfile.open(input_map_path,permissive=True) as mrc:
        chain_prob = mrc.data
    chain_prob = np.array(chain_prob)
    input_mrc = MRC(input_map_path, params['g'])
    keyword='chain'
    input_mrc.upsampling_specify_prob(keyword,chain_prob,threshold=params['threshold'],filter_array=None)
    map_data, mapc, mapr, maps, origin, nxstart, nystart, nzstart = process_map_data(input_map_path)
    map_info_list=[mapc, mapr, maps, origin, nxstart, nystart, nzstart]
    construct_LDP(input_mrc,input_mrc.__dict__['%s_dens'%keyword], input_mrc.__dict__['%s_Nact'%keyword],
                                                input_map_path,output_pdb_path,params,
                                                map_info_list,keyword)#both works same with relax_LDP=True/False

def construct_LDP(input_mrc,sugar_density, sugar_Nact,
                  origin_map_path,save_path,
                  params,map_info_list,ext_name):
    mkdir(save_path)
    mean_shift_path = os.path.join(save_path, ext_name+'_mean_shift')
    sugar_point = Points(params, sugar_Nact)
    if not os.path.exists(mean_shift_path + '_cd.txt') \
            or not os.path.exists(mean_shift_path + '_dens.txt'):
        input_mrc.general_mean_shift(sugar_density,sugar_point, mean_shift_path)
    else:
        input_mrc.load_general_mean_shift(sugar_density,sugar_point, mean_shift_path)
    #change the density to common value without dividing the Nori
    #sugar_point.recover_density()
    sugar_point_path = os.path.join(save_path, ext_name+'_point.txt')
    # init_id, x,y,z,density, merged_to_id
    # can use the x,y,z here to assign detailed probability for each LDP points.
    if not os.path.exists(sugar_point_path):
        sugar_point.Merge_point(input_mrc, sugar_point_path)  # You will get a merged point file here.
    else:
        sugar_point.load_merge(input_mrc, sugar_point_path)
    merged_cd_dens = np.loadtxt(sugar_point_path[:-4] + 'onlymerged.txt')
    LDP_save_path = os.path.join(save_path,ext_name+"_LDP.mrc")
    save_LDP_map(LDP_save_path, merged_cd_dens, origin_map_path)
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    All_location = Extract_LDP_coord(merged_cd_dens,mapc, mapr, maps, origin, nxstart, nystart, nzstart)
    graph_path = os.path.join(save_path, ext_name+"_LDP.pdb")
    Show_Graph_Connect(All_location, [], graph_path)
    #Get each LDP's sum probability values from its neighbors
    sugar_point = calculate_merge_point_density(sugar_point,sugar_density)
    #plot in b factor
    ldp_prob_path = os.path.join(save_path, ext_name+"_LDPdens.cif")
    Show_Bfactor_cif(ext_name+"_dens",All_location,ldp_prob_path,sugar_point.merged_cd_dens[:,3])
    #turns out density showed very good results its correlation with real phosphate positions
    return sugar_point
