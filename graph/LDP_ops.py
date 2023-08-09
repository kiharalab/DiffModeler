import numpy as np
from ops.map_utils import permute_ns_coord_to_pdb,permute_map_coord_to_pdb


def Extract_LDP_coord(merged_cd_dens, mapc, mapr, maps, origin, nxstart, nystart, nzstart):
    #1st filter out all the coordinates in this cluster
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    for node_id in range(len(merged_cd_dens)):
        node_id = int(node_id)
        location = merged_cd_dens[node_id,:3]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
    return All_location
def Convert_LDPcoord_To_Reallocation(Coordinate_List, map_info_list):
    #1st filter out all the coordinates in this cluster
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    for node_id in range(len(Coordinate_List)):
        node_id = int(node_id)
        location =Coordinate_List[node_id]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
    return All_location

def calculate_merge_point_density(point,prob_array):
    before_merge_data = point.merged_data #[merge_to_id,x,y,z, density,merge_status]# 1st column is the id points to the merged points, which is the member id
    after_merge_data = point.merged_cd_dens
    Output_Prob = np.zeros(len(after_merge_data))
    count_isolate=0
    for k in range(len(before_merge_data)):
        merge_to_id,x,y,z,density,merge_status=before_merge_data[k]
        x,y,z = int(x),int(y),int(z)
        if merge_to_id==k and merge_status==-1:
            count_isolate+=1
            continue
        while merge_status==-1:
            merge_to_id = int(merge_to_id)
            merge_to_id,_,_,_,_,merge_status=before_merge_data[merge_to_id]
        final_id = int(merge_status)
        Output_Prob[final_id]+=prob_array[x,y,z]
    point.merge_prob = Output_Prob
    print("in total %d isolated grid points"%count_isolate)
    Output_Prob = np.zeros(len(after_merge_data))
    for k in range(len(after_merge_data)):
        x,y,z,_ = after_merge_data[k]
        x,y,z = int(x),int(y),int(z)
        Output_Prob[k]=prob_array[x,y,z]
    point.point_prob = Output_Prob
    return point
