import numpy as np

def permute_ns_coord_to_pdb(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[2]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[0]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[1]
        out_z = input_coord[2]
    else:
        exit()
    return [out_x, out_y, out_z]

def permute_pdb_coord_to_map(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[2]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[0]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[1]
        out_z = input_coord[2]
    else:
        exit()
    return out_x, out_y, out_z


def permute_map_coord_to_pdb(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[2]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[0]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[1]
        out_z = input_coord[2]
    else:
        exit()
    return out_x, out_y, out_z

def Filter_Map_Location(map_data,  map_info_list,label):
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    #1st filter out all the coordinates in this cluster
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    data= np.array(map_data)
    pure_coord_list = np.argwhere(data==label)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    for k in range(len(pure_coord_list)):
        location = pure_coord_list[k]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
    return All_location

def Convert_Location_between_Maps(location1,origin1,map_info1,origin2,map_info2):
    """

    :param location1: location in the first map
    :param origin1:
    :param map_info1:
    :param origin2:
    :param map_info2:
    :return:
    """
    mapc, mapr,maps = map_info1
    mapc2, mapr2,maps2 = map_info2
    actual_location = permute_map_coord_to_pdb(location1, mapc, mapr,maps)
    actual_location = list(actual_location)
    for k in range(3):
        actual_location[k] = actual_location[k]+origin1[k]-origin2[k]
    new_x, new_y, new_z = permute_pdb_coord_to_map(actual_location,mapc2,mapr2,maps2)
    return [new_x,new_y,new_z]




