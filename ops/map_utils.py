import mrcfile
import numpy as np
from ops.map_coord_utils import permute_ns_coord_to_pdb
import os
def process_map_data(input_map_path):
    with mrcfile.open(input_map_path, permissive=True) as mrc:
        orig = mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(","")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]))
        print("Origin:", new_origin)
        data = mrc.data
        mapc = mrc.header.mapc
        mapr = mrc.header.mapr
        maps = mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nxstart, nystart, nzstart = mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart
        #mapc to x, mapr to y maps to z
    return data, mapc, mapr, maps,new_origin,nxstart,nystart,nzstart

def save_label_map(input_map_path,atom_label_path,output_atom_label):
    """
    :param input_map_path:
    :param atom_label_path: save path
    :param output_atom_label: new array with same shape as input map data
    :return:
    """
    with mrcfile.open(input_map_path, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        output_atom_label = np.float32(output_atom_label)
        mrc_new = mrcfile.new(atom_label_path, data=output_atom_label, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nx = nx
        mrc_new.header.ny = ny
        mrc_new.header.nz = nz
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mx = mx
        mrc_new.header.my = my
        mrc_new.header.mz = mz
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.origin = orig
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()
def save_threshold_map(input_map_path,atom_label_path,threshold):
    """
    :param input_map_path:
    :param atom_label_path: save path
    :param output_atom_label: new array with same shape as input map data
    :return:
    """
    with mrcfile.open(input_map_path, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        output_atom_label = np.array(mrc.data,dtype=np.float32)
        output_atom_label[output_atom_label<=threshold]=0
        mrc_new = mrcfile.new(atom_label_path, data=output_atom_label, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nx = nx
        mrc_new.header.ny = ny
        mrc_new.header.nz = nz
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mx = mx
        mrc_new.header.my = my
        mrc_new.header.mz = mz
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.origin = orig
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()


def extract_origin(input_map_path):
    with mrcfile.open(input_map_path, permissive=True) as map_mrc:
        orig = map_mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(", "")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        nxstart, nystart, nzstart = map_mrc.header.nxstart, \
                                            map_mrc.header.nystart, \
                                            map_mrc.header.nzstart
        nstart = [nxstart, nystart, nzstart]
        mapc = map_mrc.header.mapc
        mapr = map_mrc.header.mapr
        maps = map_mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nstart = permute_ns_coord_to_pdb(nstart, mapc, mapr, maps)
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k])+float(nstart[k]))

    print("Map Origin:", new_origin)
    return new_origin

def extract_naive_origin(input_map_path):
    with mrcfile.open(input_map_path, permissive=True) as map_mrc:
        orig = map_mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(", "")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]))

    print("Map Origin:", new_origin)
    return new_origin

def Write_Local_Map(actual_location,save_path, origin_map_path,segment_map_voxel):
    with mrcfile.open(origin_map_path, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        #orig = mrc.header.origin
        output_atom_label = np.float32(segment_map_voxel)
        mrc_new = mrcfile.new(save_path, data=output_atom_label, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.origin = actual_location
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()


def find_top_density(map_data,threshold):
    use_density=map_data[map_data>0]

    hist,bin_edges=np.histogram(use_density, bins=200)
    #print(hist)
    log_hist = [np.log(x) if x>0 else 0 for x in hist]
    sum_cutoff=np.sum(log_hist)*threshold
    cumulative=0
    for j in range(len(log_hist)):
        cumulative+=log_hist[j]
        if cumulative>sum_cutoff:
            return bin_edges[j]
    return bin_edges[-1]



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


def SimuMap1A_pdb2vol(input_path,save_path):

    split_lists=os.path.split(save_path)
    pdb_name=split_lists[1]
    pdb_id=pdb_name.replace(".pdb","")
    map_dir = split_lists[0]
    reso= 1#random.uniform(5, 10)
    tmp_file_name=os.path.join(map_dir,pdb_id+'_help.txt')
    with open(tmp_file_name,'w') as wFile:
        wFile.write("1\n1\n1\n-" + str(reso) + "\n1\n1\n1\n")

    gen_map_name=save_path
    #os.chdir(new_map_path)
    #gen_map_name=pdb_id.upper()+'.mrc'
    #tmp_file_name=pdb_id+'_help.txt'
    #new_pdb_path=pdb_name
    command_line='cat '+tmp_file_name+" | "+'pdb2vol ' + input_path + " " +gen_map_name
    os.system(command_line)
    return gen_map_name
from TEMPy import *
from TEMPy.protein.structure_parser import *
from TEMPy.protein.structure_blurrer import *
from TEMPy.protein.scoring_functions import *
from TEMPy.maps.map_parser import *
def SimuMap1A(input_path,save_path):
    #tempy based simulated map generation to save time
    sb = StructureBlurrer()
    pdb1 = PDBParser.read_PDB_file('PDB1',input_path)
    sim_map = sb.gaussian_blur_real_space(prot=pdb1,resolution=1)
    sim_map.write_to_MRC_file(save_path)
def save_dens_map(save_map_path,new_dens,
                              origin_map_path):

    with mrcfile.open(origin_map_path, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        print("Origin:", orig)
        print("Previous voxel size:", prev_voxel_size)
        print("nx, ny, nz", nx, ny, nz)
        print("nxs,nys,nzs", nxs, nys, nzs)
        print("mx,my,mz", mx, my, mz)

        data_new = np.float32(new_dens)
        mrc_new = mrcfile.new(save_map_path, data=data_new, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.origin = orig
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()
        del data_new
def increase_map_density(input_path,output_path,add_contour):
    add_contour=abs(add_contour)
    with mrcfile.open(input_path,permissive=True) as mrc:
        data=mrc.data
        data=np.float32(data)
    data = data+add_contour
    save_dens_map(output_path,data,add_contour)
    return output_path
