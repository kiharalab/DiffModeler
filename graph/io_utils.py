
import mrcfile
import numpy as np
def save_LDP_map(save_map_path,LDP_array,origin_map_path):

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
        data= mrc.data
        prediction = np.zeros(data.shape)
        for k in range(len(LDP_array)):
            x,y,z,density = LDP_array[k]
            x, y, z = int(x), int(y),int(z)
            prediction[x,y,z]=density

        data_new = np.float32(prediction)
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


def append_cif_info(cur_entry_id,cur_frag_path,fragment_all_path):
    with open(fragment_all_path,'a+') as wfile:
        with open(cur_frag_path,'r') as rfile:
            wfile.write("#\n")
            wfile.write("data_"+str(cur_entry_id)+"\n")
            wfile.write("_entry.id "+str(cur_entry_id)+"\n")
            wfile.write("#\n")
            for j,line in enumerate(rfile):
                if j>=2:
                    wfile.write(line)
