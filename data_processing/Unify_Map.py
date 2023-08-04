from pathlib import Path

import mrcfile
import numpy as np


def Unify_Map(input_map_path, new_map_path):
    # Read MRC file
    mrc = mrcfile.open(input_map_path, permissive=True)

    data = mrc.data.copy()
    voxel_size = np.asarray(mrc.voxel_size.tolist(), dtype=np.float32)
    origin = np.array(mrc.header.origin.tolist(), dtype=np.float32)
    nstart = np.asarray([mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart], dtype=np.float32)
    cella = np.array(mrc.header.cella.tolist(), dtype=np.float32)
    mapcrs = np.asarray([mrc.header.mapc, mrc.header.mapr, mrc.header.maps], dtype=int)
    if np.sum(nstart) == 0:
        return input_map_path
    mrc.print_header()
    mrc.close()

    # Reorder
    sort = np.asarray([0, 1, 2], dtype=np.int64)
    for i in range(3):
        sort[mapcrs[i] - 1] = i
    nstart = np.asarray([nstart[i] for i in sort])
    data = np.transpose(data, axes=2 - sort[::-1])

    # Move offsets from nstart to origin
    origin = origin + nstart * voxel_size

    # Save the unified map
    mrc_new = mrcfile.new(new_map_path, data=data, overwrite=True)
    (mrc_new.header.origin.x, mrc_new.header.origin.y, mrc_new.header.origin.z) = origin
    (mrc_new.header.cella.x, mrc_new.header.cella.y, mrc_new.header.cella.z) = cella
    mrc_new.update_header_stats()
    mrc_new.print_header()
    mrc_new.close()

    # Validate the new MRC file
    mrcfile.validate(new_map_path)
    return new_map_path


if __name__ == "__main__":
    input_map_path = Path('../example/21051.mrc')
    new_map_path = Path('../example/21051_unified.mrc')
    new_map_path = Unify_Map(input_map_path, new_map_path)
    print(f"New map path is {new_map_path}")
