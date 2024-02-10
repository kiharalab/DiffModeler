import argparse
import mrcfile
import numpy as np


def permute_map_coord_to_pdb(input_coord, mapc, mapr, maps):
    transformations = {
        (1, 2, 3): (input_coord[2], input_coord[1], input_coord[0]),
        (1, 3, 2): (input_coord[2], input_coord[0], input_coord[1]),
        (2, 1, 3): (input_coord[1], input_coord[2], input_coord[0]),
        (2, 3, 1): (input_coord[0], input_coord[2], input_coord[1]),
        (3, 1, 2): (input_coord[1], input_coord[0], input_coord[2]),
        (3, 2, 1): (input_coord[0], input_coord[1], input_coord[2]),
    }

    if (mapc, mapr, maps) in transformations:
        return transformations[(mapc, mapr, maps)]
    else:
        raise ValueError("Invalid mapc, mapr, maps values")


def segment_map(input_map, output_map, contour=0):
    """
    Segments a 3D map based on a given contour level and saves the resulting map to a new file.

    Args:
        input_map (str): Path to the input map file.
        output_map (str): Path to the output map file.
        contour (float, optional): Contour level to use for segmentation. Defaults to 0.

    Returns:
        str: Path to the output map file.
    """
    with mrcfile.open(input_map, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size["x"])
        prev_voxel_size_y = float(prev_voxel_size["y"])
        prev_voxel_size_z = float(prev_voxel_size["z"])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = (
            mrc.header.nx,
            mrc.header.ny,
            mrc.header.nz,
            mrc.header.nxstart,
            mrc.header.nystart,
            mrc.header.nzstart,
            mrc.header.mx,
            mrc.header.my,
            mrc.header.mz,
        )
        orig = mrc.header.origin
        # check the useful density in the input
        input_data = mrc.data
        useful_index = np.argwhere(input_data > contour)
        min_x = int(np.min(useful_index[:, 0]))
        max_x = int(np.max(useful_index[:, 0]))
        min_y = int(np.min(useful_index[:, 1]))
        max_y = int(np.max(useful_index[:, 1]))
        min_z = int(np.min(useful_index[:, 2]))
        max_z = int(np.max(useful_index[:, 2]))
        new_data = input_data[min_x:max_x, min_y:max_y, min_z:max_z]
        # map to new axes
        shift_start = permute_map_coord_to_pdb([min_x, min_y, min_z], int(mrc.header.mapc), int(mrc.header.mapr),
                                               int(mrc.header.maps))
        origin = np.array(mrc.header.origin.tolist(), dtype=np.float32)
        origin = np.array(origin) + np.array(shift_start)
        mrc_new = mrcfile.new(output_map, data=new_data, overwrite=True)
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

        (mrc_new.header.origin.x, mrc_new.header.origin.y, mrc_new.header.origin.z) = origin
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()
    return output_map


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Segment a 3D map based on a given contour level.")
    parser.add_argument("input_map", type=str, help="Path to the input map file.")
    parser.add_argument("output_map", type=str, help="Path to the output map file.")
    parser.add_argument("-c", "--contour", type=float, default=0, help="Contour level to use for segmentation.")
    args = parser.parse_args()
    segment_map(args.input_map, args.output_map, args.contour)
