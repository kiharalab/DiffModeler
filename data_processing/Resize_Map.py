import mrcfile
import numpy as np
import torch
import torch.nn.functional as F
import argparse

def my_reform_1a(input_mrc, output_mrc, use_gpu=False):

    with torch.no_grad() and torch.cuda.amp.autocast(enabled=use_gpu):

        with mrcfile.open(input_mrc, permissive=True) as orig_map:
            voxel_size = np.array([orig_map.voxel_size.x, orig_map.voxel_size.y, orig_map.voxel_size.z])

            # orig_data = torch.from_numpy(orig_map.data.copy().transpose((2, 1, 0))).unsqueeze(0).unsqueeze(0)
            orig_data = torch.from_numpy(orig_map.data.copy()).unsqueeze(0).unsqueeze(0)

            orig_data = orig_data.cuda() if use_gpu else orig_data

            print("Previous shape: ", orig_data.shape)
            # orig = np.array([orig_map.header.origin.x, orig_map.header.origin.y, orig_map.header.origin.z])

            new_grid_size = np.floor(np.array(orig_data.shape[2:]) * voxel_size).astype(np.int32)

            # for compatibility with torch 1.9 and below
            kwargs = {"indexing": "ij"} if (torch.__version__.split(".")[0] >= "2" or torch.__version__.split(".")[0] >= "10") else {}

            new_grid = torch.stack(
                torch.meshgrid(
                    torch.arange(0, new_grid_size[2], device="cuda" if use_gpu else "cpu") / voxel_size[2] / (orig_data.shape[4] - 1) * 2
                    - 1,
                    torch.arange(0, new_grid_size[1], device="cuda" if use_gpu else "cpu") / voxel_size[1] / (orig_data.shape[3] - 1) * 2
                    - 1,
                    torch.arange(0, new_grid_size[0], device="cuda" if use_gpu else "cpu") / voxel_size[0] / (orig_data.shape[2] - 1) * 2
                    - 1,
                    **kwargs
                ),
                dim=-1,
            )

            new_grid = new_grid.unsqueeze(0)
            new_data = F.grid_sample(orig_data, new_grid, mode="bilinear", align_corners=True).cpu().numpy()[0, 0]

            new_voxel_size = np.array((1.0, 1.0, 1.0))
            print("Real voxel size: ", new_voxel_size)
            print("New shape: ", new_data.shape)

            new_data = new_data.transpose((2, 1, 0))

            with mrcfile.new(output_mrc, data=new_data.astype(np.float32), overwrite=True) as mrc:
                vox_sizes = mrc.voxel_size
                vox_sizes.flags.writeable = True
                vox_sizes.x = new_voxel_size[0]
                vox_sizes.y = new_voxel_size[1]
                vox_sizes.z = new_voxel_size[2]
                mrc.voxel_size = vox_sizes
                mrc.update_header_from_data()
                mrc.header.nxstart = orig_map.header.nxstart * voxel_size[0]
                mrc.header.nystart = orig_map.header.nystart * voxel_size[1]
                mrc.header.nzstart = orig_map.header.nzstart * voxel_size[2]
                mrc.header.origin = orig_map.header.origin
                mrc.header.mapc = orig_map.header.mapc
                mrc.header.mapr = orig_map.header.mapr
                mrc.header.maps = orig_map.header.maps
                mrc.update_header_stats()
                mrc.flush()

def Resize_Map(input_map_path,new_map_path):
    try:
        my_reform_1a(input_map_path, new_map_path, use_gpu=True)
    except:
        print("GPU reform failed, falling back to CPU")
        my_reform_1a(input_map_path, new_map_path, use_gpu=False)
    return new_map_path

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("-i", "--input_map_path", type=str, default=None)
    args.add_argument("-o", "--output_map_path", type=str, default=None)
    args = args.parse_args()

    Resize_Map(args.input_map_path,args.output_map_path)