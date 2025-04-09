import mrcfile
import numpy as np
import torch
import torch.nn.functional as F
import argparse
import os

def my_reform_1a(input_mrc, output_mrc, use_gpu=False):
    with torch.no_grad() and torch.cuda.amp.autocast(enabled=use_gpu):
        with mrcfile.open(input_mrc, permissive=True) as orig_map:
            orig_voxel_size = np.array([orig_map.voxel_size.x, orig_map.voxel_size.y, orig_map.voxel_size.z])
            if orig_voxel_size[0] == 1 and orig_voxel_size[1] == 1 and orig_voxel_size[2] == 1:
                if os.path.exists(output_mrc):
                    os.remove(output_mrc)
                os.symlink(input_mrc, output_mrc)
                print(f"No resizing needed for {input_mrc}. Copied to {output_mrc}.")
                return
            orig_data = torch.from_numpy(orig_map.data.astype(np.float32).copy()).unsqueeze(0).unsqueeze(0)
            orig_data = orig_data.cuda() if use_gpu else orig_data

            print(f"[{input_mrc}] Original Shape (ZYX): {orig_data.shape[2:]}")
            print(
                f"[{input_mrc}] Original Voxel Size (ZYX): {orig_map.voxel_size.z:.4f}, {orig_map.voxel_size.y:.4f}, {orig_map.voxel_size.x:.4f}"
            )

            new_grid_size = np.array(orig_data.shape[2:]) * np.array([orig_map.voxel_size.z, orig_map.voxel_size.y, orig_map.voxel_size.x])
            print(f"[{input_mrc}] New Grid Size (ZYX): {new_grid_size[0]}, {new_grid_size[1]}, {new_grid_size[2]}")
            new_grid_size = np.floor(new_grid_size).astype(np.int32)  # ZYX
            print(f"[{input_mrc}] New Grid Size (int) (ZYX): {new_grid_size[0]}, {new_grid_size[1]}, {new_grid_size[2]}")

            kwargs = {"indexing": "ij"} if (int(torch.__version__.split(".")[0]) >= 2 or int(torch.__version__.split(".")[1]) >= 10) else {}

            z = (
                torch.arange(0, new_grid_size[0], device="cuda" if use_gpu else "cpu") / orig_voxel_size[2] / (orig_data.shape[2] - 1) * 2
                - 1
            )
            y = (
                torch.arange(0, new_grid_size[1], device="cuda" if use_gpu else "cpu") / orig_voxel_size[1] / (orig_data.shape[3] - 1) * 2
                - 1
            )
            x = (
                torch.arange(0, new_grid_size[2], device="cuda" if use_gpu else "cpu") / orig_voxel_size[0] / (orig_data.shape[4] - 1) * 2
                - 1
            )

            new_grid = torch.stack(
                torch.meshgrid(x, y, z, **kwargs),
                dim=-1,
            )

            new_grid = new_grid.unsqueeze(0)
            new_data = F.grid_sample(orig_data, new_grid, mode="bilinear", align_corners=True).cpu().numpy()[0, 0]

            new_voxel_size = np.array((1.0, 1.0, 1.0))

            new_data = new_data.transpose((2, 1, 0))

            print(f"[{input_mrc}] New Shape (ZYX): {new_data.shape}")
            print(f"[{input_mrc}] New Voxel Size (ZYX): {new_voxel_size[2]:.4f}, {new_voxel_size[1]:.4f}, {new_voxel_size[0]:.4f}")

            with mrcfile.new(output_mrc, data=new_data.astype(np.float32), overwrite=True) as mrc:
                vox_sizes = mrc.voxel_size
                vox_sizes.flags.writeable = True
                vox_sizes.x = new_voxel_size[0]
                vox_sizes.y = new_voxel_size[1]
                vox_sizes.z = new_voxel_size[2]
                mrc.voxel_size = vox_sizes
                mrc.update_header_from_data()
                mrc.header.nxstart = 0
                mrc.header.nystart = 0
                mrc.header.nzstart = 0
                mrc.header.origin.x = orig_map.header.origin.x + orig_map.header.nxstart * orig_voxel_size[0]
                mrc.header.origin.y = orig_map.header.origin.y + orig_map.header.nystart * orig_voxel_size[1]
                mrc.header.origin.z = orig_map.header.origin.z + orig_map.header.nzstart * orig_voxel_size[2]
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
