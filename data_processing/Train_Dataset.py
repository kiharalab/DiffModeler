import os
import numpy as np
import torch
import torch.utils.data
import random
from collections import defaultdict

class Train_Dataset(torch.utils.data.Dataset):
    def __init__(self, data_path,train_id_list,box_size,
                 input_channel, output_channel,
                 is_train=True):
        """
        :param data_path: training data path
        :param training_id: specify the id that will be used for training
        """
        self.is_train=is_train
        self.input_channel = input_channel #should be a list
        self.output_channel = output_channel #should be a list
        self.input_path=[]
        self.output_path=[]
        self.id_list=[]
        self.map_dict= {}
        self.map_box_list=defaultdict(list)
        prev_index=0
        listfiles = os.listdir(data_path)
        for map_name in listfiles:
            #print("collecting data for %s"%map_name)
            cur_dir = os.path.join(data_path,map_name)
            if not os.path.isdir(cur_dir):
                print("expected dir %s for training is empty!"%map_name)
                continue
            map_name = map_name.replace("_new","")
            if map_name not in train_id_list:
                continue
            cur_input_list = [os.path.join(cur_dir, x) for x in os.listdir(cur_dir)
                                  if "input" in x and ".npy" in x]

            cur_nuc_list = [os.path.join(cur_dir, x) for x in os.listdir(cur_dir)
                                if "output" in x and ".npy" in x]
            cur_input_list.sort()
            cur_nuc_list.sort()
            cur_id_list = [map_name+"_"+str(os.path.split(x)[1].replace("output_","").replace(".npy","")) for x in cur_nuc_list]
            assert  len(cur_input_list)==len(cur_nuc_list)

            self.input_path += cur_input_list
            self.output_path += cur_nuc_list
            self.id_list += cur_id_list
            #this is to map the index in the dataset to the map_name and the index in the map_name
            for kk in range(prev_index,len(self.input_path)):
                self.map_dict[kk]=map_name
                self.map_box_list[map_name].append(kk)
            print("accumulating %d training cubes"%len(self.input_path))
            prev_index = len(self.input_path)
        print("in total we have %d training cubes"%len(self.input_path))
    def __len__(self):
        return len(self.input_path)

    def __getitem__(self, idx):
        """
        :param idx:
        :return:
        use index to map its ids.
        """
        inputfile=self.input_path[idx]
        atomfile = self.output_path[idx]
        input = np.load(inputfile)
        atom_label = np.load(atomfile)
        rot = random.randint(0, 3)
        axis = random.randint(0, 2)
        if self.is_train:
            if len(input.shape)==3:
                input = rot90(input, rot, axis)
            else:
                for k in range(len(input)):
                    input[k] = rot90(input[k], rot, axis)
            for k in range(len(atom_label)):
                atom_label[k] = rot90(atom_label[k], rot, axis)

        #change aimset to 0, 1 filling by adding channels
        if len(input.shape)==3:
            input = input[np.newaxis, :]
        input = torch.from_numpy(np.array(input, np.float32, copy=True))
        atom_final_input=[]
        for k in range(len(self.input_channel)):
            atom_final_input.append(input[int(self.input_channel[k]):int(self.input_channel[k])+1])
        atom_final_input = torch.cat(atom_final_input,dim=0)


        atom_label[atom_label<0]=0
        if len(atom_label.shape)==3:
            atom_label = atom_label[np.newaxis,:]
        atom_label = torch.from_numpy(np.array(atom_label, np.float32, copy=True))

        atom_final_label = []
        for k in range(len(self.output_channel)):
            atom_final_label.append(atom_label[int(self.output_channel[k]):int(self.output_channel[k]+1)])
        atom_final_label = torch.cat(atom_final_label,dim=0)

        return {'density': atom_final_input, 'backbone': atom_final_label,  'index': idx}
    



def basic_rot_ax(m, ax=0):
    """
    :param m: 3d matrix
    :return: rotate the cube around axis ax, perpendicular to the face [[0,1],[2,3]]
    """
    ax %= 3
    if ax == 0:
        return np.rot90(m[:, ::-1, :].swapaxes(0, 1)[::-1, :, :].swapaxes(0, 2), 3)
    if ax == 1:
        return np.rot90(m, 1)
    if ax == 2:
        return m.swapaxes(0, 2)[::-1, :, :]
def axial_rotations(m, rot=1, ax=2):


    """
    :param m: 3d matrix
    :param rot: number of rotations
    :param ax: axis of rotation
    :return: m rotate rot times around axis ax, according to convention.
    """

    if len(m.shape)!=3:

        assert IOError



    rot %= 4



    if rot == 0:
        return m
    for _ in range(rot):
        m = basic_rot_ax(m, ax=ax)
    return m
def rot90(m, k=1, axis=2):
    """Rotate an array k*90 degrees in the counter-clockwise direction around the given axis"""
    m = np.swapaxes(m, 2, axis)
    m = np.rot90(m, k)
    m = np.swapaxes(m, 2, axis)
    return m


