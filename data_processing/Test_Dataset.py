
import os
import numpy as np
import torch
import torch.utils.data
import random

class Single_Dataset(torch.utils.data.Dataset):
    def __init__(self, data_path, search_key="input_"):
        """
        :param data_path: training data path
        :param training_id: specify the id that will be used for training
        """
        listfiles = [x for x in os.listdir(data_path) if search_key in x]
        listfiles.sort()
        self.input_path = []
        self.id_list = []
        for x in listfiles:
            self.input_path.append(os.path.join(data_path,x))
            cur_id = int(x.replace(search_key,"").replace(".npy",""))
            self.id_list.append(cur_id)


    def __len__(self):
        return len(self.input_path)

    def __getitem__(self, idx):
        inputfile = self.input_path[idx]
        cur_id = self.id_list[idx]
        input = np.load(inputfile)
        input = input[np.newaxis, :]
        input = torch.from_numpy(np.array(input, np.float32, copy=True))
        return {'density': input, 'Index': idx,"ID":cur_id}
