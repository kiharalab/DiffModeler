import os.path
from collections import defaultdict
def read_structure_txt(input_dir,input_file_path):
    structure_dict={}
    with open(input_file_path,'r') as rfile:
        for line in rfile:
            line = line.strip("\n")
            split_info= line.split()
            input_file_path = os.path.join(input_dir,split_info[0])
            structure_dict[input_file_path]=split_info[1:]
    print("structure waiting to be fitted: ",structure_dict)
    return structure_dict
import pickle
def load_pickle(path):
    with open(path,'rb') as file:
        data=pickle.load(file)
    return data

def write_pickle(data,path):
    with open(path,'wb') as file:
        pickle.dump(data, file)


def download_file(url,file):
    import requests
    response = requests.get(url)

    with open(file, 'wb') as f:
      f.write(response.content)
