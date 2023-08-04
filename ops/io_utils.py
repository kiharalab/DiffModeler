import os.path
from collections import defaultdict
def read_structure_txt(input_file_path):
    structure_dict={}
    with open(input_file_path,'r') as rfile:
        for line in rfile:
            line = line.strip("\n")
            split_info= line.split()
            input_file_path = os.path.abspath(split_info[0])
            structure_dict[input_file_path]=split_info[1:]
    print("structure waiting to be fitted: ",structure_dict)
    return structure_dict
