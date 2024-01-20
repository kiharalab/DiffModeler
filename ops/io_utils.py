import os.path
import os
from collections import defaultdict
import shutil
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
def run_code_remote(command_line,save_dir):
    root_dir=os.getcwd()
    os.chdir(save_dir)
    os.system(command_line)
    os.chdir(root_dir)
import time
def download_file(url,file):
    try:
        import requests
        response = requests.get(url)

        with open(file, 'wb') as f:
          f.write(response.content)
    except:
        try:
            time.sleep(60)
            dir_name = os.path.dirname(os.path.realpath(file))
            file_name = os.path.basename(os.path.realpath(file))
            run_code_remote("wget %s"%url,dir_name)
            origin_file = os.path.join(dir_name,file_name)
            shutil.copy(origin_file,file)
        except:
            time.sleep(60)
            try:
                import requests
                response = requests.get(url)

                with open(file, 'wb') as f:
                  f.write(response.content)
            except:
                return False
    return True
