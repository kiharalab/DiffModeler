import os 
import pickle
import shutil
from ops.io_utils import load_json
from ops.pdb_utils import filter_pdb_by_residues
#call sword2 to get different domains and then run diffmodeler with domains
def prepare_domain_input(fitting_dict,save_dir,num_cpu=8):
    """
    :param fitting_dict:    [pdb_path]:[list of chain_id]
    :param save_dir:
    return: 
    fitting_dict
    """
    domain_fitting_dict_path = os.path.join(save_dir,"domain_fitting_dict.pkl")
    if os.path.exists(domain_fitting_dict_path):
        fitting_dict = pickle.load(open(domain_fitting_dict_path,"rb"))
        #check every path exists
        verify_fitting_dict = True
        for pdb_path in fitting_dict:
            if not os.path.exists(pdb_path):
                verify_fitting_dict = False
                break
        if verify_fitting_dict:
            return fitting_dict
        
    current_dir = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.path.join(current_dir,"SWORD2")
    sword_script_path = os.path.join(current_dir,"SWORD2.py")
    #clean the save_dir
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    os.makedirs(save_dir,exist_ok=True)
    domain_fitting_dict = {}

    for pdb_path in fitting_dict:
        chain_ids = fitting_dict[pdb_path]
        pdb_name = os.path.basename(pdb_path).replace(".pdb","")
        pdb_dir = os.path.join(save_dir,pdb_name)
        os.makedirs(pdb_dir,exist_ok=True)
        #split pdb into different domains
        command_line=f"python3 {sword_script_path} -i {pdb_path} -o {pdb_dir} -c A --cpu {num_cpu}"
        os.system(command_line)
        cur_output_dir=os.path.join(pdb_dir,f"{pdb_name}_A")
        if not os.path.exists(cur_output_dir):
            print(f"Error: SWORD2 failed to split {pdb_path}")
            continue
        cur_output_json = os.path.join(cur_output_dir,"SWORD2_summary.json")
        domain_info = load_json(cur_output_json)
        boundary_info = domain_info['Optimal partition']['BOUNDARIES']
        print(f"Domain info for {pdb_path}: {boundary_info}")
        for domain_name in boundary_info:
            current_domain = boundary_info[domain_name]
            #filter all residues
            residue_index_list=[]
            for tmp_range in current_domain:
                residue_index_list+=list(range(tmp_range[0],tmp_range[1]+1))
            #output the domain pdb
            #    
            domain_pdb_path = os.path.join(save_dir,f"{pdb_name}_{domain_name}.pdb")
            filter_pdb_by_residues(pdb_path,domain_pdb_path,residue_index_list)
            current_chain_ids = [x+str(domain_name) for x in chain_ids]
            domain_fitting_dict[domain_pdb_path]=current_chain_ids
            print(f"Domain {domain_name} for {pdb_path} is saved at {domain_pdb_path}")
    pickle.dump(domain_fitting_dict,open(domain_fitting_dict_path,"wb"))

    return domain_fitting_dict
