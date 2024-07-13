import os
from ops.os_operation import mkdir,functToDeleteItems
from ops.io_utils import download_file
from ops.pdb_utils import count_atom_line,filter_chain_cif,cif2pdb,filter_chain_pdb,count_residues
from ops.fasta_searchdb import download_pdb

def parse_template_info(input_file):
    candidate_list=[]
    with open(input_file,'r') as rfile:
        for kk in range(10):
            line=rfile.readline()
            line = line.strip("\n")
            split_info = line.split(":")
            database = split_info[0]
            pdb_id = split_info[1]
            candidate_list.append([database,pdb_id])
    return candidate_list
def fasta2pool(params,save_path):
    max_file_name_limit=250#default is 255, to be safe
    single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
    mkdir(single_chain_pdb_dir)

    final_pdb_dir = os.path.join(single_chain_pdb_dir,"PDB")
    mkdir(final_pdb_dir)
    fasta_path = os.path.abspath(params['P'])
    from ops.fasta_utils import read_fasta,write_fasta
    chain_dict = read_fasta(fasta_path)
    search_script= os.path.join(os.getcwd(),"ops")
    search_script = os.path.join(search_script,"fasta_to_similar_pdb.py")
    print("start fetching pdb from the database with fasta sequence information.")
    from multiprocessing import Pool
    from ops.os_operation import run_command
    pool = Pool(min(params["fasta_thread"],len(chain_dict)))
    for chain_name_list in chain_dict:
        fasta_list = chain_dict[chain_name_list]
        chain_name_list = chain_name_list.replace(",","-")
        final_pdb_path = os.path.join(final_pdb_dir,chain_name_list[:max_file_name_limit]+".pdb")
        if os.path.exists(final_pdb_path) and count_atom_line(final_pdb_path)>=50:
            continue
        current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list[:max_file_name_limit]))
        mkdir(current_chain_dir)
        input_fasta_path = os.path.join(current_chain_dir,"input.fasta")
        use_chain_name = chain_name_list.split("-")[0]
        write_fasta(fasta_list,use_chain_name,input_fasta_path)
        command_line="cd %s; python %s --email %s --program fasta " \
                     "--stype protein --database pdb,afdb " \
                     "--sequence %s"%(current_chain_dir,search_script,params['email'],input_fasta_path)
        pool.apply_async(run_command,args=(command_line,))
    pool.close()
    pool.join()


    #after blocking finished, extract the top 1 id and fetch corressponding pdb from pdb/afdb
    fitting_dict={}
    for chain_name_list in chain_dict:

        chain_name_list = chain_name_list.replace(",","-")
        final_pdb_path = os.path.join(final_pdb_dir,chain_name_list[:max_file_name_limit]+".pdb")
        if os.path.exists(final_pdb_path) and count_atom_line(final_pdb_path)>=50:
            final_chain_list = chain_name_list.split("-")
            fitting_dict[final_pdb_path]=final_chain_list
            continue
        current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list[:max_file_name_limit]))
        listfiles = [x for x in os.listdir(current_chain_dir) if ".ids.txt" in x]
        if len(listfiles)==0:
            print("fail to find search results for chain %s"%chain_name_list)
            print("-"*20+" search again "+"-"*20)
            input_fasta_path = os.path.join(current_chain_dir,"input.fasta")
            command_line="cd %s; python %s --email %s --program fasta " \
                     "--stype protein --database pdb,afdb " \
                     "--sequence %s"%(current_chain_dir,search_script,
                                      params['email'],input_fasta_path)
            run_command(command_line)
        cur_file = os.path.join(current_chain_dir,listfiles[0])
        candidate_list = parse_template_info(cur_file)
        listfiles_acc = [x for x in os.listdir(current_chain_dir) if ".accs.txt" in x]
        if len(listfiles_acc)==0:
            print("fail to find search results for chain %s *.accs.txt from EBI-Search!"%chain_name_list)
            exit()
        cur_acc_file = os.path.join(current_chain_dir,listfiles_acc[0])
        pdb_candidate_list = parse_template_info(cur_acc_file)
        
        find_flag=False
        for candidate_index,candidate in enumerate(candidate_list):
            database,pdb_id=candidate
            print("candidate %s: %s"%(database,pdb_id))
            tmp_cur_chain_dir = os.path.join(current_chain_dir,str(candidate_index))
            mkdir(tmp_cur_chain_dir)
            if database=="PDB":
                _,pdb_candidate = pdb_candidate_list[candidate_index]
                download_pdb(pdb_candidate,tmp_cur_chain_dir,final_pdb_path)

            else:
                #alphafold db
                download_link = "https://alphafold.ebi.ac.uk/files/%s-model_v4.pdb"%pdb_id
                download_file(download_link,final_pdb_path)
            expected_seq_length = len(chain_dict[chain_name_list.replace("-",",")])*params['search']['length_ratio']
            actual_structure_length = count_residues(final_pdb_path)
            max_allow_length = len(chain_dict[chain_name_list.replace("-",",")])*params['search']['max_length_ratio']
            if actual_structure_length>=expected_seq_length and \
                    actual_structure_length<=max_allow_length:
                find_flag=True
                break
            else:
                print("structure length %d is not in the range of %d and %d"%(actual_structure_length,expected_seq_length,max_allow_length))
                functToDeleteItems(tmp_cur_chain_dir)
                os.remove(final_pdb_path)
        if find_flag is False:
            database,pdb_id=candidate_list[0]
            _,pdb_candidate = pdb_candidate_list[0]
            download_pdb(pdb_candidate,current_chain_dir,final_pdb_path)
        final_chain_list = chain_name_list.split("-")
        fitting_dict[final_pdb_path]=final_chain_list
    print("collecting finish: fitting dict: ",fitting_dict)
    return fitting_dict
