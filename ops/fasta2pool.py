import os
from ops.os_operation import mkdir
from ops.io_utils import download_file
from ops.pdb_utils import count_atom_line
def fasta2pool(params,save_path):
    single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
    mkdir(single_chain_pdb_dir)
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
        fasta_list = chain_name_list[chain_dict]
        chain_name_list = chain_name_list.replace(",","-")
        current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list))
        mkdir(current_chain_dir)
        input_fasta_path = os.path.join(current_chain_dir,"input.fasta")
        use_chain_name = chain_name_list.split("-")[0]
        write_fasta(fasta_list,use_chain_name,input_fasta_path)
        command_line="python %s --email %s --program fasta " \
                     "--stype protein --database pdb,afdb " \
                     "--sequence %s"%(search_script,params['email'],input_fasta_path)
        pool.apply_async(run_command,args=(command_line,))
    pool.close()
    pool.join()

    final_pdb_dir = os.path.join(single_chain_pdb_dir,"PDB")
    mkdir(final_pdb_dir)
    #after blocking finished, extract the top 1 id and fetch corressponding pdb from pdb/afdb

    for chain_name_list in chain_dict:

        chain_name_list = chain_name_list.replace(",","-")
        current_chain_dir = os.path.join(single_chain_pdb_dir,str(chain_name_list))
        listfiles = [x for x in os.listdir(current_chain_dir) if ".ids.txt" in x]
        if len(listfiles)==0:
            print("fail to find search results for chain %s"%chain_name_list)
            print("-"*20+" search again "+"-"*20)
            input_fasta_path = os.path.join(current_chain_dir,"input.fasta")
            command_line="python %s --email %s --program fasta " \
                     "--stype protein --database pdb,afdb " \
                     "--sequence %s"%(search_script,params['email'],input_fasta_path)
            run_command(command_line)
        cur_file = os.path.join(current_chain_dir,listfiles[0])
        with open(cur_file,'r') as rfile:
            line=rfile.readline()
            line = line.strip("\n")
            split_info = line.split(":")
            database = split_info[0]
            pdb_id = split_info[1]
        if database=="PDB":
            pdb = pdb_id.split("_")[0]
            chain_id = pdb_id.split("_")[1]
            download_link = "https://files.rcsb.org/download/%s.cif"%pdb
            cif_file = os.path.join(current_chain_dir,"input.cif")
            download_file(download_link,cif_file)
            if os.path.exists(cif_file) and count_atom_line(cif_file)>=50:
                #sgement the specific chains

            else:

