
import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
import shutil
import sis

def init_save_path(origin_map_path):
    save_path = os.path.join(os.getcwd(), 'Predict_Result')
    mkdir(save_path)
    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    map_name = map_name.replace(".map", "")
    map_name = map_name.replace("(","").replace(")","")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    return save_path,map_name

def set_up_envrionment(params):
    if params['resolution']>20:
        print("maps with %.2f resolution is not supported! We only support maps with resolution 0-20A!"%params['resolution'])
        exit()
    gpu_id = params['gpu']
    if gpu_id is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
    cur_map_path = os.path.abspath(params['F'])
    if cur_map_path.endswith(".gz"):
        from ops.os_operation import unzip_gz
        cur_map_path = unzip_gz(cur_map_path)

    if params['output'] is None:
        save_path,map_name = init_save_path(cur_map_path)
    else:
        save_path=params['output']
        map_name="input_diffmodeler" #to avoid server same name bugs
        mkdir(save_path)
    try:
        print("pre-compile VESPER to accelerate!")
        running_dir = os.path.dirname(os.path.abspath(__file__))
        os.system(f"cd {running_dir}; python -O -m compileall VESPER_CUDA")
    except:
        print("pre-compile VESPER failed! No impact to main scripts!")
    save_path = os.path.abspath(save_path)
    from data_processing.Unify_Map import Unify_Map
    cur_map_path = Unify_Map(cur_map_path,os.path.join(save_path,map_name+"_unified.mrc"))
    from data_processing.Resize_Map import Resize_Map
    cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
    if params['contour']<0:
        #change contour level to 0 and increase all the density
        from ops.map_utils import increase_map_density
        cur_map_path = increase_map_density(cur_map_path,os.path.join(save_path,map_name+"_increase.mrc"),params['contour'])
        params['contour']=0
    from modeling.map_utils import segment_map
    new_map_path = os.path.join(save_path,map_name+"_segment.mrc")
    segment_map(cur_map_path,new_map_path,contour=0)
    return save_path,new_map_path

def diffusion_trace_map(save_path,cur_map_path,params):
    if params['resolution']>=2:
        from predict.infer_diffusion import infer_diffem
        diffusion_dir = os.path.join(save_path,"infer_diffusion")
        diff_trace_map=infer_diffem(cur_map_path,diffusion_dir,params)
    else:
        print("skip diffusion with very high resolution map %f"%params['resolution'])
        diff_trace_map = cur_map_path
    print("Diffusion process finished! Traced map saved here %s"%diff_trace_map)
    # segment this difftrace map to save time
    from modeling.map_utils import segment_map
    diff_new_trace_map = os.path.join(save_path,"diffusion.mrc")
    segment_map(diff_trace_map,diff_new_trace_map,contour=0)
    return diff_new_trace_map
from ops.io_utils import delete_dir
def construct_single_chain_candidate(params,save_path):
    #first build a dict from the input text configure file
    single_chain_pdb_input = os.path.abspath(params['P'])
    single_chain_pdb_dir = os.path.join(save_path,"single_chain_pdb")
    from ops.pdb_utils import cif2pdb
    #if os.path.exists(single_chain_pdb_dir):
        #shutil.rmtree(single_chain_pdb_dir)
    #delete_dir(single_chain_pdb_dir)
    if not os.path.isdir(single_chain_pdb_input):
        os.makedirs(single_chain_pdb_dir,exist_ok=True)
        from ops.os_operation import extract_compressed_file
        single_chain_pdb_dir=extract_compressed_file(single_chain_pdb_input,single_chain_pdb_dir)
    else:
        if os.path.exists(single_chain_pdb_dir):
            shutil.rmtree(single_chain_pdb_dir)
        from ops.os_operation import copy_directory
        single_chain_pdb_dir = copy_directory(single_chain_pdb_input,single_chain_pdb_dir)
    #for every .cif files in the single_chain_pdb_dir, convert them to pdb
    for file in os.listdir(single_chain_pdb_dir):
        if file.endswith(".cif"):
            cur_cif_path = os.path.join(single_chain_pdb_dir,file)
            cur_pdb_path = os.path.join(single_chain_pdb_dir,file.replace(".cif",".pdb"))
            cif2pdb(cur_cif_path,cur_pdb_path)
    from ops.io_utils import read_structure_txt
    fitting_dict = read_structure_txt(single_chain_pdb_dir,os.path.abspath(params['M']))
    return fitting_dict

def fix_cif_for_coot(input_cif, output_cif):
    """
    Reads a CIF file and ensures:
    - '_atom_site.auth_asym_id' is copied from '_atom_site.label_asym_id'.
    - '_atom_site.auth_seq_id' is copied from '_atom_site.label_seq_id'.
    - All atom entries align properly in the CIF format for Coot.
    """

    with open(input_cif, 'r') as infile:
        lines = infile.readlines()

    new_lines = []
    in_atom_site = False
    headers = []
    modified_headers = False

    for line in lines:
        stripped_line = line.strip()

        # Detect start of _atom_site loop
        if stripped_line.startswith("loop_"):
            in_atom_site = False  # Reset detection
        if stripped_line.startswith("_atom_site."):
            in_atom_site = True
            headers.append(stripped_line)

        # Modify header to include '_atom_site.auth_asym_id' and '_atom_site.auth_seq_id'
        if in_atom_site and not modified_headers and "_atom_site.label_seq_id" in stripped_line:
            if "_atom_site.auth_asym_id" not in headers:
                headers.append("_atom_site.auth_asym_id")
            if "_atom_site.auth_seq_id" not in headers:
                headers.append("_atom_site.auth_seq_id")
            modified_headers = True
            continue  # Skip writing this line since we'll rewrite the headers later

        # Ensure atom data has correct column count
        if in_atom_site and stripped_line.startswith("ATOM"):
            parts = stripped_line.split()
            if len(parts) == len(headers) - 2:  # Missing two columns
                label_asym_id = parts[5]  # Chain ID
                label_seq_id = parts[6]  # Residue number
                parts.insert(7, label_asym_id)  # Insert _atom_site.auth_asym_id
                parts.insert(8, label_seq_id)  # Insert _atom_site.auth_seq_id
            elif len(parts) != len(headers):  # If still inconsistent, print warning
                print(f"WARNING: Inconsistent CIF loop at line:\n {stripped_line}")

            new_line = "  ".join(parts) + "\n"
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    # Write fixed headers at the correct place
    with open(output_cif, 'w') as outfile:
        for line in new_lines:
            if line.strip().startswith("_atom_site."):
                # Write all headers once at the correct place
                if headers:
                    outfile.write("\n".join(headers) + "\n")
                    headers = []  # Clear headers so they aren't written again
                continue
            outfile.write(line)

    print(f"Reformatted CIF file saved as: {output_cif}")



if __name__ == "__main__":
    params = argparser()
    save_path,cur_map_path = set_up_envrionment(params)
    running_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(running_dir)
    if not os.path.isabs(params['model']['path']):
        params['model']['path'] =os.path.join(running_dir,params['model']['path'])
    if not os.path.isabs(params['db_exp_path']):
        params['db_exp_path'] = os.path.join(running_dir,params['db_exp_path'])
    if not os.path.isabs(params['db_path']):
        params['db_path'] = os.path.join(running_dir,params['db_path'])
    if params['mode']==0:
        #first build a dict from the input text configure file
        fitting_dict = construct_single_chain_candidate(params,save_path)
    elif params['mode']==1:
        from ops.fasta2pool import fasta2pool
        fitting_dict = fasta2pool(params,save_path)
    elif params['mode']==2:
        from ops.fasta_searchdb import fasta_searchdb
        fitting_dict = fasta_searchdb(params,save_path)
        if params['seq_search']:
            print("sequence search finished!")
            exit()
    elif params['mode']==3:
        #support fasta+template mixed mode. Will first check if template is provided, will use the provided one first. Then for remained results, do db search
        fitting_dict = construct_single_chain_candidate(params,save_path)
        params['fasta_path'] = os.path.abspath(params['fasta_path'])
        from ops.fasta_utils import read_fasta
        chain_dict = read_fasta(params['fasta_path'])
        refined_fasta_path = os.path.join(save_path,"refined_input.fasta")
        from ops.fasta_utils import refine_fasta_input
        fitting_dict = refine_fasta_input(chain_dict,fitting_dict,refined_fasta_path)
        print("updated fitting_dict:",fitting_dict)
        params['P']=refined_fasta_path
        from ops.fasta_searchdb import fasta_searchdb
        additional_fitting_dict = fasta_searchdb(params,save_path)
        #merge two fitting dict
        for key in additional_fitting_dict:
            fitting_dict[key]=additional_fitting_dict[key]
        print("final fitting dict:",fitting_dict)
    else:
        print("mode %d is not supported!"%params['mode'])
        exit()

    if len(fitting_dict)==0:
        print("Empty Template candiate, DiffModeler can not run!!!")
        exit()
    #clean fitting dict to avoid strange pdb cause entire program fail
    final_template_dir = os.path.join(save_path,"final_template_input")
    from ops.pdb_utils import clean_pdb_template
    final_fitting_dict=clean_pdb_template(fitting_dict,final_template_dir)

    if params['domain']:
        from ops.domain_utils import prepare_domain_input
        domain_template_dir = os.path.join(save_path,"domain_template_input")
        final_fitting_dict = prepare_domain_input(final_fitting_dict,
                                            domain_template_dir,
                                            num_cpu=params['SWORD_thread'])
        print("Domain split finished!",final_fitting_dict)
    #diffusion inference
    diff_trace_map = diffusion_trace_map(save_path,cur_map_path,params)

    #VESPER singl-chain fitting process
    fitting_dir = os.path.join(save_path,"structure_modeling")
    from modeling.fit_structure_chain import fit_structure_chain
    fit_structure_chain(diff_trace_map,final_fitting_dict,fitting_dir,params)

    #VESPER assembling
    modeling_dir = os.path.join(save_path,"structure_assembling")
    from modeling.assemble_structure import assemble_structure
    source_cif = assemble_structure(diff_trace_map,final_fitting_dict,fitting_dir,modeling_dir,params)
    output_cif = os.path.join(save_path,"DiffModeler_alpha.cif")
    shutil.copy(source_cif,output_cif)
    # convert the cif format
    input_cif = output_cif
    output_cif = os.path.join(save_path, "DiffModeler.cif")
    fix_cif_for_coot(input_cif, output_cif)
    
    #generate a cif file to save the fitting score to b-factor field for easier visualization
    #for server visualization on server
    from modeling.pdb_utils import swap_cif_occupancy_bfactor
    score_specific_path = os.path.join(save_path,"DiffModeler_fitscore.cif")
    swap_cif_occupancy_bfactor(output_cif,score_specific_path)

    print(f"Please check DiffModeler's output structure in {output_cif}")





