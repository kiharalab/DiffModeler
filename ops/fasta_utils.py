
from collections import  defaultdict
def read_fasta(input_fasta_path):
    #format should be
    #>chain_id, chain_id2, chain_id3
    #sequence
    chain_dict=defaultdict(list)#key: chain_list, value: nuc sequence
    current_id=None
    dna_rna_set= set(['A','U','T','C','G'])
    tmp_chain_list=[chr(i) for i in range(ord('A'), ord('Z') + 1)]  # uppercase letters
    tmp_chain_list.extend([chr(i) for i in range(ord('a'), ord('z') + 1)])  # lowercase letters
    tmp_chain_list+=['1','2','3','4','5','6','7','8','9','0']
    extend_chain_list=[]
    for i in range(len(tmp_chain_list)):
        for j in range(len(tmp_chain_list)):
            newid=tmp_chain_list[i]+tmp_chain_list[j]
            extend_chain_list.append(newid)
    tmp_chain_list+=extend_chain_list
    use_set=set()
    with open(input_fasta_path,'r') as file:
        for line in file:
            if line[0]==">":
                current_id = line.strip("\n")
                current_id = current_id.replace(">","")
                current_id = current_id.replace(" ","")
                #for users who do not follow our formats
                split_id_list = current_id.split(",")
                final_id=""
                for tmp_item in split_id_list:
                    if len(tmp_item)>2:
                        for tmp_chain_name in tmp_chain_list:
                            if tmp_chain_name not in use_set:
                                final_id+=tmp_chain_name+","
                                use_set.add(tmp_chain_name)
                                break
                        # final_id+=tmp_item+","
                        # use_set.add(tmp_item)
                    else:
                        if tmp_item not in use_set:
                            final_id+=tmp_item+","
                            use_set.add(tmp_item)
                        else:
                            for tmp_chain_name in tmp_chain_list:
                                if tmp_chain_name not in use_set:
                                    final_id+=tmp_chain_name+","
                                    use_set.add(tmp_chain_name)
                                    break

                    # if len(tmp_item)==1:
                    #     for tmp_chain_name in tmp_chain_list:
                    #         if tmp_chain_name not in use_set:
                    #             final_id+=tmp_chain_name+","
                    #             use_set.add(tmp_chain_name)
                    #             break
                    # elif len(tmp_item)==1:
                    #     final_id+=tmp_item+","
                    #     use_set.add(tmp_item)
                current_id =final_id[:-1]#remove last ","

            else:
                line=line.strip("\n").replace(" ","")
                tmp_resid_list=[]
                for item in line:
                    tmp_resid_list.append(item)
                    #chain_dict[current_id].append(item)
                tmp_resid_set = set(tmp_resid_list)
                tmp_useful_set = tmp_resid_set-dna_rna_set
                if len(tmp_useful_set)>0:
                    #may not be protein residues, skip DNA/RNA/ligand
                    chain_dict[current_id]+=tmp_resid_list
    print("read chain info from fasta:",chain_dict)
    return chain_dict

def write_fasta(fasta_list,use_chain_name,input_fasta_path):
    with open(input_fasta_path,'w') as file:
        file.write(">%s\n"%use_chain_name)
        for item in fasta_list:
            file.write(item)

def write_all_fasta(chain_dict,output_path):
    with open(output_path,'w') as wfile:
        for key in chain_dict:
            wfile.write(">%s\n"%key)
            fasta_list = chain_dict[key]
            for item in fasta_list:
                wfile.write(item)
            wfile.write("\n")

def refine_fasta_input(chain_dict,fitting_dict,refined_fasta_path):
    #check if chain in chain_dict is in fitting_dict
    #if not, remove it from chain_dict
    #if yes, write it to refined_fasta_path
    input_chain_list=list(chain_dict.keys())
    all_chain_list=[]
    for chain in input_chain_list:
        cur_chain_list = chain.split(",")
        for item in cur_chain_list:
            if len(item)>0:
                all_chain_list.append(item)
    print("all input chain list:",all_chain_list)
    use_chain_list=[]
    for fitting_path in fitting_dict:
        current_chain_list = fitting_dict[fitting_path] 
        for chain_id in current_chain_list:
            if chain_id in all_chain_list:
                use_chain_list.append(chain_id)
    print("use_chain_list:",use_chain_list)
    final_fitting_dict= {k:v for k,v in fitting_dict.items()}
    refined_chain_dict=defaultdict(list)
    for chain_list in chain_dict:
        cur_chain_list = chain_list.split(",")
        current_use_chain_list = []
        for item in cur_chain_list:
            if item in use_chain_list:
                current_use_chain_list.append(item)
        if len(current_use_chain_list)==0:
            refined_chain_dict[chain_list]=chain_dict[chain_list]
        else:
            if len(current_use_chain_list)==len(cur_chain_list):
                print("chain have been put into template",cur_chain_list)
            else:
                #some chains are indicated in the sequence but only one put into template, automatic match
                #add that missed chain into the template dict
                extra_chain_info = []
                for item in cur_chain_list:
                    if item not in current_use_chain_list:
                        extra_chain_info.append(item)
                print("extra unmatched chain info:",extra_chain_info)
                match_flag =False
                #update fitting dict
                for fitting_path in fitting_dict:
                    cur_fit_chain_list = fitting_dict[fitting_path]
                    cur_match_chain_list=[]
                    for item in cur_fit_chain_list:
                        if item in current_use_chain_list:
                            cur_match_chain_list.append(item)
                    #allow template to have more matches.
                    if len(cur_match_chain_list)==len(current_use_chain_list):
                        print("successfully match ",cur_match_chain_list,current_use_chain_list)
                        combine_chain_list = cur_chain_list+cur_fit_chain_list
                        combine_chain_list = list(set(combine_chain_list))
                        final_fitting_dict[fitting_path]=combine_chain_list
                        match_flag=True
                        break
                if match_flag is False:
                    print("can not find match for ",current_use_chain_list)
                    refined_chain_dict[chain_list]=chain_dict[chain_list]
    write_all_fasta(refined_chain_dict,refined_fasta_path)
    return final_fitting_dict

