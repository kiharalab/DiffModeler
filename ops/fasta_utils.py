
from collections import  defaultdict
def read_fasta(input_fasta_path):
    #format should be
    #>chain_id, chain_id2, chain_id3
    #sequence
    chain_dict=defaultdict(list)#key: chain_list, value: nuc sequence
    current_id=None

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
                        final_id+=tmp_item+","
                        use_set.add(tmp_item)
                    else:
                        if tmp_item not in use_set:
                            final_id+=tmp_item+","
                        else:
                            for tmp_chain_name in tmp_chain_list:
                                if tmp_chain_name not in use_set:
                                    final_id+=tmp_chain_name+","
                                    use_set.add(tmp_chain_name)
                                    break
                        use_set.add(tmp_item)
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
                for item in line:
                    chain_dict[current_id].append(item)
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
