
from collections import  defaultdict
def read_fasta(input_fasta_path):
    #format should be
    #>chain_id, chain_id2, chain_id3
    #sequence
    chain_dict=defaultdict(list)#key: chain_list, value: nuc sequence
    current_id=None

    tmp_chain_list=[chr(i) for i in range(ord('A'), ord('Z') + 1)]  # uppercase letters
    tmp_chain_list.extend([chr(i) for i in range(ord('a'), ord('z') + 1)])  # lowercase letters

    with open(input_fasta_path,'r') as file:
        for line in file:
            if line[0]==">":
                current_id = line.strip("\n")
                current_id = current_id.replace(">","")
                current_id = current_id.replace(" ","")

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
