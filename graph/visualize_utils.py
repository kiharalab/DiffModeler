import numpy as np
from ops.os_operation import mkdir


def Show_Graph_Connect(coord_list,edge_list,save_path):
    Natm =1
    with open(save_path,'w') as file:
        file.write('MODEL\n')
        for i in range(len(coord_list)):
            line = ''
            tmp=coord_list[i]
            tmp_chain='A'
            line += "ATOM%7d  %3s %3s%2s%4d    " % (Natm, "CA ", "ALA", " " + tmp_chain, 1)
            line += "%8.3f%8.3f%8.3f%6.2f%6.8f\n" % (
             tmp[0], tmp[1], tmp[2], 1.0, 1.0)
            Natm += 1
            file.write(line)
        for i in range(len(edge_list)):
            nid1,nid2=edge_list[i]
            line = "BOND %d %d\n" % (nid1,nid2)
            file.write(line)
        file.write('ENDMDL\n')

def Show_Bfactor_cif(name,coord_list,save_path,base_prob):
    Natm =1
    chain_dict ={0:"A",1:"B",2:"C",3:"D",4:"E",5:"F",6:"G",7:"H"}
    with open(save_path,'w') as file:
        line = 'data_%s\n'%name
        line += "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" \
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"\
                    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"\
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"\
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"\
            "_atom_site.pdbx_PDB_model_num\n"
        file.write(line)
        for i in range(len(coord_list)):
            tmp=coord_list[i]
            current_prob = base_prob[i]
            tmp_chain="A"
            line =""
            line += "ATOM %-10d C %-3s . %-3s %-2s . %-10d  .  " % (Natm, "CA ", "ALA", " " + tmp_chain, Natm)
            line += "%-8.3f %-8.3f %-8.3f %-6.2f %-6.8f %-10d %-2s 1 \n" % (
             tmp[0], tmp[1], tmp[2], 1.0, current_prob,Natm,tmp_chain)
            Natm += 1
            file.write(line)
import os
def Show_Coord_cif(name,coord_list,save_path):
    Natm =1
    with open(save_path,'w') as file:
        line = 'data_%s\n'%name
        line += "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" \
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"\
                    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"\
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"\
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"\
            "_atom_site.pdbx_PDB_model_num\n"
        file.write(line)
        for i in range(len(coord_list)):
            tmp=coord_list[i]
            tmp_chain='A'
            line =""
            line += "ATOM %-10d C %-3s . %-3s %-2s . %-10d  .  " % (Natm, "CA ", "ALA", " " + tmp_chain, Natm)
            line += "%-8.3f %-8.3f %-8.3f %-6.2f %-6.8f %-10d %-2s 1 \n" % (
             tmp[0], tmp[1], tmp[2], 1.0, 1.0,Natm,tmp_chain)
            Natm += 1
            file.write(line)
def visualize_graph(save_path,ext_name, coordinate_list, edge_pairs):
    graph_path = os.path.join(save_path, ext_name+".cif")
    Show_Coord_cif( ext_name,coordinate_list,graph_path)
    if edge_pairs is None:
        return
    graph_path = os.path.join(save_path, ext_name+".pml")
    with open(graph_path,'w') as file:
        file.write("set connect_mode=1\n")
        file.write("load %s, TRACE\n"%(ext_name+".cif"))
        for k in range(len(edge_pairs)):
            id1,id2= edge_pairs[k]
            file.write("bond resi %d and TRACE, resi %d and TRACE\n"%(id1,id2))
        file.write("hide everything\n")
        file.write("show sticks, TRACE\n")
        file.write("set connect_mode=0\n")
def Show_BaseCoord_cif(name,coord_list,save_path,base_label,base_prob):
    Natm =1
    chain_dict ={0:"A",1:"U",2:"C",3:"G"}
    with open(save_path,'w') as file:
        line = 'data_%s\n'%name
        line += "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" \
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"\
                    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"\
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"\
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"\
            "_atom_site.pdbx_PDB_model_num\n"
        file.write(line)
        for i in range(len(coord_list)):
            tmp=coord_list[i]
            current_label = int(base_label[i])
            current_prob = base_prob[i]
            tmp_chain=chain_dict[current_label]
            line =""
            line += "ATOM %-10d C %-3s . %-3s %-2s . %-10d  .  " % (Natm, "CA ", "ALA", " " + tmp_chain, Natm)
            line += "%-8.3f %-8.3f %-8.3f %-6.2f %-6.8f %-10d %-2s 1 \n" % (
             tmp[0], tmp[1], tmp[2], 1.0, current_prob,Natm,tmp_chain)
            Natm += 1
            file.write(line)

def visualize_base(save_path,ext_name,coordinate_list,base_prediction):
    #base prediction will be [N,4] array each indicates the probability.
    path_path = os.path.join(save_path, ext_name+".cif")
    base_label = np.argmax(base_prediction,axis=1)
    base_prob = np.max(base_prediction,axis=1)
    Show_BaseCoord_cif(ext_name,coordinate_list,path_path,base_label,base_prob)


def Visualize_LDP_Path(save_path,ext_name,Path_ID_List,merged_cd_dens,map_info_list):
    from graph.LDP_ops import Convert_LDPcoord_To_Reallocation
    current_location_list = [merged_cd_dens[int(kk)] for kk in Path_ID_List]
    current_coord_list = Convert_LDPcoord_To_Reallocation(current_location_list, map_info_list)
    Show_Coord_cif(ext_name,current_coord_list,save_path)


def Visualize_assign_DPbase(save_dir,ext_name,path_node_id_list,base_assign_list,
               path_base_prob_list, merged_cd_dens,map_info_list):
    from graph.LDP_ops import Convert_LDPcoord_To_Reallocation#to avoid circular import
    current_location_list = [merged_cd_dens[int(kk)] for kk in path_node_id_list]
    coordinate_list= Convert_LDPcoord_To_Reallocation(current_location_list, map_info_list)
    new_coord_list=[]
    map_dict={"A":0,"T":1,"C":2,"G":3}
    for j in range(len(coordinate_list)-1):
        cur_coord = coordinate_list[j]
        next_coord = coordinate_list[j+1]
        new_coord =[(cur_coord[kk]+next_coord[kk])/2 for kk in range(3)]
        new_coord_list.append(new_coord)
    final_coord_list=[]
    base_label=[]
    base_prob_list=[]
    for j in range(len(new_coord_list)):
        current_chr = base_assign_list[j]
        if current_chr=="-":
            continue
        current_label= map_dict[current_chr]
        final_coord_list.append(new_coord_list[j])
        base_label.append(current_label)
        base_prob_list.append(path_base_prob_list[j][current_label])
    tmp_save_path = os.path.join(save_dir,ext_name+".cif")
    Show_BaseCoord_cif(ext_name,final_coord_list,tmp_save_path,base_label,base_prob_list)
