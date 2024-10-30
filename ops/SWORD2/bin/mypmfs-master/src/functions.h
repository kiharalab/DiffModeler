/*
 * functions.h
 *
 * Functions used by both training.cpp and scoring.cpp
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2017 Guillaume Postic (guillaume.postic@upmc.fr)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */



std::map<std::string,std::string> create_map(){
    std::map<std::string,std::string> three2one;
    three2one["ALA"]="A"; three2one["ARG"]="R"; three2one["ASN"]="N"; three2one["ASP"]="D";
    three2one["CYS"]="C"; three2one["GLU"]="E"; three2one["GLN"]="Q"; three2one["GLY"]="G";
    three2one["HIS"]="H"; three2one["ILE"]="I"; three2one["LEU"]="L"; three2one["LYS"]="K";
    three2one["MET"]="M"; three2one["PHE"]="F"; three2one["PRO"]="P"; three2one["SER"]="S";
    three2one["THR"]="T"; three2one["TRP"]="W"; three2one["TYR"]="Y"; three2one["VAL"]="V";
    return three2one;
}


std::set<std::string> set_allatom(){
    // The 167 atom types
    std::set<std::string> atypes = {"AC", "ACA", "ACB", "AN", "AO", "CC", "CCA", "CCB", "CN", "CO", "CSG", "DC", "DCA", "DCB", "DCG", "DN", "DO", "DOD1", "DOD2", "EC", "ECA", "ECB", "ECD", "ECG", "EN", "EO", "EOE1", "EOE2", "FC", "FCA", "FCB", "FCD1", "FCD2", "FCE1", "FCE2", "FCG", "FCZ", "FN", "FO", "GC", "GCA", "GN", "GO", "HC", "HCA", "HCB", "HCD2", "HCE1", "HCG", "HN", "HND1", "HNE2", "HO", "IC", "ICA", "ICB", "ICD1", "ICG1", "ICG2", "IN", "IO", "KC", "KCA", "KCB", "KCD", "KCE", "KCG", "KN", "KNZ", "KO", "LC", "LCA", "LCB", "LCD1", "LCD2", "LCG", "LN", "LO", "MC", "MCA", "MCB", "MCE", "MCG", "MN", "MO", "MSD", "NC", "NCA", "NCB", "NCG", "NN", "NND2", "NO", "NOD1", "PC", "PCA", "PCB", "PCD", "PCG", "PN", "PO", "QC", "QCA", "QCB", "QCD", "QCG", "QN", "QNE2", "QO", "QOE1", "RC", "RCA", "RCB", "RCD", "RCG", "RCZ", "RN", "RNE", "RNH1", "RNH2", "RO", "SC", "SCA", "SCB", "SN", "SO", "SOG", "TC", "TCA", "TCB", "TCG2", "TN", "TO", "TOG1", "VC", "VCA", "VCB", "VCG1", "VCG2", "VN", "VO", "WC", "WCA", "WCB", "WCD1", "WCD2", "WCE2", "WCE3", "WCG", "WCH2", "WCZ2", "WCZ3", "WN", "WNE1", "WO", "YC", "YCA", "YCB", "YCD1", "YCD2", "YCE1", "YCE2", "YCG", "YCZ", "YN", "YO", "YOH"};
    return atypes;
}


std::set<std::string> set_allatomCG(){
    // The 49 coarse-grained atom types
    std::set<std::string> atypes = {"ABB", "CBB", "CSC1", "DBB", "DSC1", "EBB", "ESC1", "FBB", "FSC1", "FSC2", "FSC3", "GBB", "HBB", "HSC1", "HSC2", "HSC3", "IBB", "ISC1", "KBB", "KSC1", "KSC2", "LBB", "LSC1", "MBB", "MSC1", "NBB", "NSC1", "PBB", "PSC1", "QBB", "QSC1", "RBB", "RSC1", "RSC2", "SBB", "SSC1", "TBB", "TSC1", "VBB", "VSC1", "WBB", "WSC1", "WSC2", "WSC3", "WSC4", "YBB", "YSC1", "YSC2", "YSC3"};
    return atypes;
}


std::set<std::string> set_backbone(){
    // The 80 backbone atom types
    std::set<std::string> atypes = {"AN", "ACA", "AC", "AO", "CN", "CCA", "CC", "CO", "DN", "DCA", "DC", "DO", "EN", "ECA", "EC", "EO", "FN", "FCA", "FC", "FO", "GN", "GCA", "GC", "GO", "HN", "HCA", "HC", "HO", "IN", "ICA", "IC", "IO", "KN", "KCA", "KC", "KO", "LN", "LCA", "LC", "LO", "MN", "MCA", "MC", "MO", "NN", "NCA", "NC", "NO", "PN", "PCA", "PC", "PO", "QN", "QCA", "QC", "QO", "RN", "RCA", "RC", "RO", "SN", "SCA", "SC", "SO", "TN", "TCA", "TC", "TO", "VN", "VCA", "VC", "VO", "WN", "WCA", "WC", "WO", "YN", "YCA", "YC", "YO"};
    return atypes;
}


std::set<std::string> set_sidechains(){
    // The 88 side chain atom types
    std::set<std::string> atypes = {"ACB", "CCB", "CSG", "DCB", "DCG", "DOD1", "DOD2", "ECB", "ECD", "ECG", "EOE1", "EOE2", "FCB", "FCD1", "FCD2", "FCE1", "FCE2", "FCG", "FCZ", "GCA", "HCB", "HCD2", "HCE1", "HCG", "HND1", "HNE2", "ICB", "ICD1", "ICG1", "ICG2", "KCB", "KCD", "KCE", "KCG", "KNZ", "LCB", "LCD1", "LCD2", "LCG", "MCB", "MCE", "MCG", "MSD", "NCB", "NCG", "NND2", "NOD1", "PCB", "PCD", "PCG", "QCB", "QCD", "QCG", "QNE2", "QOE1", "RCB", "RCD", "RCG", "RCZ", "RNE", "RNH1", "RNH2", "SCB", "SOG", "TCB", "TCG2", "TOG1", "VCB", "VCG1", "VCG2", "WCB", "WCD1", "WCD2", "WCE2", "WCE3", "WCG", "WCH2", "WCZ2", "WCZ3", "WNE1", "YCB", "YCD1", "YCD2", "YCE1", "YCE2", "YCG", "YCZ", "YOH"};
    return atypes;
}


std::set<std::string> set_sidechainsCG(){
    // The 31 coarse-grained side chain atom types
    std::set<std::string> atypes = {"ABB", "CSC1", "DSC1", "ESC1", "FSC1", "FSC2", "FSC3", "GBB", "HSC1", "HSC2", "HSC3", "ISC1", "KSC1", "KSC2", "LSC1", "MSC1", "NSC1", "PSC1", "QSC1", "RSC1", "RSC2", "SSC1", "TSC1", "VSC1", "WSC1", "WSC2", "WSC3", "WSC4", "YSC1", "YSC2", "YSC3"};
    return atypes;
}


std::set<std::string> set_CA(){
    // The 20 CA atom types
    std::set<std::string> atypes = {"ACA", "CCA", "DCA", "ECA", "FCA", "GCA", "HCA", "ICA", "KCA", "LCA", "MCA", "NCA", "PCA", "QCA", "RCA", "SCA", "TCA", "VCA", "WCA", "YCA"};
    return atypes;
}


std::set<std::string> set_CB(){
    // The 20 CB atom types
    std::set<std::string> atypes = {"ACB", "CCB", "DCB", "ECB", "FCB", "GCA", "HCB", "ICB", "KCB", "LCB", "MCB", "NCB", "PCB", "QCB", "RCB", "SCB", "TCB", "VCB", "WCB", "YCB"};
    return atypes;
}


std::set<std::string> set_BB(){
    // The 20 BB atom types
    std::set<std::string> atypes = {"ABB", "CBB", "DBB", "EBB", "FBB", "GBB", "HBB", "IBB", "KBB", "LBB", "MBB", "NBB", "PBB", "QBB", "RBB", "SBB", "TBB", "VBB", "WBB", "YBB"};
    return atypes;
}


std::set<std::string> set_SC1(){
    // The 20 SC1 atom types
    std::set<std::string> atypes = {"ABB", "CSC1", "DSC1", "ESC1", "FSC1", "GBB", "HSC1", "ISC1", "KSC1", "LSC1", "MSC1", "NSC1", "PSC1", "QSC1", "RSC1", "SSC1", "TSC1", "VSC1", "WSC1", "YSC1"};
    return atypes;
}


std::vector<std::string> file2vec(std::string filename, std::map<std::string,std::string> three2one, std::set<std::string> atypes){
    // Extract from the PDB file the lines that are useful for other functions
    std::ifstream fh(filename);
    if(!fh){ std::cerr << "Error while opening file " << filename << std::endl; exit(1); }

    std::vector<std::string> filevec;
    std::string line;
    while(getline(fh, line) and line.substr(0,6) != "ENDMDL"){ // will only consider the first MODEL
        if (line.size() > 53 and line.substr(0,4) == "ATOM"){
            std::string atomname = line.substr(12,4);
            atomname.erase(std::remove(atomname.begin(),atomname.end(),' '),atomname.end());
            std::string threeletter = line.substr(17,3);

            // if not unusual amino acid
            if (three2one.count(threeletter))
                // if is one of the atom types
                if (atypes.find(three2one[threeletter]+atomname) != atypes.end())
                    filevec.push_back(line);
        }
    }

    return filevec;
}


std::vector<std::string> findchains(std::vector<std::string>& filevec){
    // Identifies the protein chain(s) in the PDB file
    // Chains are stored in a std::set<std::string>
    std::set<std::string> chains;
    for (auto line : filevec){
        chains.insert(line.substr(21,1));
    }

    // Store chains into a std::vector
    std::vector<std::string> chainlist;
    for (auto chain : chains){
        chainlist.push_back(chain);
    }

    return chainlist;
}


std::string exec(std::string cmd) {
    // Returns the STDOUT as a std::string, including '\n'!
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}


bool fexists(std::string filename){
  std::ifstream ifile(filename);
  return (bool)ifile;
}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

