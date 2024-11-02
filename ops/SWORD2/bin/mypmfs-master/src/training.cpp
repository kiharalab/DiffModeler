/*
 * training.cpp
 *
 * The training part of the MyPMFs suite
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


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <getopt.h>
#include <cmath>
#include <iomanip> // setprecision
#include <sys/stat.h>
#include "functions.h"
#include <cctype>

using namespace std;

void squaredist(string pdbcode, vector<string>& filevec, vector<float>& sqdist, vector<string>& atompairs, vector<string>& pairID,
                  double distmax, double distmin, string chain1, string chain2, int diffmin, int diffmax, map<string,string> three2one, bool allatom);
bool sortbysec(const pair<string,double> &a, const pair<string,double> &b);
string makeupper(string s);
string makelower(string s);
string findfile(string code, string inputdir);
void f_writeRefFreq(map<string, vector<double>>& frequencies, string jobdir);

// Max file size (distance distribution) on which R KDE can be applied
const int MAX_SIZE_KDE = 1000; // Recommended value for the average system = 1000
// Same for the xx distribution (i.e. when not accounting for the residue type)
const int MAX_SIZE_KDEXX = 4500; // Recommended value for the average system = 4500

int main(int argc, char** argv){
    string optlist =
        "    Usage:\n"
        "    ./training {-l INPUT_LIST_NATIVE | -L INPUT_LIST} -d PDB_FILES_DIR [-o OUTPUT_DIR] [-m MAX_DIST] [-n MIN_DIST] [-w BIN_WIDTH]\n"
        "               [-x|-y] [-i MIN_SEQ_SEPARATION] [-j MAX_SEQ_SEPARATION] [-t ENERGY_THRESHOLD] [-r REPRESENTATION] [-g] [-p]\n"
        "               [-W|-X|-R REF_STATE_FILE] [-Y|-U|-Z REF_STATE_FILES] [-f R_SCRIPTS_DIR] [-k KERNEL] [-b BANDWIDTH] [-c CUT]\n"
        "               [-s ALTERNATIVE_BANDWIDTH2] [-a ADJUST] [-e ALTERNATIVE_BANDWIDTH1]\n"
        "\n"
        "    Options:\n"
        "    -l    string    List of PDB codes (each made of 4 to 5 characters; case sensitive)\n"
        "    -L    string    Alternative to -l: List of coordinate files; This allows training on non-native protein structures\n"
        "    -d    string    Directory containing the coordinate files (format: xxxx.pdb, XXXX.pdb, or pdbxxxx.ent; useless with -L)\n"
        "    -o    string    Output directory; will be created if not existing\n"
        "    -m    float     Maximum interatomic distance (Å) (default=15.0)\n"
        "    -n    float     Minimum interatomic distance (Å) (default=0.0)\n"
        "    -w    float     Width of each bin in the distance distributions (Å) (default=0.1)\n"
        "    -x              Inter-chain interactions only\n"
        "    -y              Intra-chain interactions only\n"
        "    -i    int       Minimum number of positions separating the residue pair (default=3)\n"
        "    -j    int       Maximum number of positions separating the residue pair\n"
        "    -t    float     Threshold for the pseudo-energy (default=10.0)\n"
        "    -r    string    Representation: CA (Calpha), CB (Cbeta), BB (backbone beads), backbone, sidechains, sidechainsCG,\n"
        "                                    SC1 (side chain beads), allatom, or allatomCG (default=CA)\n"
        "    -g              Force write the distance files (*.dat) even when not using KDE\n"
        "    -p              Plot the pseudo-energy profiles\n"
        "    -f    string    Directory containing the R scripts (default=./src/)\n"
        "    -A              For processing large input datasets; not compatible with Kernel Density Estimations\n"
        "\n"
        "    Separate computation of the reference state (use -Y/-U and -Z when training the reference state on decoys):\n"
        "    -W              Write into a file ('frequencies.ref') the frequencies of the reference distribution (see also -R option below)\n"
        "    -X              Same as -W, but for large input datasets (Note: -X or -W option implies an early stoppage of the program)\n"
        "    -R    string    Use a pre-computed reference distribution (e.g. -R /home/me/frequencies.ref)\n"
        "    -Y              Write the frequencies into '*.frq' files for every atomic pair (see also -Z option below)\n"
        "    -U              Same as -Y, but for large input datasets (Note: -U or -Y option implies an early stoppage of the program)\n"
        "    -Z    string    Directory containing the pre-computed frequencies (*.frq) to use as reference frequencies\n"
        "\n"
        "    Related to the density() function from the R standard library (Kernel Density Estimations):\n"
        "    -k    string    Kernel: g (gaussian), e (epanechnikov), r (rectangular), t (triangular),\n"
        "                            b (biweight), c (cosine), or o (optcosine); (default=g)\n"
        "    -b    string    Bandwidth selection method: nrd, nrd0, bcv, ucv, SJ-dpi, or SJ-ste (default=nrd0)\n"
        "    -e    string    Alternative bandwidth in case of error with SJ methods: bcv, ucv, nrd0, or nrd (default=bcv)\n"
        "    -s    string    Alternative bandwidth in case of slowness: nrd0 or nrd (default=nrd0)\n"
        "    -a    float     Adjust (default=1.0)\n"
        "    -c    float     Cut (default=3.0)\n"
        "\n"
        "    -h              Help\n";

    int diffmin = 3;
    int diffmax = 5100; // there is no protein *chain* longer than 5100
    double binwidth = 0.1;
    double emax = 10.0;
    double distmax = 15.0;
    double distmin = 0.0;
    string listpdb;
    string listfiles;
    string inputdir;
    string kernel;
    string bandwidth;
    string jobdir;
    string altbwselector1; // For R errors
    string altbwselector2; // For R slowness
    string cacbbb = "CA";
    string adjust;
    string cut;
    string refdistrib;
    string rscripts = "./src/";
    bool interchain = true;
    bool intrachain = true;
    bool stillwrite = false;
    bool writeRefFreq = false;
    bool xwriteRefFreq = false;
    bool plotprofiles = false;
    bool writeFreq = false;
    bool uwriteFreq = false;
    string dirFreqFiles;
    bool largeinput =false;

    int opt;
    while ((opt = getopt(argc,argv,"hxygWXYUApl:L:d:m:n:o:w:k:b:a:e:s:i:j:t:r:c:f:R:Z:")) != EOF){
        switch(opt){
            case 'l': listpdb = optarg; break;
            case 'L': listfiles = optarg; break;
            case 'd': inputdir = optarg; break;
            case 'm': distmax = atof(optarg); break;
            case 'n': distmin = atof(optarg); break;
            case 'w': binwidth = atof(optarg); break;
            case 'o': jobdir = optarg; break;
            case 'k': kernel = optarg; break;
            case 'b': bandwidth = optarg; break;
            case 'a': adjust = optarg; break;
            case 'e': altbwselector1 = optarg; break;
            case 's': altbwselector2 = optarg; break;
            case 'f': rscripts = optarg; break;
            case 'x': intrachain = false; break;
            case 'y': interchain = false; break;
            case 'g': stillwrite = true; break;
            case 'W': writeRefFreq = true; break;
            case 'X': xwriteRefFreq = true; break;
            case 'R': refdistrib = optarg; break;
            case 'Y': writeFreq = true; break;
            case 'U': uwriteFreq = true; break;
            case 'Z': dirFreqFiles = optarg; break;
            case 'A': largeinput = true; break;
            case 'p': plotprofiles = true; break;
            case 'i': diffmin = atoi(optarg); break;
            case 'j': diffmax = atoi(optarg); break;
            case 't': emax = atof(optarg); break;
            case 'r': cacbbb = optarg; break;
            case 'c': cut = optarg; break;
            case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
        }
    }

    vector<string> possibleb { "nrd", "nrd0", "bcv", "ucv", "SJ-dpi", "SJ-ste" };
    vector<string> possiblek { "g", "e", "r", "t", "b", "c", "o" };
    vector<string> possibler { "CA", "CB", "BB", "SC1", "allatom", "allatomCG", "backbone", "sidechains", "sidechainsCG" };
    vector<string> possiblee { "nrd", "nrd0", "bcv", "ucv" };
    vector<string> possibles { "nrd", "nrd0" };

    if (argc == 1){ fprintf(stderr, "%s", optlist.c_str()); return 1; }
    if (listpdb.empty() and listfiles.empty()){ cerr << "\nError: Missing -l or -L argument\n" << endl; return 1; }
    if (!listpdb.empty() and !listfiles.empty()){ cerr << "\nError: -l and -L are incompatible options\n" << endl; return 1; }
    if (writeRefFreq and !refdistrib.empty()){ cerr << "\nError: -W and -R are incompatible options\n" << endl; return 1; }
    if (xwriteRefFreq and !refdistrib.empty()){ cerr << "\nError: -X and -R are incompatible options\n" << endl; return 1; }
    if (xwriteRefFreq and writeRefFreq){ cerr << "\nError: -X and -W are incompatible options\n" << endl; return 1; }
    if (writeFreq and !dirFreqFiles.empty()){ cerr << "\nError: -Y and -Z are incompatible options\n" << endl; return 1; }
    if (uwriteFreq and !dirFreqFiles.empty()){ cerr << "\nError: -U and -Z are incompatible options\n" << endl; return 1; }
    if (uwriteFreq and writeFreq){ cerr << "\nError: -U and -Y are incompatible options\n" << endl; return 1; }
    if (uwriteFreq and writeRefFreq){ cerr << "\nError: -U and -W are incompatible options\n" << endl; return 1; }
    if (uwriteFreq and xwriteRefFreq){ cerr << "\nError: -U and -X are incompatible options\n" << endl; return 1; }
    if (writeFreq and writeRefFreq){ cerr << "\nError: -Y and -W are incompatible options\n" << endl; return 1; }
    if (writeFreq and xwriteRefFreq){ cerr << "\nError: -Y and -X are incompatible options\n" << endl; return 1; }
    if (uwriteFreq and !refdistrib.empty()){ cerr << "\nError: -U and -R are incompatible options\n" << endl; return 1; }
    if (writeFreq and !refdistrib.empty()){ cerr << "\nError: -Y and -R are incompatible options\n" << endl; return 1; }
    if (writeRefFreq and !dirFreqFiles.empty()){ cerr << "\nError: -W and -Z are incompatible options\n" << endl; return 1; }
    if (xwriteRefFreq and !dirFreqFiles.empty()){ cerr << "\nError: -X and -Z are incompatible options\n" << endl; return 1; }
    if (distmin > distmax){ cerr << "\nError: -m arg must be greater than -n arg\n" << endl; return 1; }
    if (diffmin > diffmax){ cerr << "\nError: -j arg must be greater than -i arg\n" << endl; return 1; }
    if (diffmin < 0){ cerr << "\nError: -i argument cannot be negative\n" << endl; return 1; }
    if (diffmax < 0){ cerr << "\nError: -j argument cannot be negative\n" << endl; return 1; }
    if (distmax < 0){ cerr << "\nError: -m argument cannot be negative\n" << endl; return 1; }
    if (distmin < 0){ cerr << "\nError: -n argument cannot be negative\n" << endl; return 1; }
    if (atof(adjust.c_str()) < 0){ cerr << "\nError: -a argument cannot be negative\n" << endl; return 1; }
    if (binwidth < 0){ cerr << "\nError: -w argument cannot be negative\n" << endl; return 1; }

    if(!bandwidth.empty() and !(find(possibleb.begin(), possibleb.end(), bandwidth) != possibleb.end())) {
        cerr << "\nError: Invalid -b argument\n" << endl; return 1;
    }
    if(!kernel.empty() and !(find(possiblek.begin(), possiblek.end(), kernel) != possiblek.end())) {
        cerr << "\nError: Invalid -k argument\n" << endl; return 1;
    }
    if(!(find(possibler.begin(), possibler.end(), cacbbb) != possibler.end())) {
        cerr << "\nError: Invalid -r argument\n" << endl; return 1;
    }
    if(!altbwselector1.empty() and !(find(possiblee.begin(), possiblee.end(), altbwselector1) != possiblee.end())) {
        cerr << "\nError: Invalid -e argument\n" << endl; return 1;
    }
    if(!altbwselector2.empty() and !(find(possibles.begin(), possibles.end(), altbwselector2) != possibles.end())) {
        cerr << "\nError: Invalid -s argument\n" << endl; return 1;
    }
    if (inputdir.empty() and listfiles.empty()){ cerr << "\nError: Missing -d argument\n" << endl; return 1; }
    if (!inputdir.empty() and !listfiles.empty()){ cerr << "\nWarning: -d argument useless with -L argument\n" << endl; }
    if (stillwrite and xwriteRefFreq){ cerr << "\nWarning: -g argument useless with -X argument\n" << endl; }

    bool allatom = (cacbbb == "allatom" or cacbbb == "backbone" or cacbbb == "sidechains") ? true : false;

    // Create output directory
    bool defaultdir = false;
    if (jobdir.empty()){
        defaultdir = true;
        string timer = to_string(time(NULL));
        string xopt = (!intrachain) ? "_x" : "";
        string yopt = (!interchain) ? "_y" : "";
        stringstream stream;
        stream << "k" << kernel
               << "_b" << bandwidth
               << "_e" << altbwselector1
               << "_s" << altbwselector2
               << "_c" << cut
               << "_a" << adjust
               << "_m"  << fixed << setprecision(2) << distmax
               << "n" << fixed << setprecision(2) << distmin
               << "w" << fixed << setprecision(2) << binwidth
               << "i" << fixed << setprecision(0) << diffmin
               << "j" << fixed << setprecision(0) << diffmax
               << "t" << fixed << setprecision(2) << emax
               << xopt << yopt
               << "_r" << cacbbb;
        jobdir = stream.str();
    }

    struct stat sb;
    if (stat(jobdir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
        if (defaultdir){ // No warning for user-defined names
            cout << "\nWarning: Directory \"" << jobdir << "\" already exists" << endl;
            cout << "Press ENTER to continue, or CTRL+C to cancel..." << endl;
            getchar();
        }
    }
    else{
        cout << "\nCreating output directory: \"" << jobdir << "\" ..." << endl;
        mkdir(jobdir.c_str(), ACCESSPERMS);
        cout << "Done" << endl;
    }

    if (!intrachain and !interchain){
        cerr << "\nError: -x and -y are incompatible options" << endl;
        return 1;
    }

    bool bwSJ = (bandwidth == "SJ-ste" or bandwidth == "SJ-dpi") ? true : false;
    bool bwCV = (bandwidth == "bcv" or bandwidth == "ucv") ? true : false;
    bool useR = false;

    if (!kernel.empty() or !bandwidth.empty() or !adjust.empty() or !cut.empty()){
        if (largeinput){cerr << "\nError: -A option is not compatible with KDE\n"; return 1;}
        cout << "\nKernel Density Estimations activated" << endl;
        kernel = (!kernel.empty()) ? " "+kernel : " g";
        bandwidth = (!bandwidth.empty()) ? " "+bandwidth : " nrd0";
        adjust = (!adjust.empty()) ? " "+adjust : " 1";
        cut = (!cut.empty()) ? " "+cut : " 3"; // 3 is the default value in R
        useR = true;
        if (xwriteRefFreq or uwriteFreq)
            cerr << "Warning: KDE are not applied when -X or U option is activated" << endl;
        if (stillwrite){cerr << "\nError: No need for -g option when KDE activated:\n"
                             << "The interatomic distance files are written anyway" << endl; return 1;}
    }

    if (useR){
        if(!fexists(rscripts+"/plotpot.R")){cerr << "\nError: R script \"plotpot.R\" not found in \"" << rscripts << "\"" << endl; return 1;}
        if(!fexists(rscripts+"/kde.R")){cerr << "\nError: R script \"kde.R\" not found in \"" << rscripts << "\"" << endl; return 1;}
    }

    if (!altbwselector1.empty() and !bwSJ){
        cerr << "Error: -e option useful only for SJ bandwidth selectors" << endl; exit(1);}

    if (!altbwselector2.empty() and !bwSJ and !bwCV){
        cerr << "Error: -s option useful only for SJ, BCV, and UCV bandwidth selectors" << endl; exit(1);}

    if (!altbwselector1.empty())
        altbwselector1 = " "+altbwselector1;
    else if (bwSJ)
        altbwselector1 = " bcv";

    if (!altbwselector2.empty())
        altbwselector2 = " "+altbwselector2;
    else if (bwSJ or bwCV)
        altbwselector2 = " nrd0";

    int total = 0; // (x²+x)/2
    if (cacbbb == "allatom"){total = 14028;}
    else if (cacbbb == "allatomCG"){total = 1225;}
    else if (cacbbb == "sidechains"){total = 3916;}
    else if (cacbbb == "sidechainsCG"){total = 496;}
    else if (cacbbb == "backbone"){total = 3240;}
    else {total = 210;} // CA, CB, BB, SC1

    string dontplot = (allatom) ? " 1" : ""; // Prevent the creation of 14028 (allatom), or 3240 (backbone), or 3916 (sidechains) plots

    /*****************************
    *                            *
    *   END OF OPTION HANDLING   *
    *                            *
    ******************************/



    map<string,string> three2one = create_map();
    set<string> atypes;
    if (cacbbb == "allatom") atypes = set_allatom();
    else if (cacbbb == "allatomCG") atypes = set_allatomCG();
    else if (cacbbb == "sidechains") atypes = set_sidechains();
    else if (cacbbb == "sidechainsCG") atypes = set_sidechainsCG();
    else if (cacbbb == "backbone") atypes = set_backbone();
    else if (cacbbb == "CA") atypes = set_CA();
    else if (cacbbb == "CB") atypes = set_CB();
    else if (cacbbb == "BB") atypes = set_BB();
    else atypes = set_SC1();

    // Listing all atom pairs
    vector<string> atomcombos;
    for (auto it1 : atypes)
        for (auto it2 : atypes)
            if (it1.compare(it2) <= 0)
                atomcombos.push_back(it1+it2);

    // Open the list of PDB files
    vector<string> inputpdbs; // List of 4-character PDB codes (if -l), or a list of filenames (if -L)
    vector<string> inputchains; // List of chains

    if (!listpdb.empty()){ // If -l argument
        ifstream fh(listpdb);
        if(!fh){ cerr << "Error while opening input list \"" << listpdb << "\"" << endl; exit(1); }
        string inputcode;
        while(getline(fh, inputcode)){
            if (inputcode.size() > 3){
                inputpdbs.push_back(inputcode.substr(0,4));
                if (inputcode.size() > 4)
                    inputchains.push_back(inputcode.substr(4,1));
                else
                    inputchains.push_back("NA");
            }
            else{
                cerr << "Error: " << inputcode << " is not a proper PDB code" << endl;
                exit(1);
            }
        }
        fh.close();
    }
    else{ // If -L argument
        ifstream fh(listfiles);
        if(!fh){ cerr << "Error while opening input list \"" << listfiles << "\"" << endl; exit(1); }
        string inputcode;
        while(getline(fh, inputcode)){
            if (find(inputcode.begin(), inputcode.end(), ' ') != inputcode.end()){ // If chain name provided
                vector<string> vec = split(inputcode, ' ');
                inputpdbs.push_back(vec[0]);
                if (vec[1].size() == 1 and isalnum(vec[1][0])) // If chain name is valid
                    inputchains.push_back(vec[1]);
                else
                    inputchains.push_back("NA");
            }
            else{
                inputpdbs.push_back(inputcode);
                inputchains.push_back("NA");
            }
        }
        fh.close();
    }

    map<string, vector<float>> allsqdist; // Key: atom pair (e.g. CAICBL); value: vector of distances
    vector<pair<string,double>> stat1; // string: pairID; double: distance
    cout << "\nProcessing the input PDB files:" << endl;
    int monochains = 0; int multichains = 0;

    ofstream tfh; // To write either frequencies.ref or xx.dat



    // The following variable is used only when -X option is activated
    map<double, unsigned long int> xxcounts; // For large datasets (-X option)

    // Initializing xxcounts map: each bin is defined by its upper limit...
    // ...while defining the bin limits into a vector
    vector<double> binlimits;
    for (double i=distmin+binwidth; i<=distmax; i+=binwidth){
        // for loops using double (or float) variables are dangerous:
        stringstream istream;
        istream << fixed << setprecision(3) << i;
        istream >> i;

        binlimits.push_back(i);
        xxcounts[i] = 0;
    }
    int binlimitssize = binlimits.size();

    unsigned long int xxtotal = 0; // For large datasets (-X option; also used for the total number of distances)

    // Used only when -U option is activated
    map<string, map<double, unsigned long int>> nncounts; // key = atom pair; value = a map similar to xxcounts
    map<string, unsigned long int> nntotal;

    // Initializing nncounts
    for (auto atomcombo : atomcombos)
        for (auto binlimit : binlimits)
            nncounts[atomcombo][binlimit] = 0;

 
    // For each PDB file of the list
    for (size_t k=0; k<inputpdbs.size(); ++k){
        string filename = listpdb.empty() ? inputpdbs[k] : findfile(inputpdbs[k], inputdir);

        vector<string> filevec = file2vec(filename, three2one, atypes);

        vector<string>& reffilevec(filevec);
        // Declare three vectors
        vector<float> sqdist; vector<string> atompairs; vector<string> pairID;
        // Declare their three references
        vector<float>& refsqdist(sqdist); vector<string>& refatompairs(atompairs); vector<string>& refpairID(pairID);
        // For the current PDB file: find the chains
        vector<string> chains = findchains(reffilevec);

        // If chain name in input
        if (inputchains[k] != "NA"){
            // And if it actually exists
            if (find(chains.begin(), chains.end(), inputchains[k]) != chains.end()){
                chains.clear();
                chains.push_back(inputchains[k]);
            }
            else{
                cerr << " Warning: There is no chain " << inputchains[k] << " in " << filename << endl;
                cerr << "Will process every chain found in " << filename << endl;
            }
        }

        if (chains.size() == 1)
            monochains++;
        else if (chains.size() > 1)
            multichains++;
        else
            cerr << "Warning: No chain found in " << filename << endl;

        // For the current PDB file:
        // For each combination of chain:
        for (size_t i=0; i<chains.size(); ++i){
            for (size_t j=i; j<chains.size(); ++j){
                // Fill the three vectors (their references: refsqdist, refatompairs, refpairID)
                if (interchain)
                    if (chains[i] != chains[j])
                        squaredist(inputpdbs[k], reffilevec, refsqdist, refatompairs, refpairID,
                                     distmax, distmin, chains[i], chains[j], diffmin, diffmax, three2one, allatom);
                if (intrachain) // NOT ELSE
                    if (chains[i] == chains[j])
                        squaredist(inputpdbs[k], reffilevec, refsqdist, refatompairs, refpairID,
                                     distmax, distmin, chains[i], chains[j], diffmin, diffmax, three2one, allatom);
            }
        }

        // If -X, -U, or -A options activated: on-the-fly computing of the histogram
        if (xwriteRefFreq or uwriteFreq or largeinput){
            sort(sqdist.begin(), sqdist.end());
            int binindex = 0;
            for(size_t j=0; j<sqdist.size(); ++j){
                while(binlimits[binindex] < sqdist[j]){
                    binindex++;
                }

                nncounts[atompairs[j]][binlimits[binindex]]++; // -U option
                nntotal[atompairs[j]]++; // -U option

                xxcounts[binlimits[binindex]]++; // -X option
                xxtotal++; // -X option
            }
        }
        // Else: fill the map (Note: sqdist, atompairs, and pairID have the same size)
        else{
            for(size_t i=0; i<sqdist.size(); ++i){
                allsqdist[atompairs[i]].push_back(sqdist[i]);
                if(!allatom) stat1.emplace_back(pairID[i],sqdist[i]);
                xxtotal++;
            }
        }
        cout.flush();
        cout << "\r" << k+1 << "/" << inputpdbs.size();
    }
    cout << "\n" << xxtotal << " distances have been computed" << endl;
    cout << "Done\n" << endl;
    /* All the coordinate files have been processed! */



    cout << "*Structures detected as monomeric:  n = " << monochains << endl;
    cout << "*Structures detected as multimeric: n = " << multichains << "\n" << endl;
    if (monochains == 0 and multichains == 0){cerr << "Error: No chain found. Program stopped." << endl; exit(1);}



    // If -X option activated: now compute the frequencies and write them into a file
    if (xwriteRefFreq){
        cout << "Writting the frequencies of the reference distribution (filename: frequencies.ref)..." << endl;
        tfh.open(jobdir+"/frequencies.ref");
        for(auto& it : xxcounts){
            tfh << double(it.second)/double(xxtotal) << endl;
        }
        tfh.close();
        cout << "Done" << endl;

        return 0;
    }

    if (uwriteFreq){
        cout << "Writing the frequencies for every atomic pair: " << total << " *.frq files..." << endl;
        for(auto& it1 : nncounts){
            tfh.open(jobdir+"/"+it1.first+".frq");
            for(auto& it2 : it1.second){
                tfh << double(it2.second)/double(nntotal[it1.first]) << endl;
            }
            tfh.close();
        }
        cout << "Done" << endl;

        return 0;
    }



    // Listing atom pairs that have not been found
    if (!largeinput){ // Requires -A option not activated (otherwise, allsqdist is empty)
        vector<string> missingpairs;
        for (auto atomcombo : atomcombos){
            int found = 0;
            for(auto& it : allsqdist)
                if (it.first == atomcombo)
                    found++;
            if (found == 0)
                missingpairs.push_back(atomcombo);
        }
    
        if (missingpairs.size() > 0){
            ofstream missfh;
            missfh.open(jobdir+"/missing_pairs.dat");
            for (auto missingpair : missingpairs)
                missfh << missingpair << endl;
            missfh.close();
            cerr << "Warning: " << missingpairs.size() << " atom pairs (out of " << total << ") have not been found." << endl;
            cerr << "    They have been listed in " << jobdir << "/missing_pairs.dat\n" << endl;
        }
    }



    bool writedat = (stillwrite or useR) ? true : false;

    if(writedat) cout << "Writing the interatomic distance files..." << endl;
    vector<pair<string,double>> stat3; // string: atom pair; double: counts

    vector<int> totalcounts(binlimitssize);
    int totalsize = 0;

    map<string, vector<double>> frequencies; // key = atom pair; value = vector of frequencies per distance bin

    int counter = 0;
    cout << counter << "/" << total;
    if(writedat) tfh.open (jobdir+"/xx.dat");


    if (largeinput){
        for (auto atomcombo : atomcombos)
            for (auto binlimit : binlimits)
                frequencies[atomcombo].push_back(double(nncounts[atomcombo][binlimit])/double(nntotal[atomcombo]));
 
        if(refdistrib.empty())
            for (auto binlimit : binlimits)
                frequencies["xx"].push_back(double(xxcounts[binlimit])/double(xxtotal));
    }
    else{
        // For each key of the map (= atom pair): create a file that contains distances (e.g. "ICALCA.dat", "LCAVCA.dat" etc.)
        for(auto& it : allsqdist){
            if(writedat){
                ofstream ofh;
                if (!writeRefFreq)
                    ofh.open (jobdir+"/"+it.first+".dat");
                // Write all the distances in the file
                for(auto distance : it.second){
                    if (!writeRefFreq)
                        ofh << distance << endl;
                    if (!writeFreq)
                        tfh << distance << endl; // and in the total file xx.dat
                }
                ofh.close();
            }
            stat3.emplace_back(it.first,(double)it.second.size()); // Cast to double because of the sortbysec function
    
            if (!useR){
                sort(it.second.begin(), it.second.end());
                int binindex = 0;
                vector<int> counts(binlimitssize);
                for (double distance : it.second){
     
                    while(binlimits[binindex] < distance){
                        binindex++;
                    }
                    counts[binindex]++;
                    totalcounts[binindex]++;
                }
                for (int count : counts){
                    frequencies[it.first].push_back((double)count/(double)it.second.size());
                }
                totalsize+=it.second.size();
            }
    
            counter++;
            cout.flush();
            cout << "\r" << counter << "/" << total;
        }
    }
    if(writedat) tfh.close();
    string note1 = (writeRefFreq and writedat) ? " (-W options activated: Only xx.dat has been written)" : "";
    cout << "\nDone" << note1 << "\n" << endl;

    if (!useR){
        if (refdistrib.empty()){
            for (int totalcount : totalcounts)
                frequencies["xx"].push_back((double)totalcount/(double)totalsize);
        }
        if (writeRefFreq){
            map<string, vector<double>>& reffrequencies(frequencies);
            f_writeRefFreq(reffrequencies, jobdir);
            return 0;
        }
    }
    else{ // if (useR)
        if (writeRefFreq)
            allsqdist.clear();

        if (refdistrib.empty() and !writeFreq) // if no -R and no -Y
            allsqdist["xx"]; // Add the "total atom pair"

        counter = 0;
        string plot = (dontplot.empty()) ? " and plotting" : "";
        cout << "Computing" << plot << " the distance frequencies for every pair of atoms:" << endl;
        cout << counter << "/" << total;
    
        // For each atom pair: call R script to compute distance frequencies
        for(auto& it : allsqdist){
            string cmd;
            string cmdout;
            string infilename = jobdir+"/"+it.first+".dat";
 
            // If input size too large:
            // Use another bandwidth selector instead
            string sizecmd = "du -k "+infilename+" | cut -f1";
            string sizeout = exec(sizecmd);
            sizeout.erase(remove(sizeout.begin(), sizeout.end(), '\n'), sizeout.end());
            int maxsize = (it.first == "xx") ? MAX_SIZE_KDEXX : MAX_SIZE_KDE; // One can wait longer if it's only for the xx pair!
 
            if (atoi(sizeout.c_str()) > maxsize and (bwSJ or bwCV)){
                cmd = "Rscript "+rscripts+"/kde.R "+infilename+" "+to_string(binwidth)+" "+to_string(distmax)+kernel+altbwselector2+adjust+cut+dontplot+" 2> /dev/null";
                cmdout = exec(cmd);
                cerr << " Warning: Input too large:" << altbwselector2 << " bandwidth selector used for " << it.first << " atom pair" << endl;
            }
            else{
                cmd = "Rscript "+rscripts+"/kde.R "+infilename+" "+to_string(binwidth)+" "+to_string(distmax)+kernel+bandwidth+adjust+cut+dontplot+" 2> /dev/null";
                cmdout = exec(cmd);
            }
 
            // If the R implementation of SJ bandwidth selection fails (may happen for large samples):
            // Use another bandwidth selector instead
            if (cmdout.empty()){
                if (bwSJ){
                    cmd = "Rscript "+rscripts+"/kde.R "+infilename+" "+to_string(binwidth)+" "+to_string(distmax)+kernel+altbwselector1+adjust+cut+dontplot+" 2> /dev/null";
                    cmdout = exec(cmd);
                    if (cmdout.find("breaks") && cmdout.find("hist")){
                        cerr << "\nError: The distance range (defined by -m and -n args) is too short for the chosen bin width (-w arg)" << endl;
                        cout << "Program stopped." << endl; exit(1);
                    }
                    else{
                        cerr << " Warning: Problem with R:" << altbwselector1 << " bandwidth selector used for " << it.first << " atom pair" << endl;
                    }
                }
                else{
                    cerr << "Unexpected error with atom pair " << it.first << endl;
                    cout << "Program stopped." << endl; exit(1);
                }
            }

            stringstream iss(cmdout);
            string freq;
            while (iss >> freq){
                if (freq == "NA"){
                    freq = "0";
                }
                frequencies[it.first].push_back(atof(freq.c_str()));
            }
            counter++;
            cout.flush();
            if (it.first != "xx")
                cout << "\r" << counter << "/" << total;
            else
                // The "allsqdist" map is ordered, thus the xx pair is last
                cout << "\nDone\n\nComputing the " << total+1 << "th pair: two atoms of any type...";
        }
        cout << "\nDone\n" << endl;

        if (writeRefFreq){
            map<string, vector<double>>& reffrequencies(frequencies);
            f_writeRefFreq(reffrequencies, jobdir);
            return 0;
        }
    } // END OF if useR

    // If using a pre-computed reference distribution
    if (!refdistrib.empty()){
        cout << "Using reference distribution from input file named \"" << refdistrib << "\"" << endl;
        ifstream fh(refdistrib);
        if(!fh){ cerr << "Error while opening reference distribution file \"" << refdistrib << "\"" << endl; exit(1); }
        string freqvalue;
        while(getline(fh, freqvalue)){
            frequencies["xx"].push_back(atof(freqvalue.c_str()));
        }
        fh.close();

        size_t providedsize = frequencies["xx"].size();
        size_t expectedsize = (frequencies.begin()->second).size();

        // Checking whether the input file has the correct number of bins
        if (providedsize != expectedsize){
            cerr << "Error: The number of bins in your reference distribution is not compatible with the other parameters that you have selected." << endl;
            cerr << "       Number of bins (= lines) in your reference file = " << providedsize << endl;
            cerr << "       Number of bins expected = " << expectedsize << endl;
            cerr << "Program stopped." << endl;
            exit(1);
        }
    }



    if (writeFreq){
        cout << "Writing the frequencies for every atomic pair: " << total << " *.frq files..." << endl;
        for(auto& it : frequencies){
            ofstream freqfh;
            freqfh.open(jobdir+"/"+it.first+".frq");
            for (size_t i=0; i<it.second.size(); ++i){
                freqfh << it.second[i] << endl;
            }
            freqfh.close();
        }
        cout << "Done" << endl;
        return 0;
    }



    vector<pair<string,double>> stat2; // string: atom pair; double: energy
    cout << "Computing the pseudo-energies..." << endl;

    // For each atom pair
    bool useFreqFiles = (dirFreqFiles.empty()) ? false : true; // -Z option
    for(auto& it : frequencies){
        if (it.first != "xx"){
            vector<double> reffreq;
            ifstream freqfh;

            if (useFreqFiles){
                freqfh.open(dirFreqFiles+"/"+it.first+".frq");
                string frequency;
                while(getline(freqfh, frequency)){
                    reffreq.push_back(stod(frequency));
                }
            }

            ofstream ofh;
            ofh.open (jobdir+"/"+it.first+".nrg");
            vector<double> energies;
            // For each distance bin
            for (size_t i=0; i<it.second.size(); ++i){
                double energy = 0;
                double obs = it.second[i]; // Observed frequency
                double ref = (useFreqFiles) ? reffreq[i] : frequencies["xx"][i]; // Reference frequency
                if (ref > 0){
                    energy = -1*log(obs/ref);
                }
                if (std::isinf(energy) or energy == 0){ // Leave std::
                    energy = emax;
                }
                energies.push_back(energy);
            }

            // Correction procedure by interpolation
            for (size_t i=0; i<energies.size(); ++i){
                if (energies[i] == emax and energies[i-1] != emax and i>0 and i<(energies.size()-1)){
                    if (energies[i+1] != emax)
                        energies[i] = (energies[i-1]+energies[i+1])/2;
                    else
                        energies[i-1] = emax;
                }
            }

            // Correction related to the use of the 'cut' parameter in R density()
            for (size_t i=0; i<energies.size(); ++i){
                if (energies[i] == emax and energies[i-1] != emax and i>0)
                    energies[i] = 0;
            }

            // Write the pseudo-energies
            for (auto energy : energies){
                ofh << energy << endl;
                // Fill the map that will be used to rank the lowest pseudo-energies
                stat2.emplace_back(it.first,energy);
            }
            ofh.close();
        }
    }

    cout << "Done\n" << endl;

    // Write the bin centers into a file
    ofstream myfh;
    myfh.open (jobdir+"/xvector.dat");
    for (double i=distmin+(binwidth/2); i<=distmax-(binwidth/2); i+=binwidth)
        myfh << i << endl;
    myfh.close();

    // Record the parameters
    myfh.open (jobdir+"/parameters.log");
    myfh << "DIFFMIN=" << diffmin << endl;
    myfh << "DIFFMAX=" << diffmax << endl;
    myfh << "DISTMAX=" << distmax << endl;
    myfh << "DISTMIN=" << distmin << endl;
    myfh << "REPRES=" << cacbbb << endl;
    myfh.close();

    // Sorting the vectors
    sort(stat1.begin(), stat1.end(), sortbysec); // distances
    sort(stat2.begin(), stat2.end(), sortbysec); // energies
    sort(stat3.begin(), stat3.end(), sortbysec); // occurrences

    // Compute some stats
    cout << "Writing files";

    int count = 0;
    if (!allatom and !largeinput){
        myfh.open (jobdir+"/"+"top_distances.tsv");
        cout << " \"top_distances.tsv\",";
        while (count<100){
            myfh << stat1[count].first << "\t" << stat1[count].second << endl;
            count++;
        }
        myfh.close();
    }

    myfh.open (jobdir+"/"+"top_energies.tsv");
    cout << " \"top_energies.tsv\",";
    set<string> printed;
    for (size_t i=0; i<stat2.size(); ++i){
        // if atom pair has not been printed
        if (printed.find(stat2[i].first) == printed.end()){
            myfh << stat2[i].first << "\t" << stat2[i].second << endl;
            printed.insert(stat2[i].first);
        }
    }
    myfh.close();

    if (!largeinput){
        myfh.open (jobdir+"/"+"top_occurrences.tsv");
        count = total-1;
        cout << " and \"top_occurrences.tsv\"";
        while (count>=0){
            myfh << stat3[count].first << "\t" << stat3[count].second << endl;
            count--;
        }
        myfh.close();
    }

    cout << "\nDone\n" << endl;

    if (plotprofiles){
        // For each atom pair: plot the pseudo-energy profile (SVG file)
        cout << "Plotting the pseudo-energy profiles..." << endl;
        count = 0;
        cout << count << "/" << total;
        for(auto& it : frequencies){
            if (it.first != "xx"){
                string plotcmd = "Rscript "+rscripts+"/plotpot.R "+jobdir+"/"+it.first+".nrg "+to_string(binwidth)+" "+to_string(distmax)+" "+it.first+" "+jobdir+"/"+" > /dev/null";
                exec(plotcmd);
                count++;
                cout.flush();
                cout << "\r" << count << "/" << total;
            }
        }
        cout << "\nDone\n" << endl;
    }

    return 0;
}


/************
* FUNCTIONS *
************/

void f_writeRefFreq(map<string, vector<double>>& frequencies, string jobdir){
    ofstream myfh;
    cout << "Writting the frequencies of the reference distribution (filename: frequencies.ref)..." << endl;
    myfh.open(jobdir+"/"+"frequencies.ref");
    for (size_t i=0; i<frequencies["xx"].size(); ++i){
        myfh << frequencies["xx"][i] << endl;
    }
    myfh.close();
    cout << "Done\n" << endl;
    return;
}


void squaredist(string pdbcode, vector<string>& filevec, vector<float>& sqdist, vector<string>& atompairs, vector<string>& pairID,
                  double distmax, double distmin, string chain1, string chain2, int diffmin, int diffmax, map<string,string> three2one, bool allatom){
/* This functions
   - reads the vector containing (part of) the PDB file
   - calculates all the interatomic distances for the one or two chains
   - push these distances into a referenced vector
   - there are 2 other referenced vectors: 1 for the atoms names, 1 for the atoms names+numbers
*/

    vector<int> num1; vector<int> num2;
    vector<double> x1; vector<double> x2;
    vector<double> y1; vector<double> y2;
    vector<double> z1; vector<double> z2;
    vector<string>atom1; vector<string>atom2;

    for (auto line : filevec){
        string oneletter = three2one[line.substr(17,3)];
        string chainfound = line.substr(21,1);

        string atomname = line.substr(12,4);
        atomname.erase(remove(atomname.begin(),atomname.end(),' '),atomname.end());
        atomname = oneletter+atomname;

        string alternate = line.substr(16,1);
        int number = atoi(line.substr(22,4).c_str());
        double xfound = atof(line.substr(30,8).c_str());
        double yfound = atof(line.substr(38,8).c_str());
        double zfound = atof(line.substr(46,8).c_str());

        if (alternate == " " or alternate == "A"){ // to avoid residues that have two versions (e.g. ALYS and BLYS)
            if (chainfound == chain1){
                atom1.push_back(atomname);
                num1.push_back(number);
                x1.push_back(xfound);
                y1.push_back(yfound);
                z1.push_back(zfound);
            }
            if (chainfound == chain2){ // NOT else if
                atom2.push_back(atomname);
                num2.push_back(number);
                x2.push_back(xfound);
                y2.push_back(yfound);
                z2.push_back(zfound);
            }
        }
    }

    bool samechain = (chain1 == chain2) ? true : false;
    float sqd; // float precision for distances to save memory: important for large (N>1000) *all-atom* data sets
    const int atom1size = atom1.size(); const int atom2size = atom2.size();
    for(int i=0; i<atom1size; ++i){
        int start = (samechain) ? i+1 : 0;
        for(int j=start; j<atom2size; ++j){
            int numdiff = num2[j]-num1[i];
            if ( (samechain and numdiff > diffmin and numdiff < diffmax) or (!samechain) ){
                sqd = sqrt((x1[i]-x2[j])*(x1[i]-x2[j]) + (y1[i]-y2[j])*(y1[i]-y2[j]) + (z1[i]-z2[j])*(z1[i]-z2[j]));
     
                if (sqd < distmax and sqd > distmin){
                    sqdist.push_back(sqd);
                    if (atom1[i].compare(atom2[j]) < 0){ // if alphabet order
                        atompairs.push_back(atom1[i]+atom2[j]);
                        if(!allatom) pairID.push_back(pdbcode+"\t"+atom1[i]+to_string(num1[i])+chain1+"\t"+atom2[j]+to_string(num2[j])+chain2);
                    }
                    else{
                        atompairs.push_back(atom2[j]+atom1[i]);
                        if(!allatom) pairID.push_back(pdbcode+"\t"+atom2[j]+to_string(num2[j])+chain2+"\t"+atom1[i]+to_string(num1[i])+chain1);
                    }
                }
            }
        }
    }
    return;
}


// Sort the vector elements by second element of pairs
bool sortbysec(const pair<string,double> &a, const pair<string,double> &b){
    return (a.second < b.second);
}


// The functions below are aimed at handling the name formats of PDB files
string makeupper(string s){
    for (size_t i=0; i<s.length(); ++i){
        s[i] = toupper(s[i]);
    }
    return s;
}


string makelower(string s){
    for (size_t i=0; i<s.length(); ++i){
        s[i] = tolower(s[i]);
    }
    return s;
}


string findfile(string code, string inputdir){
    string filename1 = inputdir+"/"+makeupper(code)+".pdb";
    string filename2 = inputdir+"/"+makelower(code)+".pdb";
    string filename3 = inputdir+"/pdb"+makelower(code)+".ent";
    // This 4th format is tested, even though it is not supposed to exist
    string filename4 = inputdir+"/pdb"+makeupper(code)+".ent";

    if(fexists(filename1)) return filename1;
    else if(fexists(filename2)) return filename2;
    else if(fexists(filename3)) return filename3;
    else if(fexists(filename4)) return filename4;
    else{cerr << "\nError: There is no file for PDB code \"" << code << "\" in directory \"" << inputdir << "\"" << endl;
         exit(1);
    }
}
