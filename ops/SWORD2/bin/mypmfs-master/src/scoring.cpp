/*
 * scoring.cpp
 *
 * The scoring part of the MyPMFs suite
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
#include "spline.h"
#include <ctime>
#include <numeric>
#include <sys/stat.h>
#include "functions.h"

using namespace std;

void createdir(string jobdir);
double linear_interpol(double x, string atompair, vector<double>& xvector, map<string, map<double, double>>& nestedmap, int lastbinindex);
void sequence(vector<string>& filevec, vector<string>& resID, vector<string>& aavec, map<string,string> three2one);
void distances(vector<string>& filevec, vector<double>& sqdist, vector<string>& atompairs, vector<string>& resIDl, vector<string>& resIDr,
               double distmax, double distmin, string chain1, string chain2, int diffmin, int diffmax, map<string,string> three2one);

string dirnaming(string filename, string indir, string chainname);
vector<string> residuefilter(vector<string> filevec, string residueslist);


int main(int argc, char** argv){
    srand (time(NULL)); // Initialize random seed

    string optlist =
        "    Usage:\n"
        "    ./scoring {-i INPUT_STRUCTURE | -l INPUT_LIST} -d ENERGY_FILES_DIR [-a CHAIN] [-m MAX_DIST] [-n MIN_DIST] [-x|-y]\n"
        "              [-k MIN_SEQ_SEPARATION] [-j MAX_SEQ_SEPARATION] [-q RESIDUE_LIST] [-r REPRESENTATION] [-c] [-p] [-w]\n"
        "              [-b WINDOW_SIZE] [-o OUTPUT_DIR] [-f R_SCRIPTS_DIR] [-z] [-M] [-s N_DECOYS] [-t SEQ_ID]\n"
        "\n"
        "    Options:\n"
        "    -i    string    Input structure file\n"
        "    -a    string    Chain name; e.g. AbC will process chains A, b and C (default: all chains)\n"
        "    -l    string    List of input structure files (including path);\n"
        "                    (faculative) Space-separated second column: chain names\n"
        "    -d    string    Directory containing the pseudo-energy files (*.nrg)\n"
        "    -o    string    Output directory; will be created if not existing\n"
        "    -m    float     Maximum interatomic distance (Å) (default: parameters.log)\n"
        "    -n    float     Minimum interatomic distance (Å) (default: parameters.log)\n"
        "    -x              Inter-chain interactions only\n"
        "    -y              Intra-chain interactions only\n"
        "    -j    int       Maximum number of positions separating the atom pair (default: parameters.log)\n"
        "    -k    int       Minimum number of positions separating the atom pair (default: parameters.log)\n"
        "    -q    string    Comma-separated list of residues (number+chain) that will be processed (default: all residues)\n"
        "    -r    string    Representation: CA (Calpha), CB (Cbeta), BB (backbone beads), backbone, sidechains, sidechainsCG,\n"
        "                                    SC1 (side chain beads), allatom, or allatomCG (default: parameters.log)\n"
        "    -c              Cubic spline interpolation (by default: linear interpolation)\n"
        "    -p              Plot the interpolated pseudo-energy profiles\n"
        "    -w              Write files \"energy_[window size].tsv\" (energy /position) and \"data.tsv\" (energy and distance /atom pair)\n"
        "    -b    int       Window size for calculating the energy per position (default=1)\n"
        "    -f    string    Directory containing the R scripts (default=./src/)\n"
        "    -z              Z-score computation (only for Calpha and backbone beads representations)\n"
        "    -M              (Z-score related) Mute the counter of random sequence decoys (important when redirecting the output)\n"
        "    -s    int       (Z-score related) Number of random sequence decoys for the Z-score computation (default=1000)\n"
        "    -t    float     (Z-score related) For each random sequence decoy: max seq id with the query structure (default=0.5)\n"
        "    -h              Help\n";

    // The 5 variables below have their default value defined by parameters.log
    int diffmin = -1;
    int diffmax = -1;
    string cacbbb;;
    double distmax = -1;
    double distmin = -1;

    bool interchain = true;
    bool intrachain = true;
    string inputfile;
    string chainname;
    string indir;
    string jobdir;
    bool cubicspline = false;
    double simthres = 0.5;
    int decoys = 1000;
    bool zopt = false;
    bool plotinterpol = false;
    bool write = false;
    string rscripts = "./src/";
    string inputlist;
    string residueslist;
    int windowsize = 12345;
    bool mute = false;

    int opt;
    while ((opt = getopt(argc,argv,"hxyzcpwMm:n:d:i:l:a:j:k:r:t:s:o:f:q:b:")) != EOF){
        switch(opt){
            case 'i': inputfile = optarg; break;
            case 'l': inputlist = optarg; break;
            case 'o': jobdir = optarg; break;
            case 'a': chainname = optarg; break;
            case 'm': distmax = atof(optarg); break;
            case 'n': distmin = atof(optarg); break;
            case 'd': indir = optarg; break;
            case 'x': intrachain = false; break;
            case 'y': interchain = false; break;
            case 'z': zopt = true; break;
            case 'c': cubicspline = true; break;
            case 'p': plotinterpol = true; break;
            case 'w': write = true; break;
            case 'M': mute = true; break;
            case 'b': windowsize = atoi(optarg); break;
            case 'j': diffmax = atoi(optarg); break;
            case 'k': diffmin = atoi(optarg); break;
            case 'r': cacbbb = optarg; break;
            case 'q': residueslist = optarg; break;
            case 'f': rscripts = optarg; break;
            case 't': simthres = atof(optarg); break;
            case 's': decoys = atoi(optarg); break;
            case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
        }
    }

    if (argc == 1){ fprintf(stderr, "%s", optlist.c_str()); return 1; }

    // Use recorded parameters from "parameters.log"
    if (indir.empty()){ cerr << "\nError: -d argument is missing\n" << endl; return 1; }
    string paramfile = indir+"/parameters.log";
    ifstream pfh(paramfile);
    if(!pfh){
        cerr << "\nWarning: File \"" << paramfile << "\" not found" << endl;
        cerr << "Cannot define default values for options -m, -n, -j, -k, and -r" << endl;
    }
    else{
        cout << "\nUsing parameters from " << paramfile << endl;
        string line;
        int linecount = 0;
        while(getline(pfh, line)){
    	    vector<string> vec = split(line, '=');
            if (linecount==0 and diffmin==-1) diffmin=atoi(vec[1].c_str());
            if (linecount==1 and diffmax==-1) diffmax=atoi(vec[1].c_str());
            if (linecount==2 and distmax==-1) distmax=atof(vec[1].c_str());
            if (linecount==3 and distmin==-1) distmin=atof(vec[1].c_str());
            if (linecount==4 and cacbbb.empty()) cacbbb=vec[1];
            linecount++;
        }
    }
    pfh.close();

    // Check options
    if (!intrachain and !interchain){ cerr << "\nError: -x and -y are incompatible options\n" << endl; return 1; }
    if (!inputlist.empty() and plotinterpol){ cerr << "\nError: -l and -p are incompatible options\n" << endl; return 1; }
    if (plotinterpol and !cubicspline){ cerr << "\nError: -p option requires -c option\n" << endl; return 1; }
    if (distmin > distmax){ cerr << "\nError: -m arg must be greater than -n arg\n" << endl; return 1; }
    if (diffmin > diffmax){ cerr << "\nError: -j arg must be greater than -k arg\n" << endl; return 1; }
    if (diffmin < 0){ cerr << "\nError: -k argument must be positive\n" << endl; return 1; }
    if (diffmax < 0){ cerr << "\nError: -j argument must be positive\n" << endl; return 1; }
    if (distmax < 0){ cerr << "\nError: -m argument must be positive\n" << endl; return 1; }
    if (distmin < 0){ cerr << "\nError: -n argument must be positive\n" << endl; return 1; }
    if (decoys < 0){ cerr << "\nError: -s argument must be positive\n" << endl; return 1; }
    if (cacbbb.empty()){ cerr << "\nError: -r argument is missing\n" << endl; return 1; }
    if (!inputfile.empty() and !inputlist.empty()){
        cerr << "\nError: -i and -l arguments are incompatible\n" << endl; return 1;
    }
    if (!chainname.empty() and !inputlist.empty()){
        cerr << "\nError: -a and -l arguments are incompatible\n" << endl; return 1;
    }
    if (inputfile.empty() and inputlist.empty()){
        cerr << "\nError: -i or -l argument is missing\n" << endl; return 1;
    }
    if (simthres <= 0 or simthres >= 1){ cerr << "\nError: -t argument must belong to ]0,1[\n" << endl; return 1; }

    vector<string> possibler { "CA", "CB", "BB", "SC1", "allatom", "allatomCG", "sidechains", "sidechainsCG", "backbone" };
    if(find(possibler.begin(), possibler.end(), cacbbb) == possibler.end()) {
        cerr << "\nError: Invalid -r argument\n" << endl; return 1;
    }

    if(!fexists(rscripts+"/plotspline.R") and plotinterpol){ // Besides plotting, R script is not needed
        cerr << "\nError: R script \"plotspline.R\" not found in \"" << rscripts << "\"" << endl; return 1;
    }

    if (!jobdir.empty() and !plotinterpol and !write)
        cerr << "\nWarning: -o option useless without -p or -w options" << endl;

    if (zopt and cacbbb != "CA" and cacbbb != "BB"){
        cerr << "\nError: Z-score calculations allowed only for CA and BB representations" << endl; return 1;
    }

    if (windowsize != 12345){ // i.e. if user-defined windowsize (assuming that no user would type '-b 12345'...)
        if (windowsize > 0 and !write){ cerr << "\nError: -b argument useless without -w option\n" << endl; return 1; }
        if (windowsize < 1 and write){ cerr << "\nError: -b argument must be greater than 0\n" << endl; return 1; }
    }

    if (jobdir.empty())
        jobdir = (!inputfile.empty()) ? dirnaming(inputfile, indir, chainname) : dirnaming(inputlist, indir, chainname);
    /*************************
    * END OF OPTION HANDLING *
    *************************/


    map<string,string> three2one = create_map();

    // Extract the distance bins and push them into a vector
    vector<double> xvector;
    string xvecfile = indir+"/"+"xvector.dat";
    ifstream xfh(xvecfile);
    if(!xfh){ cerr << "\nError while opening file \"" << xvecfile << "\"" << endl; exit(1); }
    string distbin;
    while(getline(xfh, distbin)){
        xvector.push_back(atof(distbin.c_str()));
    }

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

    // Create a vector with every atom type
    vector<string> vecatypes;
    for (auto atype : atypes){
        vecatypes.push_back(atype);
    }
    sort(vecatypes.begin(), vecatypes.end());
    int atypessize = vecatypes.size();

    // For every atom pair (two loops)
    map<string, vector<double>> potentials; // This map is used for the cubic spline
    map<string, map<double, double>> nestedmap; // This (nested) map is used for the linear interpolation
    for (int i=0; i<atypessize; ++i){
        for (int j=i; j<atypessize; ++j){
            string pair = vecatypes[i]+vecatypes[j];
            string infile = indir+"/"+pair+".nrg";
            // Open the energy file
            ifstream fh(infile);
            if(!fh){ cerr << "\nError while opening potential file \"" << infile << "\"" << endl; exit(1); }
            string energy;
            int k = 0;
            while(getline(fh, energy)){
                // Push all the energies into a vector corresponding to a key (e.g. LV) in the map
                potentials[pair].push_back(atof(energy.c_str()));
                double mybin = xvector[k];
                nestedmap[pair][mybin] = atof(energy.c_str());
                k++;
            }
        }
    }

    // LINEAR INTERPOLATION
    // References for the linear_interpol() function
    vector<double>& refxvector(xvector);
    map<string, map<double, double>>& refnestedmap(nestedmap);
    int lastbinindex = xvector.size()-1;
    double binwidth = xvector[1]-xvector[0];
    double lastbin = xvector[lastbinindex];

    // CUBIC SPLINE
    // For each key of the map (= each atom pair)
    map<string, tk::spline> splinemap;// Will have the same keys as the "potentials" map
    if (cubicspline){
        if (plotinterpol){
            // Create ouput directory that will contain interpolated potentials
            createdir(jobdir);
            cout << "\nPlotting the interpolated potentials..." << endl;
        }

        int mycount = 0;
        for(auto& it : potentials){
            string aapair = it.first;
            // Compute a cubic spline based on the distance bins and corresponding energies
            tk::spline myspline;
            myspline.set_points(xvector, it.second);
            // Store the spline into a map
            splinemap[aapair] = myspline;

            if (plotinterpol){
                // Plot interpolated potential for the current AA pair
                // Cubic spline
                string outfilename = jobdir+"/"+aapair+".spline";
                ofstream ofh;
                ofh.open (outfilename);
                for (double bin=binwidth; bin<lastbin; bin+=binwidth)
                    ofh << bin << " " << myspline(bin) << endl;
                ofh.close();
                string plotcmd = "Rscript "+rscripts+"/plotspline.R "+outfilename+" "+aapair+"_spline "+jobdir+" > /dev/null";
                exec(plotcmd);

                // Linear interpolation
                outfilename = jobdir+"/"+aapair+".interpol";
                ofh.open (outfilename);
                for (double bin=binwidth; bin<lastbin; bin+=binwidth)
                    ofh << bin << " " << linear_interpol(bin, aapair, refxvector, refnestedmap, lastbinindex) << endl;
                ofh.close();
                plotcmd = "Rscript "+rscripts+"/plotspline.R "+outfilename+" "+aapair+"_interpol "+jobdir+" > /dev/null";
                exec(plotcmd);

                mycount++;
                cout.flush();
                cout << "\r" << mycount << "/210";
            }
        }
        if (plotinterpol)
            cout << endl << "Done" << endl;
    }

    vector<string> inputvec;
    if (!inputlist.empty()){
        ifstream lfh(inputlist);
        if(!lfh){ cerr << "\nError: File \"" << inputlist << "\" not found" << endl; return 1;}
        string item;
        while(getline(lfh, item))
            inputvec.push_back(item);
        lfh.close();
    }
    else
        inputvec.push_back(inputfile);

    int countall = 1;
    for (size_t l=0; l<inputvec.size(); ++l){
        string filename = inputvec[l];
        bool chainwarn = false;
        // If -l option used, check if there is a chain name
        if (!inputlist.empty()){
            chainname.erase(chainname.begin(), chainname.end());
            if (find(inputvec[l].begin(), inputvec[l].end(), ' ') != inputvec[l].end()){
                vector<string> vec = split(inputvec[l], ' ');
                if (vec[1].size() == 1 and isalnum(vec[1][0]))
                    chainname = vec[1];
                else if ( (vec[1].size() == 1 and !isalnum(vec[1][0])) or vec[1].size() > 1)
                    chainwarn = true;
                filename = vec[0];
            }
        }

        cout << "\nNow processing " << filename << endl;
        if (chainwarn) cerr << "Warning: Format of provided chain name is invalid\nWill process every chain found..." << endl;

        // If no user-defined name for output dir
        string myjobdir = (inputlist.empty()) ? jobdir : jobdir+"/"+dirnaming(filename, indir, chainname);

        // Push the file (only relevant ATOM lines) into a vector (faster when multiple reads)
        vector<string> filevec = file2vec(filename, three2one, atypes);

        // -q option: will only consider the residues selected by the user
        if (!residueslist.empty())
            filevec = residuefilter(filevec, residueslist);

        vector<string>& reffilevec(filevec);

        map<string, vector<double>> allsqdist;
        // Declare 4 vectors for the distances() function
        vector<double> sqdist; // distances
        vector<string> atompairs; // interacting pairs of atoms
        vector<string> resIDl; vector<string> resIDr; // resID: num+chain; l & r: interacting atom pairs (~left & right columns)

        // Declare 2 vectors for the sequence() function (only usefull for the Z-score calculation)
        vector<string> resID; // The protein residue ID (num+chain) sequence
        vector<string> aavec; // The protein AA (single letter code + atome name) sequence (same size as resID)

        // Declare their references
        vector<double>& refsqdist(sqdist); vector<string>& refatompairs(atompairs);
        vector<string>& refresIDl(resIDl); vector<string>& refresIDr(resIDr);
        vector<string>& refresID(resID); vector<string>& refaavec(aavec);

        // Fill the aavec and resID vectors
        sequence(reffilevec, refresID, refaavec, three2one); 

        // For the current PDB file: find the chains
        vector<string> chains = findchains(reffilevec);
 
        // If one (or several) particular chain has been chosen (-a option)
        if (!chainname.empty()){
            vector<string> validchains;
            // Check if (not) exists
            for (size_t i=0; i<chainname.size(); i++){
                string mychain(1,chainname[i]);
                if (find(chains.begin(), chains.end(), mychain) == chains.end())
                    cerr << "Warning: There is no chain " << mychain << endl;
                else
                    validchains.push_back(mychain);
            }

            if (validchains.size() == 0){
                cerr << "No valid chain in input; Will process every chain found..." << endl;
            }
            else{
                chains.clear();
                chains = validchains;
            }
        }

        cout << "\nChain(s) that will be processed: ";
        for (auto chain : chains)
            cout << chain << " ";
        cout << endl;
 
        if (chains.size() < 1)
            cerr << "Warning: No chain found in " << filename << endl;

        if (chains.size() > 1 and windowsize != 12345)
            cerr << "Warning: -b option not compatible with multichain calculs" << endl;

        if (windowsize == 12345)
            windowsize = 1;

        if (windowsize >= (int)aavec.size() and chains.size() > 1){
            cerr << "Error: -b argument must be smaller than the length of the protein sequence" << endl;
            continue;
        }

        /* For the current PDB file:
           For each combination of chain: */
        for (size_t i=0; i<chains.size(); ++i){
            for (size_t j=i; j<chains.size(); ++j){
                // Fill the six vectors
                if (interchain)
                    if (chains[i] != chains[j])
                        distances(reffilevec, refsqdist, refatompairs, refresIDl, refresIDr,
                                  distmax, distmin, chains[i], chains[j], diffmin, diffmax, three2one);
                if (intrachain) // NOT ELSE
                    if (chains[i] == chains[j])
                        distances(reffilevec, refsqdist, refatompairs, refresIDl, refresIDr,
                                  distmax, distmin, chains[i], chains[j], diffmin, diffmax, three2one);
            }
        }
    
        // Compute total energy
        double total_energy = 0;
        map<string, double> aaenergy; // Key: resID (num+chain); value: cumulated energy
        vector<string> data; // Tab-separated strings made of: distance, energy, AA pair, pair resID (num+chain)
        const int sqdistsize = sqdist.size();
        if (cubicspline)
            for(int i=0; i<sqdistsize; ++i){
                double myenergy = splinemap[atompairs[i]](sqdist[i]);
                total_energy+=myenergy;
                aaenergy[resIDl[i]]+=myenergy;
                aaenergy[resIDr[i]]+=myenergy;
                data.push_back(to_string(myenergy)+"\t"+to_string(sqdist[i])+"\t"+atompairs[i]+"\t"+resIDl[i]+"-"+resIDr[i]);
            }
        else
            for(int i=0; i<sqdistsize; ++i){
                double myenergy = linear_interpol(sqdist[i], atompairs[i], refxvector, refnestedmap, lastbinindex);
                total_energy+=myenergy;
                aaenergy[resIDl[i]]+=myenergy;
                aaenergy[resIDr[i]]+=myenergy;
                data.push_back(to_string(myenergy)+"\t"+to_string(sqdist[i])+"\t"+atompairs[i]+"\t"+resIDl[i]+"-"+resIDr[i]);
            }

        if (chains.size() > 0) 
            cout << "\nPseudo-energy = " << total_energy << endl;
        cout << "Done\n" << endl;
 
        // If write option activated
        if (write){
            if (!plotinterpol)
                createdir(myjobdir);
    
            ofstream tfh;
            cout << "Writing files \"data.tsv\", and \"energy_w" << to_string(windowsize) << ".tsv\"..." << endl;
            tfh.open (myjobdir+"/data.tsv");
            for (auto line : data)
                tfh << line << endl;
            tfh.close();


            tfh.open (myjobdir+"/energy_w"+to_string(windowsize)+".tsv");
            size_t limit = (size_t)(floor((double)windowsize/2));
            size_t start = (size_t)(ceil((double)windowsize/2)-1);

            for(size_t i=start; i<aavec.size()-limit; ++i){
                double localsum = 0;
                for(size_t j=i-start; j<=i+limit; j++){
                    localsum+=aaenergy[resID[j]];
                }
                double localaverage = localsum/(double)windowsize;
                tfh << aavec[i].substr(0,1) << "\t" << resID[i] << "\t" << localaverage << endl;
            }
            tfh.close();


            cout << "Done\n" << endl;
        }
 
 
        // If Z-score computation activated
        if (zopt and total_energy != 0){
            cout << "Computing the Z-score..." << endl;
 
            // Cannot compute a Z-score if the pseudo-energy is always 0
            if (chains.size() == 1 and !intrachain){
                cerr << "\nProgram (or iteration) stopped: No inter-chain interactions in a monomeric structure" << endl;
                continue;
            }
 
            const int resIDsize = resID.size(); // size of the protein sequence
            const double resIDsizedouble = (double)resIDsize;
            const int resIDlsize = resIDl.size(); // size of all pairwise interactions in the protein
 
            int counter = 0;
            vector<double> rand_energies;
            // For each random sequence decoys
            while(counter<decoys){
                // Shuffle the sequence of residues (which are identified by their number + chain)
                vector<string> resID_new(resID);
                vector<string> aavec_new(aavec); // This second vector will be used to calculate the similarity
                double similarity = 1.0;
                while (similarity > simthres){
                    int i = rand() % resIDsize + 0;
                    int j = rand() % resIDsize + 0;
 
                    if (i != j){
                        swap(resID_new[i], resID_new[j]);
                        swap(aavec_new[i], aavec_new[j]);
     
                        // After each residue swap, calculate the similarity with  the original sequence
                        int count = 0;
                        for (int k=0; k<resIDsize; ++k)
                            if (aavec_new[k] == aavec[k])
                                count++;
    
                        similarity = (double)count/resIDsizedouble;
                    }
                }
 
                // Shuffled map (will be used to create a new vector of residue pairs):
                // The amino acids (values) content is unchanged, but the correponding residue number+chain (keys) is different
                map<string,string> mapresID_new;
                for(int i=0; i<resIDsize; ++i)
                    mapresID_new[resID_new[i]] = aavec[i];
 
                // The new vector of residue pairs after shuffling 
                vector<string> aapairs_new;
                for (int i=0; i<resIDlsize; ++i){
                    string aa1 = mapresID_new[resIDl[i]];
                    string aa2 = mapresID_new[resIDr[i]];
 
                    if (aa1.compare(aa2) < 0) // if alphabet order
                        aapairs_new.push_back(aa1 + aa2);
                    else
                        aapairs_new.push_back(aa2 + aa1);
                }
 
                // Compute the energy of the decoys
                double rand_energy = 0;
                if (cubicspline)
                    for(int i=0; i<sqdistsize; ++i)
                        rand_energy+=splinemap[aapairs_new[i]](sqdist[i]);
                else
                    for(int i=0; i<sqdistsize; ++i)
                        rand_energy+=linear_interpol(sqdist[i], aapairs_new[i], refxvector, refnestedmap, lastbinindex);

                // Push the energy into a vector
                rand_energies.push_back(rand_energy);
                counter++;
                if (!mute){
                    cout.flush();
                    cout << "\r" << counter << "/" << decoys;
                }
            }
            cout << endl;
 
            // Compute the Z-score (requires sum, mean, and standard deviation)
            double sum = accumulate(rand_energies.begin(), rand_energies.end(), 0.0);
            double mean = sum / rand_energies.size();
            vector<double> diff(rand_energies.size());
            transform(rand_energies.begin(), rand_energies.end(), diff.begin(), [mean](double x) { return x - mean; });
            double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            double stdev = sqrt(sq_sum / rand_energies.size());
 
            double zscore = (total_energy-mean)/stdev;
 
            cout << "Z-score = " << zscore << "\nDone\n" << endl;
        }

        else if (zopt and total_energy == 0){
            cout << "Z-score = 0" << endl;
        }

        cerr.flush();
        if (inputfile.empty())
            cerr << "\r" << countall << "/" << inputvec.size();
        countall++;
    }
    cerr << endl;

    return 0;
}


/************
* FUNCTIONS *
************/

void createdir(string jobdir){
    cout << "\nCreating output directory \"" << jobdir << "\" ..." << endl;
    mkdir(jobdir.c_str(), ACCESSPERMS);
    cout << "Done" << endl;
    return;
}


double linear_interpol(double x, string atompair, vector<double>& xvector, map<string, map<double, double>>& nestedmap, int lastbinindex){

    while(xvector[lastbinindex] > x){
        lastbinindex--;
    }

    double x0 = xvector[lastbinindex];
    double x1 = xvector[lastbinindex+1];
    double y0 = nestedmap[atompair][x0];
    double y1 = nestedmap[atompair][x1];

    double interpol = y0 + (x-x0) * (y1-y0) / (x1-x0);

    return interpol;
}


void sequence(vector<string>& filevec, vector<string>& resID, vector<string>& aavec, map<string,string> three2one){
    int prevnum = -999;
    for (auto line : filevec){
        string oneletter = three2one[line.substr(17,3)];
        string chainfound = line.substr(21,1);

        string atom = line.substr(12,4);
        atom.erase(remove(atom.begin(),atom.end(),' '),atom.end());
        string atomname = oneletter+atom;

        string alternate = line.substr(16,1);
        int number = atoi(line.substr(22,4).c_str());

        if (alternate == " " or alternate == "A"){
            if (number != prevnum){
                resID.push_back(to_string(number)+chainfound);
                aavec.push_back(atomname);
                prevnum = number;
            }
        }
    }
}


void distances(vector<string>& filevec, vector<double>& sqdist, vector<string>& atompairs, vector<string>& resIDl, vector<string>& resIDr,
               double distmax, double distmin, string chain1, string chain2, int diffmin, int diffmax, map<string,string> three2one){
/* This functions
   - reads the vector containing (part of) the PDB file
   - calculates all the S̶Q̶U̶A̶R̶E̶D̶ interatomic distances for the one or two chains
   - push these distances into a referenced vector
   - there are 2 other referenced vectors: 1 for the atoms names, 1 for the atoms names+numbers
*/

    vector<int> num1; vector<int> num2;
    vector<double> x1; vector<double> x2;
    vector<double> y1; vector<double> y2;
    vector<double> z1; vector<double> z2;
    vector<string>atom1; vector<string>atom2;

    for (auto line : filevec){
        int number = atoi(line.substr(22,4).c_str());
        string chainfound = line.substr(21,1);

        string oneletter = three2one[line.substr(17,3)];
        string atomname = line.substr(12,4);
        atomname.erase(remove(atomname.begin(),atomname.end(),' '),atomname.end());
        atomname = oneletter+atomname;

        string alternate = line.substr(16,1);
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
    double sqd;
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
                        resIDl.push_back(to_string(num1[i])+chain1);
                        resIDr.push_back(to_string(num2[j])+chain2);
                    }
                    else{
                        atompairs.push_back(atom2[j]+atom1[i]);
                        resIDl.push_back(to_string(num2[j])+chain2);
                        resIDr.push_back(to_string(num1[i])+chain1);
                    }
                }
            }
        }
    }

    return;
}


string dirnaming(string filename, string indir, string chainname){
    // Output directory named after -i, -a, and -d options
    int pos1 = filename.rfind("/");
    int pos2 = filename.rfind(".");
    string name1 = filename.substr(pos1+1, pos2-pos1-1);
    string name2 = indir;
    if (indir.back() == '/'){
        //name2.pop_back(); // requires C++11
        name2.erase(name2.size()-1); // no need for C++11
    }
    int pos3 = name2.rfind("/");
    int pos4 = name2.size()-1;
    name2 = name2.substr(pos3+1, pos4-pos3);
    string dirname = name1+chainname+"_"+name2;

    return dirname;
}


vector<string> residuefilter(vector<string> filevec, string residueslist){
    vector<string> residuesvec = split(residueslist, ',');
    vector<string> outvec;

    for (auto line : filevec){
        string chain = line.substr(21,1);
        string number = line.substr(22,4);
        number.erase(remove(number.begin(),number.end(),' '),number.end());
        string id = number + chain;
        if (find(residuesvec.begin(), residuesvec.end(), id) != residuesvec.end()){
            outvec.push_back(line);
        }
    }

    return outvec;
}
