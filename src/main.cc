#include "main.h"

// ========================
// main test
//
// ========================

#define MP

namespace PARAM{
    vector<string> pdblist;
    string rsRNASP_folder;
    bool norm = false;
}
 

bool calc_rsRNASP_list(vector<string> pdblist, string rsRNASP_folder){

    // rsRNASP_calculator dc(rsRNASP_folder);
    rsRNASP_calculator * dc = new rsRNASP_calculator(rsRNASP_folder);
    rsRNASP_score * dp = new rsRNASP_score(rsRNASP_folder);
    
    // rsRNASP_dihedral ddih(rsRNASP_folder); 
    // rsRNASP_dipolar ddipo(rsRNASP_folder);
    // string rsRNASP_dih_folder = "/home/tc/GIT/rsRNASP_rna2/tmpdata";
    // dc.init_rsRNASP_dih(rsRNASP_dih_folder);

    int ne = pdblist.size();
    vector<string> sdim (ne, "");
    FILE *fp = stdout;

#if defined MP
    #pragma omp parallel for
#endif
    for (int i = 0; i < ne; i++){
        double energy = 0;

        string line = pdblist[i];
        // boost::filesystem::path p(line);



        // if (!boost::filesystem::exists(p)){
        if (!misc::file_exists(line)){
            cerr << "pdb " << line << " not exist, skipped" << endl;
            continue;
        }

        // structure *astr = new structure();
        rna *astr = new rna(line);
        // astr->readpdb(line);
        // astr->chains[0].init_dihedrals();
        // astr->rna_init_dihedrals();
        // astr->rna_init_polar();
//        std::cout << dp->get_rsRNASP(*astr)<<endl;
        // double energy_td = dc.calc_3drna(*astr);
       
	if (ne < dc->k_d){
        energy = dc->get_rsRNASP(*astr);
        }
        else{
        energy = dp->get_rsRNASP(*astr);
        }       
//        else{
//        double energy = dp->get_rsRNASP(*astr);
//        }

        double energy_dipo_pn = 0.;
        double energy_dipo_pp = 0.;
        double energy_dih = 0.;

        {
        // double energy_dih = ddih.get_rsRNASP(*astr); 
        // std::pair<double, double> energy_dipo = ddipo.get_rsRNASP(*astr);
        // double energy_dipo_pn = 0.energy_dipo.first;
        // double energy_dipo_pp = 0.energy_dipo.second;
        }






        char sout[200];
 
        sprintf (sout, "%s %8.6f \n", line.c_str(), energy);


#if defined MP
        sdim[i] = sout;

#else
        fprintf(fp,"%s",sout);
#endif
        delete astr;

    }


#if defined MP
    for (int i = 0 ; i < ne; i++){
        fprintf(fp, "%s", sdim[i].c_str());

    }
#endif

    delete dc;

    return 0;
}

void print_help(char const *argv[]){
        std::cout << "#######################################################" << endl;
        std::cout << "# Calculate rsRNASP_rna score for a pdb or a list of pdbs" << endl;
        std::cout << "#######################################################" << endl;
        std::cout << "Usage: " << argv[0] <<" pdb " << endl;
        std::cout << "   or: " << argv[0] <<" [ options ] " << endl;


        std::cout << "Options:" << endl;
        std::cout << "   pdb [ pdb2 pdb3 ...], input RNA structures in pdb format" << endl;
        std::cout << "   -d directory,         OPTIONAL, override default directory of energyfiles" << endl;
        std::cout << "                         default:" << PARAM::rsRNASP_folder  << endl;
        std::cout << "   -norm                 normalize rsRNASP score by RNA length (experimental)" << endl;
        std::cout << "   -l pdblist,           A list of absolute paths to pdb files (plain text) UTF-8 encoding" << endl;

        exit(1);
}

void read_param(int argc, char const *argv[]){

    string projfolder = getenv("rsRNASP_RNA_HOME");

    PARAM::rsRNASP_folder = projfolder+ "/data/energyfiles";

    // vector<string> pdblist;

    for(int i=1; i<argc; i++){
        string opt = argv[i];
        if (opt == "-d") {
            i++;
            PARAM::rsRNASP_folder = argv[i];
        }
        else if (opt == "-l") {
            i++;

            std::ifstream fs(argv[i]);
            std::string line;

            while(std::getline(fs, line)){
                PARAM::pdblist.push_back(line);
            }
        }
        else if (opt == "-s"){
            continue;
        }
        else if (opt == "-norm"){
            PARAM::norm = true;
            continue;
        }
        else if (argv[i][0] == '-'){
            cerr << "unknown option: " << argv[i] << endl;
        }
        else {
            if (misc::file_exists(argv[i])){
                PARAM::pdblist.push_back(argv[i]);
                // cerr << "pdb added: " << argv[i]<< endl;
            }
            else {
                cerr << "File not exist: " << argv[i] << endl;
            }
            
        }
    }

    if (argc < 2){
        print_help(argv);

    }
    if (PARAM::rsRNASP_folder.length() < 3){
        cerr << "!! rsRNASP_rna root not found: " << PARAM::rsRNASP_folder  << endl;
    }
}



int run_rsRNASP(int argc, char const *argv[])
{
    //#######################
    // find executable
    // https://stackoverflow.com/a/49227138/7503413
    //#######################

    read_param(argc, argv);

    return calc_rsRNASP_list(PARAM::pdblist, PARAM::rsRNASP_folder);
}


int main(int argc, char const *argv[])
{
    run_rsRNASP(argc, argv);

    return 0;
}
