#ifndef _rsRNASP_CALC
#define _rsRNASP_CALC

#include <unordered_map>
// #include "structure.h"
// #include "PDB.h"
#include "rsRNASP_PDB.h"
#include "pdb_utils.h" 
#include "misc.h"
#include "xvec.h"

#include <fstream>
// #include <boost/filesystem.hpp>
// #include "SMCRA_param.h" // already included in structure.h
// #include "boost/multi_array.hpp" // moved to rsRNASP_dipolar.h

// ========================
// vanilla rsRNASP atom-atom distance term for rna
//
// ========================

// #define DEBUG

class rsRNASP_calculator{
protected:
// ==============================================
// common

    // std::string rsRNASP_fpath;

    // deprecated, moved to SMCRA_param.h, use param_p
    // std::unordered_map<std::string,int> atomtype;
    // std::unordered_map<int,std::string> atomtype2;
    // void init_rna_atomtype();
    // int get_atomtype_int(std::string);
    // smcra_param * param_p;

// ==============================================
// rsRNASP vanilla
// params and etable
// modified: 181115 
// make alpha, bin_size, max_dist adjustable for better tweaking.
// ==============================================
     
    static constexpr int mbin=100, matype = 85;


    double ALPHA = 1.61;
    double bin_size = 0.3;
    double max_dist = 24.0;
    int max_bin = 80;




    double rsRNASP[matype][matype][mbin]={};
    double potential1[matype][matype][43]={};
    int count[matype][matype][mbin]={};
    double rsRNASP_c[matype][matype][mbin]={};
// ==============================================
// rsRNASP vanilla
// init
// ==============================================

    // void init_rsRNASP();
    // void init_rsRNASP(std::string);
    // void init_rsRNASP_etable();
    


// ==============================================
// rsRNASP vanilla
// training
// ==============================================


    
    void count_to_rsRNASP();
    // void count_pairs(structure &);
    void count_pairs(rna &);

    bool is_far_away(atom & atom_a, atom & atom_b, double maxdist);

    // 
    // void count_to_3drna();
    // double calc_3drna();

    // double calc_rsRNASP(chain & achain); // deprecated
    // double calc_rsRNASP(structure & astr); // deprecated, use rsRNASP_PDB::rna instead;

    double calc_rsRNASP(rna &astr);





public:

// ==============================================
// rsRNASP vanilla
// use get_rsRNASP(structure & astr) to get rsRNASP energy
// ==============================================

    // used for training
    // rsRNASP_calculator(){param_p = new smcra_param();};
    rsRNASP_calculator(){};
    int k_d=105;

    rsRNASP_calculator(const std::string &energy_dir){
        // param_p = new smcra_param();
        init_rsRNASP_etable(energy_dir);
    };
    void init_rsRNASP_etable(const std::string &energy_dir);




    ~rsRNASP_calculator(){};
    void train_rsRNASP(string pdbl_fpath, double alpha, double bin_size, double max_dist);

// alpha, bin_size, max_dist
// modified: 181115 
    void set_abm(double param_alpha, double param_bin_size, double param_max_dist){
        this->ALPHA = param_alpha;
        this->bin_size = param_bin_size;

        double residual = fmod(param_max_dist, param_bin_size);
        int mm = (int) (param_max_dist / param_bin_size);

        // adjust the last bin to make it a complete bin
        if (residual > 1e-6) {
            cerr << "# max_dist adjusted to upperbound" << endl;
            mm ++;
        }

        this->max_bin = mm;
        this->max_dist = mm * bin_size;
    }

// translate distance to bin, based on bin_size
// modified: 181115 
    inline int distance2bin(double dist, double bin_size){
        int bin = (int) (dist / bin_size);
        return bin;
    }

    void print_parameters(FILE * out);

    void print_ersRNASP(FILE *out);
    void print_count(FILE *out);


    // double get_rsRNASP(structure & astr){return calc_rsRNASP(astr);}; //deprecated, use rsRNASP_PDB
    double get_rsRNASP(rna & astr){return calc_rsRNASP(astr);};
    void get_rsRNASP_per_res(structure & astr, FILE *out);
    

    // similar to 3drna;
    void train_3drna(string pdbl_fpath);
    void count_to_3drna();
    double  calc_3drna(rna &astr);
    double get_3drna(rna & astr){return calc_3drna(astr);};



};

#endif
