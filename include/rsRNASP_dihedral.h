#ifndef _rsRNASP_CALC_1
#define _rsRNASP_CALC_1
// rsRNASP dihedral term

// #include "PDB.h"
#include "rsRNASP_PDB.h"
#include "pdb_utils.h"
#include "misc.h"
#include "xvec.h"

#include <unordered_map>

//#include "rsRNASP_PDB.h"
//#include "pdb_utils.h" 
//#include "misc.h"
//#include "xvec.h"

#include <fstream>
// #include <unordered_map>
// #include "structure.h"
// #include "xvec.h"
// #include <boost/multi_array.hpp>
// #include <boost/filesystem.hpp>

// ================================
// rsRNASP dihedral term
//
// ================================

class rsRNASP_dihedral{
protected:
    std::string rsRNASP_fpath;
    // void init_rsRNASP_dih();
    // void init_rsRNASP_dih(std::string);

    void init_rsRNASP_etable(std::string);

    // void set_dih_reftype(int reftype){dih_ref_type = reftype;}
    // int get_dih_reftype(){return dih_ref_type;}

    static constexpr int dbin=80, dtype=7;
    static constexpr double dbin_width=4.5;
    int dih_ref_type = 0;

    int count_dih[dtype][dbin]={};
    double rsRNASP_dih[dtype][dbin]={};

    std::unordered_map<std::string,int> dihtype;
    std::unordered_map<int,std::string> dihtype2;

    double dihedral(atom &aa, atom &ab, atom &ac, atom &ad);

// dihedral term

    // reference type: 0 -- average, as in 3drnascore
    //                 1 -- uniform distribution
    // used in count_to_rsRNASPdih



    // void count_dihs(chain & achain); // deprecated
    // void count_dihs(structure & astr);
    void count_dihs(rna & astr);
    void count_to_rsRNASPdih();

    // uniform distribution
    // void count_to_rsRNASPdih_unibg();

    int get_dihtype_int(std::string);
    // void init_rsRNASP_dih(std::string);


    // double calc_rsRNASP_dih(chain & achain); // deprecated

    // double calc_rsRNASP_dih(structure & astr);
    double calc_rsRNASP_dih(rna & astr);



public:

    rsRNASP_dihedral(){};
    rsRNASP_dihedral(std::string energy_dir){
        init_rsRNASP_etable(energy_dir);
    };

    ~rsRNASP_dihedral(){};
    // ==============================================
    // training
    // ==============================================

    void train_rsRNASP_dih(string pdbl_fpath, int length_limit = 500);
    void print_ersRNASP_dih(FILE * out);
    void print_count_dih(FILE * out);

    double energy=0.0;


    // double get_rsRNASP(structure & astr){return calc_rsRNASP_dih(astr);};
    double get_rsRNASP(rna & astr){return calc_rsRNASP_dih(astr);};
    void get_rsRNASP_per_res(structure & astr, FILE *out);
};


class rsRNASP_score{
protected:

     
    static constexpr int mbin=100, matype = 85;


    double bin_size = 0.3;
    double max_dist = 24.0;
    int max_bin = 80;




    double rsRNASP[matype][matype][mbin]={};
    double potential1[matype][matype][43]={};
    int count[matype][matype][mbin]={};
    double rsRNASP_c[matype][matype][mbin]={};


    



    double calc_rsRNASP(rna &astr);





public:


    rsRNASP_score(){};
    int k_d=105;
    
    rsRNASP_score(const std::string &energy_dir){
        // param_p = new smcra_param();
        init_rsRNASP_etable(energy_dir);
    };
    void init_rsRNASP_etable(const std::string &energy_dir);




    ~rsRNASP_score(){};
    void train_rsRNASP(string pdbl_fpath, double alpha, double bin_size, double max_dist);

    double fun(int n){
 return -2685.0/sqrt(n+16) + 542.0;
}




    inline int distance2bin(double dist, double bin_size){
        int bin = (int) (dist / bin_size);
        return bin;
    }

 
    double get_rsRNASP(rna & astr){return calc_rsRNASP(astr);};
    void get_rsRNASP_per_res(structure & astr, FILE *out);
    





};


#endif
