#ifndef _rsRNASP_score
#define _rsRNASP_score

#include <unordered_map>
// #include "structure.h"
// #include "PDB.h"
#include "rsRNASP_PDB.h"
#include "pdb_utils.h" 
#include "misc.h"
#include "xvec.h"

#include <fstream>


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
