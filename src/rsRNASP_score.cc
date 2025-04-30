#include "rsRNASP_score.h"
 



void rsRNASP_score::init_rsRNASP_etable(const string &energy_dir){

    std::string rsRNASP_file = energy_dir + "/dihedral/long.dat";//=============
 
    // if (! boost::filesystem::exists(rsRNASP_file)){
    if (!misc::file_exists(rsRNASP_file)){
        cout << "# File not found: " << rsRNASP_file.c_str() << endl;
        exit(1);
    }

    std::ifstream fs(rsRNASP_file);
    // std::ifstream fs(fn);
    std::string line;
    string ss[2];
    int count = 0;

    double param_bin_size = 0.3;  
    double param_max_dist = 24;

    while (std::getline(fs, line)){

        // read params from energy file

        if(line.substr(0,1) == "#"){
		
             continue;
        }
        else{
            std::istringstream iss(line);
            iss >> ss[0] >> ss[1] ;
            // atomtype is initialized by init_atomtype();
            // int id1 = param_p->atomtype[ss[0]];
            // int id2 = param_p->atomtype[ss[1]];
            int id1 = pdb_utils::get_atomtype_int(ss[0]);
            int id2 = pdb_utils::get_atomtype_int(ss[1]);

            // if (count == 5) printf("%d %d\n", id1,id2);
            for (int m = 0; m < mbin; m++)
            {
                iss >> rsRNASP[id1][id2][m];
                // printf(" %d %f\n", m, rsRNASP[id1][id2][m]);
            }
        }
    }
    fs.close();
    
    std::string short_file = energy_dir + "/dihedral/short.dat";//=============
 

    std::ifstream fs2(short_file);

//    int count = 0;


    while (std::getline(fs2, line)){


        if(line.substr(0,1) == "#"){
		
             continue;
        }
        else{
            std::istringstream iss(line);
            iss >> ss[0] >> ss[1] ;

            int id1 = pdb_utils::get_atomtype_int(ss[0]);
            int id2 = pdb_utils::get_atomtype_int(ss[1]);

            // if (count == 5) printf("%d %d\n", id1,id2);
            for (int m = 0; m < 43; m++)
            {
                iss >> potential1[id1][id2][m];

            }
        }
    }
    fs2.close();    
    
};







double rsRNASP_score::calc_rsRNASP(rna &astr){

    double etotal = 0.0;
    double energy1 = 0.0, energy2 = 0.0;
    int npair = 0;

    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & bchain : astr.chains){
            if (bchain.get_chaintype() != "NT") continue;
            for (auto & ares : achain.residues){
                for (auto & bres : bchain.residues){

                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()-4) continue;
                    }

                    if (achain.get_chainid() > bchain.get_chainid()) continue;

                    atom * c4a = ares.get_atom("C4'");
                    if (!c4a) continue;
                    atom * c4b = bres.get_atom("C4'");
                    if (!c4b) continue;
                    if (c4a->distance2(*c4b) > 9 * this->max_dist * this->max_dist) continue;
         

                    {
                        for(auto & aatom : ares.atoms){
                            // atom_num++;
                            for (auto & batom: bres.atoms){
                                double d2 = aatom.distance2(batom);

                                if (d2 >= this->max_dist * this->max_dist) continue;
                                // int b = int(d * 2);
                                
                                

                                double *x1 = aatom.get_x();
                                double *x2 = batom.get_x();

                                if (abs(x1[0] - x2[0]) >  this->max_dist) continue;
                                else if (abs(x1[1] - x2[1]) >  this->max_dist) continue;
                                else if (abs(x1[2] - x2[2]) >  this->max_dist) continue;

                                int b = this->distance2bin(sqrt(d2), this->bin_size);
                                // int b = int(d * 2);
                                // if (b >= 30) continue;

                                int atype_ia=-1, atype_ib=-1;

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());

                                // std::cout << aatom.get_type() << " " <<aatom.get_type() <<endl;

                                if (atype_ia > -1 && atype_ib > -1){
                                    energy2 += rsRNASP[atype_ia][atype_ib][b];
                                 //   cout << atype_ia << "  " << atype_ib << "    " << b << endl;
                                    npair ++;
                                // residue level energy
                                    // ares.rsRNASP_vanilla += 0.5 *this->rsRNASP[atype_ia][atype_ib][b];
                                    // bres.rsRNASP_vanilla += 0.5 *this->rsRNASP[atype_ia][atype_ib][b];

                                }
                            }
                        }

                    }
                }
            }
        }

    }


    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & bchain : astr.chains){
            if (bchain.get_chaintype() != "NT") continue;
            for (auto & ares : achain.residues){
                for (auto & bres : bchain.residues){

                    if (achain.get_chainid() == bchain.get_chainid() && (bres.get_residsd()-ares.get_residsd()==1)  ){
                       


                

                    atom * c4a = ares.get_atom("C4'");
                    if (!c4a) continue;
                    atom * c4b = bres.get_atom("C4'");
                    if (!c4b) continue;
                    if (c4a->distance2(*c4b) > 9 * this->max_dist * this->max_dist) continue;
         

                    {
                        for(auto & aatom : ares.atoms){
                            // atom_num++;
                            for (auto & batom: bres.atoms){
                                double d2 = aatom.distance2(batom);

                                if (d2 >= 13*13) continue;
                                // int b = int(d * 2);
                                
                                

                                double *x1 = aatom.get_x();
                                double *x2 = batom.get_x();

                                if (abs(x1[0] - x2[0]) >  this->max_dist) continue;
                                else if (abs(x1[1] - x2[1]) >  this->max_dist) continue;
                                else if (abs(x1[2] - x2[2]) >  this->max_dist) continue;

                                int b = this->distance2bin(sqrt(d2), this->bin_size);
                                // int b = int(d * 2);
                                // if (b >= 30) continue;

                                int atype_ia=-1, atype_ib=-1;

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());

                                // std::cout << aatom.get_type() << " " <<aatom.get_type() <<endl;

                                if (atype_ia > -1 && atype_ib > -1 && b <43){
                                    energy1 += potential1[atype_ia][atype_ib][b];
  //                                  cout << atype_ia << "  " << atype_ib << "    " << b << "     "<< potential1[atype_ia][atype_ib][b] << endl;
                                    npair ++;
                                // residue level energy
                                    // ares.rsRNASP_vanilla += 0.5 *this->rsRNASP[atype_ia][atype_ib][b];
                                    // bres.rsRNASP_vanilla += 0.5 *this->rsRNASP[atype_ia][atype_ib][b];

                                }
                            }
                        }
                    }
                    }
                }
            }
        }

    }
    int residue_num=0;
    for (auto & achain : astr.chains){
      for (auto & ares : achain.residues){
       residue_num++;
       }
       }
    

    etotal = energy1 + 120*energy2/rsRNASP_score::fun(residue_num);
    astr.n_atompairs = npair;
//std::cout << energy1 << "   " <<16*energy2/fun(24)<<"    " <<residue_num<< endl;
    return etotal;

}




