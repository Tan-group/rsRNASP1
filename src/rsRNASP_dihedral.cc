#include "rsRNASP_dihedral.h"

void rsRNASP_dihedral::print_ersRNASP_dih(FILE * out) {
    // FILE * out;
    // out = fopen(out.c_str(),"w");

    for (int i=0;i<dtype;i++){
        fprintf(out,"type_%d",i);
        for(int j=0;j<dbin;j++){
                fprintf(out,"\t%8.6f",rsRNASP_dih[i][j]);
        }
        fprintf(out,"\n");
    }
 
}


void rsRNASP_dihedral::print_count_dih(FILE * out){
    for (int i=0;i<dtype;i++){
        fprintf(out,"type_%d",i);
        for(int j=0;j<dbin;j++){
                fprintf(out,"\t%d",count_dih[i][j]);
        }
        fprintf(out,"\n");
    }
}


void rsRNASP_dihedral::train_rsRNASP_dih(string pdbl_fpath, int length_limit){


    string train_list=pdbl_fpath;

    std::ifstream fs(train_list.c_str());
    std::string line;

    while(std::getline(fs, line)){

        cout << line <<endl;
        // structure *astr = new structure();
        rna *astr = new rna(); 
        astr->readpdb(line);
        cout << astr->get_nres() <<endl;

        // remove > 500 nt structures

        if (astr->get_nres() > length_limit){
            cout << "skipped" << endl;
            continue;
        }
        astr->rna_init_dihedrals(); 

        count_dihs(*astr);
    }

    count_to_rsRNASPdih();

}

void rsRNASP_dihedral::init_rsRNASP_etable(string energy_dir){

    std::string rsRNASP_file = energy_dir + "/energy.dih.dat";

    // if (! boost::filesystem::exists(rsRNASP_file)){
    if (! misc::file_exists(rsRNASP_file)){
        cerr << "# File not found:" << rsRNASP_file.c_str() << endl;
        exit(1);
    }

    std::ifstream fs(rsRNASP_file);
    std::string line;
    string ss[2];

    int row=0;

    while (std::getline(fs, line)){
        std::istringstream iss(line);
        iss >> ss[0];
        for (int m = 0; m < dbin; m++)
        {
            iss >> rsRNASP_dih[row][m];
        }
        row++;
    }
};

double rsRNASP_dihedral::calc_rsRNASP_dih(rna &astr){ 

    double total_energy=0.0;

    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & ares : achain.residues){ 
            if (ares.has_dih()){
                double *dihs = ares.get_dihedrals(); 
                for (int i = 0; i < 7; i++){
                    double tmpdih=dihs[i];
                    if (tmpdih < 0) tmpdih = 360 + tmpdih;
                    int dihbin = (int)tmpdih/this->dbin_width;
                    total_energy += rsRNASP_dih[i][dihbin];
                }
            }
        }
    }

    return total_energy;
}

// double rsRNASP_dihedral::dihedral(atom &aa, atom &ab, atom &ac, atom &ad){
//     Xvec x1(aa.get_x());
//     Xvec x2(ab.get_x());
//     Xvec x3(ac.get_x());
//     Xvec x4(ad.get_x());

//     double tor = torsion_xvec(x1,x2,x3,x4);

//     return tor;

// }

void rsRNASP_dihedral::count_dihs(rna & astr){ 

    for (auto & achain : astr.chains){
        for (auto & ares : achain.residues){ 
            if (ares.has_dih()){
                double *dihs = ares.get_dihedrals();
                for (int i = 0; i < 7; i++){
                    double tmpdih=dihs[i];
                    if (tmpdih < 0) tmpdih = 360 + tmpdih;

                    int dihbin = (int)tmpdih/this->dbin_width;

                    count_dih[i][dihbin]++;
                }
            }
        }
    }
}



void rsRNASP_dihedral::count_to_rsRNASPdih(){

    double dih_colsum[dbin]= {};
    double dih_rowsum[dtype]={};
    double dih_total=0;

    for (int i = 0; i < dtype; i++){
        for (int j = 0; j < dbin; j++){
            dih_colsum[j] += count_dih[i][j];
            dih_rowsum[i] += count_dih[i][j];
            dih_total += count_dih[i][j];
        }
    }

    if (this->dih_ref_type == 0) {
        for (int i = 0; i < dtype; i++){
            for (int j = 0; j < dbin; j++){
                double quot, out;

                quot = 1.0 * count_dih[i][j] * dih_total / dih_rowsum[i] / dih_colsum[j];

                if (quot < 0.000001) out = 10.0;
                else out = - 0.001987 * 300 * log(quot);

                rsRNASP_dih[i][j]=out;
            }
        }
    }
    else if (this->dih_ref_type == 1){
        for (int i = 0; i < dtype; i++){
            for (int j = 0; j < dbin; j++){
                double quot, out;

                quot = 1.0 * count_dih[i][j] / dih_rowsum[i] * 1.0 / dbin;

                if (quot < 0.000001) out = 10.0;
                else out = - 0.001987 * 300 * log(quot);

                rsRNASP_dih[i][j]=out;
            }
        }
    }

}
void rsRNASP_score::init_rsRNASP_etable(const string &energy_dir){

    std::string rsRNASP_file = energy_dir + "/dihedral/nonlocal.dat";
    //=============
 
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
    
    std::string short_file = energy_dir + "/dihedral/local.dat";//=============
 

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
    rsRNASP_dihedral tor;

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
    

    etotal = energy1 + 120*energy2/rsRNASP_score::fun(residue_num) + tor.energy;
    astr.n_atompairs = npair;
//std::cout << energy1 << "   " <<16*energy2/fun(24)<<"    " <<residue_num<< endl;
    return etotal;

}

