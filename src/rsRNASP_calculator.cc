#include "rsRNASP_calculator.h"
 


// ==============================================
// rsRNASP vanilla
// init etable
// 
// modified: 181115 add three parameters: alpha, bin_size, max_dist
// ==============================================

void rsRNASP_calculator::init_rsRNASP_etable(const string &energy_dir){

    std::string rsRNASP_file = energy_dir + "/dist/nonlocal.dat";//=============
 
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

    double param_alpha = 1;
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
    
    std::string short_file = energy_dir + "/dist/local.dat";//=============
 

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

// ==============================================
// 3drna
// training
// ==============================================


void rsRNASP_calculator::train_3drna(string pdbl_fpath){ 

    // this->rsRNASP={};
    // this->count={};

    string train_list=pdbl_fpath;

    std::ifstream fs(train_list.c_str());
    std::string line;

    while(std::getline(fs, line)){

        cout << line <<endl;
        // structure *astr = new structure();
        rna * astr = new rna(line);

        // astr->readpdb();

        cout << astr->get_nres() <<endl;

        if (astr->get_nres() > 500 || astr->get_nres() <20) continue;
        // chain &a_chain=astr->chains[0];
        count_pairs(*astr);
    }

    this->count_to_3drna();
    // print_ersRNASP(out_path);
}

void rsRNASP_calculator::count_to_3drna(){

    int row_sum[matype][matype] = {0};

    int col_sum[mbin] = {0};

    int total = 0;

    for (int i=0; i < matype;i++){
        for (int j=0;j <matype; j++){
            for (int k=0; k<30; k++){
                row_sum[i][j] += count[i][j][k];
                col_sum[k] += count[i][j][k];
                total += count[i][j][k];
            }
        }
    }

    for (int i=0; i < matype;i++){
        for (int j=0;j <matype; j++){
            for (int k=0; k<30; k++){
                double out, quot;

                // refc = 1.0 * row_sum[i][j] * col_sum[k] / total;
                if (col_sum[k] < 0.0000001) quot = 0;
                else quot = 1.0 *count[i][j][k] * total / row_sum[i][j] / col_sum[k];

                if (quot < 0.000001) out = 10.0; 
                else out = -.001987*300.0* log(quot);

                rsRNASP_c[i][j][k]= out;
            }
        }
    }
}

double rsRNASP_calculator::calc_3drna(rna &astr){

    double etotal = 0.0;

    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & bchain : astr.chains){
            if (bchain.get_chaintype() != "NT") continue;
            for (auto & ares : achain.residues){
                for (auto & bres : bchain.residues){

                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()) continue;
                    }

                    if (achain.get_chainid() > bchain.get_chainid()) continue;

                    {
                        for(auto & aatom : ares.atoms){
                            // atom_num++;
                            for (auto & batom: bres.atoms){

                                double d = aatom.distance(batom);
                                int b = int(d * 2);
                                if (b >= 30) continue;

                                int atype_ia=-1, atype_ib=-1;

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());

                                // std::cout << aatom.get_type() << " " <<aatom.get_type() <<endl;

                                if (atype_ia > -1 && atype_ib > -1){

                                    if (pdb_utils::is_rna_sideatom(aatom.get_type()) || pdb_utils::is_rna_sideatom(batom.get_type()))
                                    {
                                        etotal += 2.5 * rsRNASP[atype_ia][atype_ib][b];
                                    }
                                    else {
                                        etotal += rsRNASP[atype_ia][atype_ib][b];
                                    }
                                    
                                    // etotal += rsRNASP[atype_ia][atype_ib][b];
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

    return etotal;

}



// ==============================================
// rsRNASP vanilla
// training
// ==============================================



void rsRNASP_calculator::train_rsRNASP(string pdbl_fpath, double alpha, double bin_size, double max_dist){ 

    // this->rsRNASP={};
    // this->count={};

    set_abm(alpha, bin_size, max_dist);

    string train_list=pdbl_fpath;

    std::ifstream fs(train_list.c_str());
    std::string line;

    while(std::getline(fs, line)){

        cout << line <<endl;
        // structure *astr = new structure();
        rna * astr = new rna(line);

        // astr->readpdb();

        cout << astr->get_nres() <<endl;

        // length filter is done manually.

        // if (astr->get_nres() > 500 || astr->get_nres() <20) continue;
        // if (astr->get_nres() > 500) continue;
        // chain &a_chain=astr->chains[0];
        count_pairs(*astr);

        delete astr;
    }

    count_to_rsRNASP();
    // print_ersRNASP(out_path);
}

// ==============================================
// rsRNASP vanilla
// training
// ==============================================


// void rsRNASP_calculator::count_pairs(structure &astr){ 

double fun(int n)
{
 return -2685.0/sqrt(n+16) + 542.0;
}


bool is_far_away(atom & atom_a, atom & atom_b, double maxdist){
    double *x1 = atom_a.get_x();
    double *x2 = atom_b.get_x();

    bool flag = false;

    if (x1[0] - x2[0] > maxdist) flag = true;
    else if (x1[1] - x2[1] > maxdist) flag = true;
    else if (x1[2] - x2[2] > maxdist) flag = true;

    return flag;



}


void rsRNASP_calculator::count_pairs(rna &astr){ 
// debug
    int pair_a=0;
    int pair_ab=0;
    // int pair_b = 0;


    for (auto & achain : astr.chains){
        for (auto & bchain : astr.chains){
           for (auto & ares : achain.residues){
                for (auto & bres : bchain.residues){
                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()) continue;
                    }

                    if (achain.get_chainid() > bchain.get_chainid()) continue;
                    // skip residue pairs with c4' distance > 3 * max dist,
                    atom * c4a = ares.get_atom("C4'");
                    if (!c4a) continue;
                    atom * c4b = bres.get_atom("C4'");
                    if (!c4b) continue;
                    if (c4a->distance2(*c4b) > 9 * this->max_dist * this->max_dist) continue;
         
                    {
                        for(auto & aatom : ares.atoms){
                            // atom_num++;
                            for (auto & batom: bres.atoms){

                                // skip far away atoms

                                double *x1 = aatom.get_x();
                                double *x2 = batom.get_x();

                                if (abs(x1[0] - x2[0]) >  this->max_dist) continue;
                                else if (abs(x1[1] - x2[1]) >  this->max_dist) continue;
                                else if (abs(x1[2] - x2[2]) >  this->max_dist) continue;

                                // use distance2

                                double d2 = aatom.distance2(batom);

                                // double d = aatom.distance(batom);

                                if (d2 >= this->max_dist * this->max_dist) continue;
                                // int b = int(d * 2);
                                int b = this->distance2bin(sqrt(d2), this->bin_size);

                                int atype_ia=-1, atype_ib=-1;

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());

                                // std::cout << aatom.get_type() << " " <<aatom.get_type() <<endl;

                                if (atype_ia > -1 && atype_ib > -1){
                                    // debug

                                    if (achain.get_chainid() != bchain.get_chainid()){
                                        pair_ab ++;
                                    }
                                    else pair_a ++;

                                    this->count[atype_ia][atype_ib][b]++;
                                    this->count[atype_ib][atype_ia][b]++;
                                }
                            }
                        }

                    }
                }
            }
        }

    }

    cout << "Atom pairs on one chain: " << pair_a << "Atom pairs on >2 chains:" << pair_ab << endl;
}


// ==============================================
// rsRNASP vanilla
// training
// ==============================================



void rsRNASP_calculator::count_to_rsRNASP(){
    for (int i=0; i < matype;i++){
        for (int j=0;j <matype; j++){
            for (int k=0; k<max_bin; k++){
                double r_rat, out, quot;

                r_rat = pow((k * this->bin_size + 0.5 * this->bin_size)/this->max_dist,ALPHA);

                if (r_rat < 0.000001) quot = 0.0;
                else quot = count[i][j][k]/(r_rat*count[i][j][this->max_bin - 1]); 

                if (quot < 0.000001) out = 10.0;
                else out = -.001987*300.0*log(quot);

                rsRNASP_c[i][j][k]= out;
            }
        }
    }
}

// ==============================================
// rsRNASP vanilla
// print
// modified: 181115 
// Add three lines for  alpha, bin_size, max_dist
// e.g
// # alpha 1.61
// # bin_size 0.5
// # max_dist 15
// ==============================================

void rsRNASP_calculator::print_parameters(FILE * out){
    fprintf(out,"################ Parameters #################\n");
    fprintf(out,"# alpha %f\n", this->ALPHA);
    fprintf(out,"# bin_size %f\n", this->bin_size);
    fprintf(out,"# max_dist %f\n", this->max_dist);
    fprintf(out,"#############################################\n");

}

void rsRNASP_calculator::print_ersRNASP(FILE * out){


    print_parameters(out);

    for (int i=0;i<matype;i++){
        for(int j=0;j<matype;j++){
            fprintf(out,"%s\t%s",pdb_utils::typeatom[i].c_str(),pdb_utils::typeatom[j].c_str());
            for (int k=0;k < max_bin;k++){
                fprintf(out,"\t%f",rsRNASP_c[i][j][k]);
            }
            fprintf(out,"\n");
        }
    }
}

void rsRNASP_calculator::print_count(FILE *out){

    print_parameters(out);

    for (int i=0;i<matype;i++){
        for(int j=0;j<matype;j++){
            fprintf(out,"%s\t%s",pdb_utils::typeatom[i].c_str(),pdb_utils::typeatom[j].c_str());
            for (int k=0;k < max_bin;k++){
                fprintf(out,"\t%d",count[i][j][k]);
            }
            fprintf(out,"\n");
        }
    }

}



double rsRNASP_calculator::calc_rsRNASP(rna &astr){

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

                    if (achain.get_chainid() == bchain.get_chainid() && (bres.get_residsd()-ares.get_residsd()>2) && (bres.get_residsd()-ares.get_residsd()<5) ){
                       


                

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
//                                    cout << atype_ia << "  " << atype_ib << "    " << b << "     "<< potential1[atype_ia][atype_ib][b] << endl;
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
    

    etotal = energy1 + 16*energy2/fun(residue_num);
    astr.n_atompairs = npair;
//std::cout << energy1 << "   " <<16*energy2/fun(24)<<"    " <<residue_num<< endl;
    return etotal;

}

// void rsRNASP_calculator::get_rsRNASP_per_res(structure & astr, FILE *out){
//     calc_rsRNASP(astr);

//     for (auto & achain: astr.chains) {
//         fprintf(out, "Chain: %8s\n", achain.get_chainid().c_str());
//         for (auto & ares: achain.residues){
//             fprintf(out, "%8d\t%8s\t%8.6f\n", ares.get_residsd(), ares.get_resname().c_str(),ares.rsRNASP_vanilla);
//         }
//     }
// }



