//
//  vide_decoder.cpp
//  LDPCDecoders
//
//  Created by Jon on 5/10/22.
//

#include "vide_decoder.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

vide_decoder::vide_decoder(string filename, mt19937& engine, int h_bb, int h_cc) :
HGP_code(filename, engine),
h_bb(h_bb),
h_cc(h_cc)
{
    L.resize(n_q);
    neighbors_in_R.resize(n_q);
    to_visit.resize(n_q);
    set_L.resize(n_q);
    set_R.resize(n_chk_q);
}

void vide_decoder::reset() {
    fill(L.begin(), L.end(), 0);
    fill(to_visit.begin(), to_visit.end(), 0);
    fill(neighbors_in_R.begin(), neighbors_in_R.end(), 0);
    fill(set_L.begin(), set_L.end(), 0);
    fill(set_R.begin(), set_R.end(), 0);
    n_L = 0;
    n_R = 0;
    n_to_visit = 0;
}

void vide_decoder::add_to_R(int chk, boost::multi_array<int, 2> &chk_nbhd) {
    set_R[chk] = 1;
    n_R++;

    for (int j = 0; j < dc_q; j++) {
        int bit = chk_nbhd[chk][j];
        neighbors_in_R[bit]+=1;
        int h = h_bb;
        if (bit >= n*n) {
            h = h_cc;
        }
        if (neighbors_in_R[bit] > h && !set_L[bit]) {
            to_visit[n_to_visit] = bit;
            n_to_visit++;
            set_L[bit] = 1;
            L[n_L] = bit;
            n_L++;
        }
    }
}

void vide_decoder::init(std::vector<int> &syndrome, boost::multi_array<int, 2> &chk_nbhd) {
    for (int i = 0; i < n_chk_q; i++) {
        if (syndrome[i]) {
            add_to_R(i, chk_nbhd);
        }
    }
}

void vide_decoder::visit_last(boost::multi_array<int, 2> &chk_nbhd, boost::multi_array<int, 2> &bit_nbhd) {
    int bit = to_visit[n_to_visit-1];
    n_to_visit--;
    
    for (int i = 0; i < dv_q[bit]; i++) {
        int chk = bit_nbhd[bit][i];
        if (!set_R[chk]) {
            add_to_R(chk, chk_nbhd);
        }
    }
}

bool vide_decoder::decode_helper(boost::multi_array<int, 2> &bit_nbhd, boost::multi_array<int, 2> &chk_nbhd, vector<int> &syndrome, vector<int> &errs) {
    reset();
    init(syndrome, chk_nbhd);
    
    while (n_to_visit > 0) {
        visit_last(chk_nbhd, bit_nbhd);
    }
    
    return true;
}

bool vide_decoder::L_correctable(boost::multi_array<int, 2> &bit_nbhd, std::vector<int> &syndrome) {
    //implement
    return false;
}

bool vide_decoder::decode() {
    bool x_decoded = decode_helper(z_bit_nbhd, z_chk_nbhd, z_syndrome, x_errs);
    bool z_decoded = decode_helper(x_bit_nbhd, x_chk_nbhd, x_syndrome, z_errs);
    return x_decoded && z_decoded;
}

void vide_decoder::debug() {
    int trials = 1;
    float err = .002;
    vector<float> p_err {1.0f - 3*err, err, err, err};
    
    int missed_ct = 0;
    int full_ct = 0;
    int ct = 0;
    int lmf_ct = 0;
    for (int i = 0; i < trials; i++) {
//        clear();
//        add_errors(p_err);
//        add_error(n*264, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*257, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*21, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*239, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*209, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*13, x_bit_nbhd, x_syndrome, z_errs);
//        add_error(n*n+47, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*n+262, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*n+9, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*n+161, x_bit_nbhd, x_syndrome, z_errs);
        add_error(n*n+26, x_bit_nbhd, x_syndrome, z_errs);
//        add_error(112934, x_bit_nbhd, x_syndrome, z_errs);
//        add_error(112993, x_bit_nbhd, x_syndrome, z_errs);
//        add_error(113026, x_bit_nbhd, x_syndrome, z_errs);
//        add_error(113068, x_bit_nbhd, x_syndrome, z_errs);
        
        decode_helper(x_bit_nbhd, x_chk_nbhd, x_syndrome, z_errs);
        
        
        int missed_errors = 0;
        int L_minus_F = 0;
        for (int i = 0; i < n_q; i++) {
            if (z_errs[i] && ! set_L[i]) {
                missed_errors++;
                
                cout << "Qubit: " << i << "=" << bit_coord(i).first << "," << bit_coord(i).second << endl;
                cout << "Neighbors :";
                for (int j = 0; j < dv_q[i]; j++) {
                    int chk = x_bit_nbhd[i][j];

                    cout << "   " << chk_coord(chk).first << "," << chk_coord(chk).second;
                    if (x_syndrome[chk]) {
                        cout << "x";
                    }
                    if (set_R[chk]) {
                        cout << "*";
                    }
                }
                cout << endl;
            }
            
            if (set_L[i] && !z_errs[i]) {
                L_minus_F++;
            }
        }
        
        if (missed_errors > 0) {
            missed_ct++;
            cout << ct << endl;
//            for (int i = 0; i < dv_q[137691]; i++) {
//                int neighbor = x_bit_nbhd[137691][i];
//                cout << neighbor << "=" << chk_coord(neighbor).first << "," << chk_coord(neighbor).second << "-----" << endl;
//                cout << "     ";
//                for (int j = 0; j < dc_q; j++) {
//                    int bit2 = x_chk_nbhd[neighbor][j];
//                    if (z_errs[bit2]) {
//                        if (set_L[bit2]) {
//                            cout << "L";
//                        }
//                        cout << bit2 << "=" << bit_coord(bit2).first << "," << bit_coord(bit2).second << "  ";
//                    }
//                }
//                cout << endl;
//            }
            
            cout << "L=" << n_L << " , " << "E=" << (n_errs[0] + n_errs[1]) << ", E \\ L=" << missed_errors << endl;
            cout << "R=" << n_R << endl;
            cout << "===============" << endl;
        }
        
        ct++;
        
        if (n_L == n_q) {
            full_ct++;
//            cout << "full" << endl;
        }
        
        if (missed_errors > 0  && (n_errs[0]+n_errs[1] - missed_errors - L_minus_F > 0)) {
            lmf_ct++;
        }
    }
    
    cout << "MISSED: " << missed_ct << endl;
    cout << "FULL: " << full_ct << endl;
    cout << "LMF: " << lmf_ct << endl;
}
