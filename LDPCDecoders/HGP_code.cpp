//
//  HGP_code.cpp
//  LDPCDecoders
//
//  Created by Jon on 5/9/22.
//

#include "HGP_code.hpp"
#include <iostream>

using namespace std;

string HGP_code::read_param(ifstream& file, string label) {
    string tmp;
    if (!getline(file,tmp,',')) {
        return "-1";
    }
    
    if (tmp.compare(label)) {
        throw invalid_argument("Bad file format read codes");
    }
    
    getline(file,tmp,'\n');
    return tmp;
}

void HGP_code::read_matrix(ifstream &file, boost::multi_array<int, 2> &matrix, string label) {
    int n_rows = (int)matrix.shape()[0];
    int n_cols = (int)matrix.shape()[1];
    
    string tmp;
    getline(file,tmp,'\n');
    if (tmp.compare(label)) {
    throw invalid_argument("Bad file format read codes2");
    }

    for (int i = 0; i < n_rows; ++i){
    for (int j = 0; j < n_cols; ++j){
        getline(file,tmp,',');
        matrix[i][j] = stoi(tmp);
    }
    getline(file,tmp,'\n');
    }
}

void HGP_code::read_file(string filename) {
    ifstream file;
    file.open(filename);
    
    n = stoi(read_param(file,"n"));
    m = stoi(read_param(file,"m"));
    dv = stoi(read_param(file,"dv"));
    dc = stoi(read_param(file,"dc"));
    id = read_param(file,"id");

    bit_nbhd_c.resize(boost::extents[n][dv]);
    read_matrix(file, bit_nbhd_c, "bit_nbhd");
    chk_nbhd_c.resize(boost::extents[m][dc]);
    read_matrix(file, chk_nbhd_c, "check_nbhd");
    
    file.close();
}

void HGP_code::gen_bits() {
    dv_q.resize(n_q);
    int dv_max = max(dv,dc);
    x_bit_nbhd.resize(boost::extents[n_q][dv_max]);
    z_bit_nbhd.resize(boost::extents[n_q][dv_max]);
    
    for (int b1 = 0; b1 < n; b1++) {
        for (int b2 = 0; b2 < n; b2++) {
            int ind = b1*n+b2;
            for (int i = 0; i < dv; i++) {
                int chk_x = m*b1 + bit_nbhd_c[b2][i];
                int chk_z = m*b2 + bit_nbhd_c[b1][i];
                x_bit_nbhd[ind][i] = chk_x;
                z_bit_nbhd[ind][i] = chk_z;
            }
            dv_q[ind] = dv;
        }
    }
    
    for (int c1 = 0; c1 < m; c1++) {
        for (int c2 = 0; c2 < m; c2++) {
            int ind = n*n + c1*m+c2;
            for (int i = 0; i < dc; i++) {
                int chk_x = m*chk_nbhd_c[c1][i] +c2;
                int chk_z = m*chk_nbhd_c[c2][i] +c1;
                x_bit_nbhd[ind][i] = chk_x;
                z_bit_nbhd[ind][i] = chk_z;
            }
            dv_q[ind] = dc;
        }
    }
}

void HGP_code::gen_checks() {
    x_chk_nbhd.resize(boost::extents[n_chk_q][dc_q]);
    z_chk_nbhd.resize(boost::extents[n_chk_q][dc_q]);
    
    for (int b = 0; b < n; b++) {
        for (int c = 0; c < m; c++) {
            int ind = m*b+c;
            
            for (int i_b = 0; i_b < dc; i_b++) {
                int bit_x = n*b+chk_nbhd_c[c][i_b];
                int bit_z = n*chk_nbhd_c[c][i_b]+b;
                x_chk_nbhd[ind][i_b]=bit_x;
                z_chk_nbhd[ind][i_b]=bit_z;
            }
            
            for (int i_c = 0; i_c < dv; i_c++) {
                int bit_x = n*n + m*bit_nbhd_c[b][i_c] + c;
                int bit_z = n*n + m*c + bit_nbhd_c[b][i_c];
                x_chk_nbhd[ind][dc+i_c] = bit_x;
                z_chk_nbhd[ind][dc+i_c] = bit_z;
            }
        }
    }
}

void HGP_code::gen_q_code() {
    n_q = n*n+m*m;
    n_chk_q = n*m;
    dc_q = dc+dv;
    x_syndrome.resize(n_chk_q);
    z_syndrome.resize(n_chk_q);
    x_errs.resize(n_q);
    z_errs.resize(n_q);
    
    gen_bits();
    gen_checks();
}

HGP_code::HGP_code(string filename, mt19937& engine) :
engine(engine),
dis(0.0,1.0),
n_errs(3)
{
    read_file(filename);
    gen_q_code();
}

void HGP_code::clear() {
    fill(x_syndrome.begin(), x_syndrome.end(), 0);
    fill(z_syndrome.begin(), z_syndrome.end(), 0);
    fill(x_errs.begin(), x_errs.end(), 0);
    fill(z_errs.begin(), z_errs.end(), 0);
    fill(n_errs.begin(), n_errs.end(), 0);
}

void HGP_code::add_error(int b, boost::multi_array<int,2>& bit_nbhd, std::vector<int>& syndrome, std::vector<int>& errs) {
    errs[b] ^= 1;
    for (int i = 0; i < dv_q[b]; i++) {
        int chk = bit_nbhd[b][i];
        syndrome[chk] ^= 1;
    }
}

void HGP_code::add_errors(vector<float>& error_p) {
    float cumulative[4];
    cumulative[0] = error_p[0];
    for (int i = 1; i < 4; i++) {
        cumulative[i] = cumulative[i-1] + error_p[i];
    }
    
    for (int b = 0; b < n_q; b++) {
        float p = dis(engine);

        if (p >= cumulative[0] && p < cumulative[1]) {
            n_errs[0]++;
            add_error(b, z_bit_nbhd, z_syndrome, x_errs);
        } else if (p >= cumulative[1] && p < cumulative[2]) {
            n_errs[1]++;
            add_error(b, z_bit_nbhd, z_syndrome, x_errs);
            add_error(b, x_bit_nbhd, x_syndrome, z_errs);
        } else if (p >= cumulative[2]) {
            n_errs[2]++;
            add_error(b, x_bit_nbhd, x_syndrome, z_errs);
        }
    }
}

void HGP_code::get_n_errors(std::vector<int>& errors_out) {
    errors_out = n_errs;
}

bool HGP_code::check_correct_graph() {
    for (int b = 0; b < n_q; b++) {
        for (int i = 0; i < dv_q[b]; i++) {
            int x_neighbor = x_bit_nbhd[b][i];
            int z_neighbor = z_bit_nbhd[b][i];
            
            bool found_x_neighbor = false;
            for (int i2 = 0; i2 < dc_q; i2++) {
                if (x_chk_nbhd[x_neighbor][i2] == b) {
                    found_x_neighbor = true;
                }
            }
            
            bool found_z_neighbor = false;
            for (int i2 = 0; i2 < dc_q; i2++) {
                if (z_chk_nbhd[z_neighbor][i2] == b) {
                    found_z_neighbor = true;
                }
            }
            
            if (!found_x_neighbor || !found_z_neighbor) {
                return false;
            }
        }
    }
    
    for (int c = 0; c < n_chk_q; c++) {
        for (int i = 0; i < dc_q; i++) {
            int x_neighbor = x_chk_nbhd[c][i];
            int z_neighbor = z_chk_nbhd[c][i];
            
            bool found_x_neighbor = false;
            for (int i2 = 0; i2 < dv_q[x_neighbor]; i2++){
                if (x_bit_nbhd[x_neighbor][i2] == c) {
                    found_x_neighbor = true;
                }
            }
            
            bool found_z_neighbor = false;
            for (int i2 = 0; i2 < dv_q[z_neighbor]; i2++){
                if (z_bit_nbhd[z_neighbor][i2] == c) {
                    found_z_neighbor = true;
                }
            }
            
            if (!found_x_neighbor || !found_z_neighbor) {
                return false;
            }
        }
    }
    
    return true;
}

bool HGP_code::check_correct_stabilizers() {
    for (int chk_x = 0; chk_x < n_chk_q; chk_x++) {
        for (int chk_z = 0; chk_z < n_chk_q; chk_z++) {
            int parity = 0;
            for (int i = 0; i < dc_q; i++) {
                for (int j = 0; j < dc_q; j++) {
                    if (x_chk_nbhd[chk_x][i] == z_chk_nbhd[chk_z][j]) {
                        parity++;
                    }
                }
            }
            
            if (parity % 2 != 0) {
                return false;
            }
        }
    }
    
    return true;
}

bool HGP_code::check_correctness() {
    return check_correct_graph() && check_correct_stabilizers();
}

pair<int,int> HGP_code::chk_coord(int chk) {
    int b = chk / m;
    int c = chk % m;
    return pair<int,int>(b,c);
}

pair<int,int> HGP_code::bit_coord(int bit) {
    if (bit < n*n) {
        int b1 = bit / n;
        int b2 = bit % n;
        return pair<int,int>(b1,b2);
    } else {
        int c1 = (bit - n*n) / m;
        int c2 = (bit - n*n) % m;
        return pair<int,int>(c1,c2);
    }
}
