//
//  lib.cpp
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "lib.hpp"
#include "lib_cofactor.hpp"
#include "lib_gauss.hpp"
#include "user_types.h"

std::vector<std::vector<double>> bulb1_t;
std::vector<std::vector<double>> bulb2_t;

void init_diffusivities(p_params_t & p_params, int n) {
    
    for(int j = 0; j < n; ++j) {
        for(int c = 0; c < n; ++c) {
            p_params.D[j][c] = 0.0;
        }
    }
}

void set_frac_comp3(c_data_t & comp_data) {
    int ng = comp_data.ng;
    int n = comp_data.n;
    
    for(int node = 0; node < ng; ++node) {
        double sum_loc = 0.0;
        for(int c = 0; c < n - 1; ++c) {
            sum_loc = sum_loc + comp_data.tube_fracs[node].x[c];
        }
        comp_data.tube_fracs[node].x[n - 1] = 1.0 - sum_loc;
    }
    
    double sum_loc1 = 0.0;
    for(int c = 0; c < n - 1; ++c) {
        sum_loc1 = sum_loc1 + comp_data.bulb_data.mol_fracs_bulb1.x[c];
    }
    
    double sum_loc2 = 0.0;
    for(int c = 0; c < n - 1; ++c) {
        sum_loc2 = sum_loc2 + comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }
    
    comp_data.bulb_data.mol_fracs_bulb1.x[n - 1] = 1.0 - sum_loc1;
    comp_data.bulb_data.mol_fracs_bulb2.x[n - 1] = 1.0 - sum_loc2;
}

void compute_mat_inv1(double ** mat, int n, double ** mat_inv) {

    compute_mat_inv(mat, n, mat_inv);
}

void compute_mat_inv2(double ** mat, int n, double ** mat_inv) {

    gauss_jordan(mat, n, mat_inv);
}

void mat_prod(double ** mat1, double ** mat2, int n, double ** prod) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0.0;
            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + mat1[i][k] * mat2[k][j];
            }
            prod[i][j] = sum_loc;
        }
    }
}

void print_mat(double ** mat, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }
}

void compute_coeff_matrix(int flux_node,
                          c_data_t & comp_data,
                          double * x) {
    int n = comp_data.n;
    
    for(int i = 0; i < n - 1; ++i) {
        for(int j = 0; j < n - 1; ++j) {
            if(i == j) {
                double sum_loc = x[i] / (comp_data.p_params.D[i][n - 1] + 1e-13);
                double sum_frac = 0.0;
                for(int k = 0; k < n - 1; ++k) {
                    sum_frac = sum_frac + x[k];
                }
                x[n - 1] = 1.0 - sum_frac;
                for(int k = 0; k < n; ++k) {
                    if(k != i) {
                        sum_loc = sum_loc + x[k] / (comp_data.p_params.D[i][k] + 1e-13);
                    }
                }
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = sum_loc;
            }
            if(i != j) {
                double Dij = comp_data.p_params.D[i][j];
                double Dnm1 = comp_data.p_params.D[i][n - 1];
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = -x[i] * (1.0 / (Dij + 1e-13) - 1.0 / (Dnm1 + 1e-13));
            }
        }
    }
    
    compute_mat_inv2(comp_data.coeff_mat_boundary[flux_node].A,
                    n - 1,
                    comp_data.coeff_mat_boundary[flux_node].A_inv);
}

void compute_coefficients(c_data_t & comp_data) {
    int flux_node = 0;
    int ng = comp_data.ng;
    int n = comp_data.n - 1;
    
    // Bulb 1
    compute_coeff_matrix(flux_node, comp_data, comp_data.bulb_data_inter.mol_fracs_bulb1.x);
    
    // Tube
    for(flux_node = 1; flux_node < ng; ++flux_node) {
        double * x_loc = new double [n + 1];
        double sum_l_c = 0.0;
        for(int c = 0; c < n; ++c) {
            double w_term = comp_data.tube_fracs_inter[flux_node - 1].x[c];
            double e_term = comp_data.tube_fracs_inter[flux_node].x[c];
            x_loc[c] = 0.5 * (w_term + e_term);
            sum_l_c = sum_l_c + x_loc[c];
        }
        
        x_loc[n] = 1.0 - sum_l_c;
        
        compute_coeff_matrix(flux_node, comp_data, x_loc);
    }
    
    flux_node = ng;

    // Bulb 2
    compute_coeff_matrix(flux_node, comp_data, comp_data.bulb_data_inter.mol_fracs_bulb2.x);
}

void reset_coefficients(c_data_t & comp_data) {
    int ng = comp_data.ng;
    int n = comp_data.n;
    
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = 0.0;
                comp_data.coeff_mat_boundary[flux_node].A_inv[i][j] = 0.0;
            }
        }
        
    }
}

void bulb1(c_data_t & comp_data_N) {

    int flux_node = 0;
    int n = comp_data_N.n - 1;
    double V = comp_data_N.e_params.V;
    double A = comp_data_N.e_params.A;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv = comp_data_N.coeff_mat_boundary[flux_node].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc = 0.0;

        for(int c = 0; c < n; ++c) {
            if(c != i) {
                double dxj = comp_data_N.tube_fracs[0].x[c] - comp_data_N.bulb_data.mol_fracs_bulb1.x[c];
                double dxj_dz = dxj / (0.5 * dz);
                sum_loc = sum_loc + A * A_inv[i][c] * dxj_dz;
            }
        }

        sum_loc = sum_loc + A * A_inv[i][i] * comp_data_N.tube_fracs[0].x[i] / (0.5 * dz);
        sum_loc = sum_loc + V / dt * comp_data_N.bulb_data_old.mol_fracs_bulb1.x[i];

        double ap = V / dt + A * A_inv[i][i] / (0.5 * dz);

        comp_data_N.bulb_data.mol_fracs_bulb1.x[i] = sum_loc / ap;
    }
}

void bulb2(c_data_t & comp_data_N) {
    int ng = comp_data_N.ng;
    
    int flux_node = ng;
    int n = comp_data_N.n - 1;
    double V = comp_data_N.e_params.V;
    double A = comp_data_N.e_params.A;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv = comp_data_N.coeff_mat_boundary[flux_node].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc = 0.0;

        for(int c = 0; c < n; ++c) {
            if(i != c) {
                double dxj = comp_data_N.bulb_data.mol_fracs_bulb2.x[c] - comp_data_N.tube_fracs[ng - 1].x[c];
                double dxj_dz =  dxj / (0.5 * dz);
                sum_loc = sum_loc +  A * A_inv[i][c] * dxj_dz;
            }
        }

        sum_loc = -sum_loc + A * A_inv[i][i] * comp_data_N.tube_fracs[ng - 1].x[i] / (0.5 * dz);
        sum_loc = sum_loc + V / dt * comp_data_N.bulb_data_old.mol_fracs_bulb2.x[i];

        double ap = V / dt + A * A_inv[i][i] / (0.5 * dz);

        comp_data_N.bulb_data.mol_fracs_bulb2.x[i] = sum_loc / ap;
    }
}

void tube0(c_data_t & comp_data_N) {
    int flux_node_w = 0;
    int flux_node_e = 1;
    int n = comp_data_N.n - 1;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
    double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc_w = 0.0;
        double sum_loc_e = 0.0;

        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz_w = (comp_data_N.tube_fracs[0].x[j] - comp_data_N.bulb_data_inter.mol_fracs_bulb1.x[j]) / (0.5 * dz);
                sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                double dxj_dz_e = (comp_data_N.tube_fracs[1].x[j] - comp_data_N.tube_fracs[0].x[j]) / dz;
                sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
            }
        }

        double xiW = comp_data_N.bulb_data_inter.mol_fracs_bulb1.x[i];
        double xiE = comp_data_N.tube_fracs[1].x[i];
        double old_term = dz / dt * comp_data_N.tube_fracs_old[0].x[i];

        double sum_loc_tot = -sum_loc_w + sum_loc_e + A_inv_w[i][i] * xiW / (0.5 * dz) + A_inv_e[i][i] * xiE / dz + old_term;

        double ap = (dz / dt + A_inv_w[i][i] / (0.5 * dz) + A_inv_e[i][i] / dz);

        comp_data_N.tube_fracs[0].x[i] = sum_loc_tot / ap;
    }
}

void tube_mid(c_data_t & comp_data_N) {
    
    int ng = comp_data_N.ng;
    
    for(int node = 1; node < ng - 1; ++node) {
        int flux_node_w = node;
        int flux_node_e = node + 1;
        int n = comp_data_N.n - 1;
        double dz = comp_data_N.e_params.dz;
        double dt = comp_data_N.t_params.dt;

        double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
        double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;

        for(int i = 0; i < n; ++i) {
            double sum_loc_w = 0.0;
            double sum_loc_e = 0.0;

            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    double dxj_dz_w = (comp_data_N.tube_fracs[node].x[j] - comp_data_N.tube_fracs[node - 1].x[j]) / dz;
                    sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                    double dxj_dz_e = (comp_data_N.tube_fracs[node + 1].x[j] - comp_data_N.tube_fracs[node].x[j]) / dz;
                    sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
                }
            }
            double xiW = comp_data_N.tube_fracs[node - 1].x[i];
            double xiE = comp_data_N.tube_fracs[node + 1].x[i];
            double old_term = dz / dt * comp_data_N.tube_fracs_old[node].x[i];

            double sum_loc_tot = -sum_loc_w + sum_loc_e + A_inv_w[i][i] / dz * xiW + A_inv_e[i][i] / dz * xiE + old_term;

            double ap = (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / dz);

            comp_data_N.tube_fracs[node].x[i] = sum_loc_tot / ap;
        }
    }
}

void tube_n(c_data_t & comp_data_N) {
    int ng = comp_data_N.ng;
    
    int flux_node_w = ng - 1;
    int flux_node_e = ng;
    int n = comp_data_N.n - 1;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;

    double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
    double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc_w = 0.0;
        double sum_loc_e = 0.0;

        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz_w = (comp_data_N.tube_fracs[ng - 1].x[j] - comp_data_N.tube_fracs[ng - 2].x[j]) / dz;
                sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                double dxj_dz_e = (comp_data_N.bulb_data.mol_fracs_bulb2.x[j] - comp_data_N.tube_fracs[ng - 1].x[j]) / (0.5 * dz);
                sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
            }
        }

        double xiW = comp_data_N.tube_fracs[ng - 2].x[i];
        double xiE = comp_data_N.bulb_data.mol_fracs_bulb2.x[i];
        double old_term = dz / dt * comp_data_N.tube_fracs_old[ng - 1].x[i];

        double sum_loc_tot = -sum_loc_w + sum_loc_e + A_inv_w[i][i] / dz * xiW + A_inv_e[i][i] / (0.5 * dz) * xiE + old_term;

        double ap = (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / (0.5 * dz));

        comp_data_N.tube_fracs[ng - 1].x[i] = sum_loc_tot / ap;
    }
}


void update_compositions(c_data_t & comp_data) {

    bulb1(comp_data);
    
    tube0(comp_data);
    
    tube_mid(comp_data);
    
    tube_n(comp_data);
    
    bulb2(comp_data);
    
    set_frac_comp3(comp_data);
}

void allocate_composition_data(c_data_t & comp_data) {
    int n = comp_data.n;
    int ng = comp_data.ng;

    comp_data.bulb_data.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb2.x = new double[n];

    // Allocate data for tube composition
    comp_data.tube_fracs = new node_t[ng];
    comp_data.tube_fracs_inter = new node_t[ng];
    comp_data.tube_fracs_old = new node_t[ng];

    // Allocate data for flux coefficient matrices
    comp_data.coeff_mat_boundary = new c_mat_b_t[ng + 1];

    // Allocate space for and initialize coefficient matrix
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        comp_data.coeff_mat_boundary[flux_node].A = mat2D(n);
        comp_data.coeff_mat_boundary[flux_node].A_inv = mat2D(n);
        for(int c = 0; c < n; ++c) {
            for(int c_l = 0; c_l < n; ++c_l) {
                comp_data.coeff_mat_boundary[flux_node].A[c][c_l] = 0.0;
                comp_data.coeff_mat_boundary[flux_node].A_inv[c][c_l] = 0.0;
            }
        }
    }

    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x = new double[n];
        comp_data.tube_fracs_inter[node].x = new double[n];
        comp_data.tube_fracs_old[node].x = new double[n];
        for(int c = 0; c < n; ++c) {
            comp_data.tube_fracs[node].x[c] = 1.0 / n;
            comp_data.tube_fracs_inter[node].x[c] = 1.0 / n;
            comp_data.tube_fracs_old[node].x[c] = 1.0 / n;
        }
    }
}

void deallocate_composition_data(c_data_t & comp_data) {
    int n = comp_data.n;
    int ng = comp_data.ng;

    // Free allocated space
    delete [] comp_data.bulb_data.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data.mol_fracs_bulb2.x;
    delete [] comp_data.bulb_data_old.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data_old.mol_fracs_bulb2.x;
    delete [] comp_data.bulb_data_inter.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data_inter.mol_fracs_bulb2.x;

    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        free_mat2D(comp_data.coeff_mat_boundary[flux_node].A, n);
        free_mat2D(comp_data.coeff_mat_boundary[flux_node].A_inv, n);
    }

    delete [] comp_data.coeff_mat_boundary;

    for(int node = 0; node < ng; ++node) {
        delete [] comp_data.tube_fracs[node].x;
        delete [] comp_data.tube_fracs_inter[node].x;
        delete [] comp_data.tube_fracs_old[node].x;
    }

    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

void init_bulb_composition(b_data_t & bulb_data, c_data_t & comp_data) {
    int n = comp_data.n;
    for(int c = 0; c < n; ++c) {
        comp_data.bulb_data.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
    }
}

void update_old_timestep_data(c_data_t & comp_data) {
    int n = comp_data.n;
    int ng = comp_data.ng;

    // Update bulb composition of previous time step
    for(int c = 0; c < n; ++c) {
        comp_data.bulb_data_old.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }

    // Update tube composition of previous time step
    for(int node = 0; node < ng; ++node) {
        for(int c = 0; c < n; ++c) {
            comp_data.tube_fracs_old[node].x[c] = comp_data.tube_fracs[node].x[c];
        }
    }
}

void update_intermediate_data(c_data_t & comp_data) {
    int n = comp_data.n;
    int ng = comp_data.ng;

    // Update intermediate bulb composition
    for(int c = 0; c < n; ++c) {
        comp_data.bulb_data_inter.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }

    // Update intermediate tube composition
    for(int node = 0; node < ng; ++node) {
        for(int c = 0; c < n; ++c) {
            comp_data.tube_fracs_inter[node].x[c] = comp_data.tube_fracs[node].x[c];
        }
    }
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               int n,
                               b_data_t & bulb_data) {
    // Organize data
    c_data_t comp_data;
    comp_data.ng = ng;
    comp_data.n = n;

    comp_data.p_params = p_params;
    comp_data.e_params = e_params;
    comp_data.t_params = t_params;

    comp_data.bulb_data = bulb_data;
    comp_data.bulb_data_old = bulb_data;
    comp_data.bulb_data_inter = bulb_data;
    
    allocate_composition_data(comp_data);

    init_bulb_composition(bulb_data, comp_data);

    // Compute composition
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    // Loop to time t = tf
    while(t < t_params.tf) {
        
        update_old_timestep_data(comp_data);

        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {
            
            update_intermediate_data(comp_data);
            
            reset_coefficients(comp_data);
            
            compute_coefficients(comp_data);
                 
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {
       
                // Update estimates bulb and tube composition
                update_compositions(comp_data);
                
                in_it++;
            }
            
            out_it++;
        }

        std::vector<double> bulb1_t_loc;
        std::vector<double> bulb2_t_loc;

        for(int c = 0; c < n; ++c) {
            bulb1_t_loc.push_back(comp_data.bulb_data.mol_fracs_bulb1.x[c]);
            bulb2_t_loc.push_back(comp_data.bulb_data.mol_fracs_bulb2.x[c]);
        }

        bulb1_t.push_back(bulb1_t_loc);
        bulb2_t.push_back(bulb2_t_loc);
        
        t = t + dt;
    }
    
    // Store data in files
    store_fractions(t_params, n);

    // Set bulb data
    for(int c = 0; c < n; ++c) {
        bulb_data.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
        bulb_data.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }
    
    // Free allocated space
    deallocate_composition_data(comp_data);
}

void store_fractions(t_params_t t_params, int n) {
    int nt = (int) bulb1_t.size();
    double dt = (t_params.tf - t_params.to) / nt;

    FILE * file_ptr1 = fopen("results_bulb1.txt", "w");
    FILE * file_ptr2 = fopen("results_bulb2.txt", "w");

    if(file_ptr1 == nullptr)
        printf("file1 could not be opened\n");

    if(file_ptr2 == nullptr)
        printf("file2 could not be opened\n");

    double t = t_params.to;

    for(int i = 0; i < nt; ++i) {
        t = t + dt;
        fprintf(file_ptr1, "%f\t", t);
        fprintf(file_ptr2, "%f\t", t);
        for(int c = 0; c < n; ++c) {
            fprintf(file_ptr1, "%f\t", bulb1_t[i][c]);
            fprintf(file_ptr2, "%f\t", bulb2_t[i][c]);
        }
        fprintf(file_ptr1, "\n");
        fprintf(file_ptr2, "\n");
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
}
