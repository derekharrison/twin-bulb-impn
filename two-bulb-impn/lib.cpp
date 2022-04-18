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
#include "user_types.h"

std::vector<std::vector<double>> bulb1_t;
std::vector<std::vector<double>> bulb2_t;

double ** mat2D(int n) {

    double ** mat = new double * [n];

    for(int i = 0; i < n; ++i)
        mat[i] = new double[n];

    return mat;
}

void free_mat2D(double ** mat, int n) {

    for(int i = 0; i < n; ++i)
        delete [] mat[i];

    delete [] mat;
}

void get_order(double ** mat, int n, double * order_arr) {
    for(int row = 0; row < n; ++row) {
        int order = 0;
        for(int c = 0; c < n; ++c) {
            if(mat[row][c] == 0) {
                order++;
            }
        }
        order_arr[row] = order;
    }
}

void merge(oa_elem_t A[], int p, int q, int r) {
    int size_r, size_l;
    int i, j;
    size_l = q - p + 1;
    size_r = r - q;
    oa_elem_t L[size_l + 1];
    oa_elem_t R[size_r + 1];
    i = 0;
    j = 0;
    for(int n = p; n < q + 1; ++n) {
        L[i] = A[n];
        ++i;
    }
    L[size_l].val = MAX_INT;
    for(int n = q + 1; n < r + 1; ++n) {
        R[j] = A[n];
        ++j;
    }
    R[size_r].val = MAX_INT;
    i = 0;
    j = 0;
    for(int n = p; n < r + 1; ++n) {
        if(L[i].val < R[j].val) {
            A[n] = L[i];
            ++i;
        }
        else {
            A[n] = R[j];
            ++j;
        }
    }
}

void merge_sort(oa_elem_t A[], int p, int r) {
    int q;
    if(p < r) {
        q = (p + r)/2;
        merge_sort(A, p, q);
        merge_sort(A, q + 1, r);
        merge(A, p, q, r);
    }
}

void mergesort(oa_elem_t A[], int size) {
    merge_sort(A, 0, size - 1);
}



void make_ordered_mat(double ** mat, int n, double * order_arr, double ** ordered_mat) {

    oa_elem_t * order_array = new oa_elem_t[n];

    for(int row = 0; row < n; ++row) {
        order_array[row].old_row = row;
        order_array[row].val = order_arr[row];
    }

    mergesort(order_array, n);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            int old_row = order_array[row].old_row;
            ordered_mat[row][c] = mat[old_row][c];
        }
    }

    delete [] order_array;
}

int count_leading_zeros(double ** mat, int n, int row) {

    int count_lz = 0;

    for(int c = 0; c < n; ++c) {
        if(mat[row][c] == 0) {
            count_lz++;
        }
        else {
            break;
        }
    }

    return count_lz;
}

void mat_mult_sq(double ** A, double ** A_inv, int n, double ** mat_res) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0;

            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + A[i][k] * A_inv[k][j];
            }

            mat_res[i][j] = sum_loc;
        }
    }
}

void singularity_check(double ** mat_ref, int n, bool & mat_is_singular) {
    // Check for singularity
    for(int row = 0; row < n; ++row) {
        bool all_zeros_r = true;
        for(int c = 0; c < n; ++c) {
            if(mat_ref[row][c] != 0)
                all_zeros_r = false;
        }
        if(all_zeros_r && !mat_is_singular) {
            printf("Matrix is singular\n");
            mat_is_singular = true;
        }
    }

    for(int col = 0; col < n; ++col) {
        bool all_zeros_c = true;
        for(int row = 0; row < n; ++row) {
            if(mat_ref[row][col] != 0)
                all_zeros_c = false;
        }
        if(all_zeros_c && !mat_is_singular) {
            printf("Matrix is singular\n");
            mat_is_singular = true;
        }
    }
}


void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double ** mat_ref = mat2D(n);
    double ** mat_ordered = mat2D(n);
    double ** mat_inv_ordered = mat2D(n);
    double * order_arr = new double[n];


    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            if(c == row) {
                mat_inv[row][c] = 1.0;
            }
            else {
                mat_inv[row][c] = 0.0;
            }
        }
    }

    get_order(mat, n, order_arr);

    make_ordered_mat(mat, n, order_arr, mat_ordered);

    make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat_ref[row][c] = mat_ordered[row][c];
            mat_inv[row][c] = mat_inv_ordered[row][c];
        }
    }

    // Initialize singularity flag
    bool mat_is_singular = false;

    // Check if input matrix is singular
    singularity_check(mat_ref, n, mat_is_singular);

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat_ref[c][c]) < 1e-8) {
            get_order(mat_ref, n, order_arr);

            make_ordered_mat(mat_ref, n, order_arr, mat_ordered);

            make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

            for(int row = 0; row < n; ++row) {
                for(int c = 0; c < n; ++c) {
                    mat_ref[row][c] = mat_ordered[row][c];
                    mat_inv[row][c] = mat_inv_ordered[row][c];
                }
            }
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat_ref[c][col] = mat_ref[c][col] / (mat_ref[c][c] + 1e-13);
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = mat_inv[c][col] / (mat_ref[c][c] + 1e-13);
        }

        mat_ref[c][c] = 1.0;

        // Delete elements in rows below
        for(int row = c + 1; row < n; ++row) {
            if(mat_ref[row][c] != 0) {
                for(int col = c + 1; col < n; ++col) {
                    mat_ref[row][col] = -1.0 * mat_ref[row][c] * mat_ref[c][col] + mat_ref[row][col];
                }
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }

            int num_lead_zeros = count_leading_zeros(mat_ref, n, row);

            if(num_lead_zeros == n && !mat_is_singular) {
                printf("Matrix is singular\n");
                mat_is_singular = true;
            }
        }

    }

    // Backtracing
    for(int c = n - 1; c > 0; --c) {
        for(int row = c - 1; row > -1; row--) {
            if(mat_ref[row][c] != 0) {
                for(int col = 0; col < n; col++) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }
        }
    }

    // Check if matrix is singular
    singularity_check(mat_ref, n, mat_is_singular);

    // Free allocated space
    free_mat2D(mat_ref, n);
    free_mat2D(mat_ordered, n);
    free_mat2D(mat_inv_ordered, n);
    delete [] order_arr;
}

void init_diffusivities(p_params_t & p_params, int n) {
    p_params.D = mat2D(n);
    
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

double determinant(double ** mat, int n) {
    double det = 0;

    if(n == 1) {
        return mat[0][0];
    }

    if(n == 2) {
        return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    }

    if(n > 2) {
        double ** M = mat2D(n - 1);

        for(int c = 0; c < n; ++c) {
            for(int i = 1; i < n; ++i) {
                int j_m = 0;
                for(int j = 0; j < n; ++j) {
                    if(j != c) {
                        M[i - 1][j_m] = mat[i][j];
                        j_m++;
                    }

                }
            }

            double fac = pow(-1, c + 2);

            det = det + mat[0][c] * fac * determinant(M, n - 1);
        }

        free_mat2D(M, n - 1);
    }

    return det;
}


double co_factor(double ** mat, int n, int i, int j) {

    double fac = 0;

    double ** mat_red = mat2D(n - 1);

    int i_m = 0;
    for(int r = 0; r < n; ++r) {
        int j_m = 0;
        if(r != i) {
            for(int c = 0; c < n; ++c) {
                if(c != j) {
                    mat_red[i_m][j_m] = mat[r][c];
                    j_m++;
                }
            }
            i_m++;
        }
    }

    fac = pow(-1, i + j + 2) * determinant(mat_red, n - 1);

    free_mat2D(mat_red, n - 1);

    return fac;
}

void adjoint_mat(double ** mat, int n, double ** adj_mat) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            adj_mat[i][j] = co_factor(mat, n, i, j);
        }
    }
}

void compute_mat_inv(double ** mat, int n, double ** mat_inv) {

    double ** adj_mat = mat2D(n);

    adjoint_mat(mat, n, adj_mat);

    double det = determinant(mat, n);

    // Check if matrix mat is singular
    if(det == 0)
        std::cout << "matrix is singular" << std::endl;

    // Matrix mat is non-singular
    if(det != 0) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                mat_inv[j][i] = 1.0 / det * adj_mat[i][j];
            }
        }
    }

    free_mat2D(adj_mat, n);
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
    int n = comp_data.n - 0;
    
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
            double w_term = comp_data.tube_fracs[flux_node - 1].x[c];
            double e_term = comp_data.tube_fracs[flux_node].x[c];
            x_loc[c] = 0.5 * (w_term + e_term);
            sum_l_c = sum_l_c + x_loc[c];
        }
        
        x_loc[n] = 1.0 - sum_l_c;
        
        compute_coeff_matrix(flux_node, comp_data, x_loc);
    }
    
    // Bulb 2
    compute_coeff_matrix(ng, comp_data, comp_data.bulb_data_inter.mol_fracs_bulb2.x);
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
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz = (comp_data_N.bulb_data.mol_fracs_bulb2.x[j] - comp_data_N.tube_fracs[ng - 1].x[j]) / (0.5 * dz);
                sum_loc = sum_loc + A_inv[i][j] * dxj_dz;
            }
        }
        double term_w = A * A_inv[i][i] * comp_data_N.tube_fracs[ng - 1].x[i] / (0.5 * dz);
        double old_term = V / dt * comp_data_N.bulb_data_old.mol_fracs_bulb2.x[i];
        sum_loc = - A * sum_loc + term_w + old_term;

        comp_data_N.bulb_data.mol_fracs_bulb2.x[i] = sum_loc / (V / dt + A * A_inv[i][i] / (0.5 * dz));
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
        double term_w = A_inv_w[i][i] * comp_data_N.bulb_data_inter.mol_fracs_bulb1.x[i] / (0.5 * dz);
        double term_e = A_inv_e[i][i] * comp_data_N.tube_fracs[1].x[i] / dz;
        double old_term = dz / dt * comp_data_N.tube_fracs_old[0].x[i];
        double sum_loc_tot = - sum_loc_w + sum_loc_e + term_w + term_e + old_term;

        double ap = dz / dt + A_inv_w[i][i] / (0.5 * dz) + A_inv_e[i][i] / dz;

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

            comp_data_N.tube_fracs[node].x[i] = sum_loc_tot / (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / dz);
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

        comp_data_N.tube_fracs[ng - 1].x[i] = sum_loc_tot / (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / (0.5 * dz));
    }
}

void comp_data_ref(c_data_t & comp_data) {

    bulb1(comp_data);
    
    tube0(comp_data);
    
    tube_mid(comp_data);
    
    tube_n(comp_data);
    
    bulb2(comp_data);
    
    set_frac_comp3(comp_data);
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
    
    comp_data.bulb_data.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb2.x = new double[n];
    
    for(int c = 0; c < n; ++c) {
        comp_data.bulb_data.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
    }
    
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

    // Compute composition
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    // Loop to time t = tf
    while(t < t_params.tf) {
        
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

        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {
            
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
            
            reset_coefficients(comp_data);
            
            compute_coefficients(comp_data);
                 
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {
       
                // Update estimates bulb and tube composition
                comp_data_ref(comp_data);
                
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
