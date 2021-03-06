/*
 * lib_gauss.cpp
 *
 *  Created on: Apr 19, 2022
 *      Author: dwh
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "lib_sort.hpp"
#include "mem_ops.hpp"
#include "user_types.h"

void get_order(double ** mat, int n, double * order_arr) {
    for(int row = 0; row < n; ++row) {
        int order = 0;
        while(fabs(mat[row][order]) <= SMALL_NUM && order < n) {
            mat[row][order] = 0.0;
            order++;
        }
        order_arr[row] = order;
    }
}

void sort_mat(double ** mat, int n, double * order_arr, double ** ordered_mat) {

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
            if(fabs(ordered_mat[row][c]) <= SMALL_NUM) {
                ordered_mat[row][c] = 0.0;
            }
        }
    }

    delete [] order_array;
}

int count_leading_zeros(double ** mat, int n, int row) {

    int count = 0;

    while(fabs(mat[row][count]) <= SMALL_NUM && count < n) {
        count++;
    }

    return count;
}

void singularity_check(double ** mat_ref, int n, bool & is_singular) {

    for(int row = 0; row < n; ++row) {
        bool all_zeros_c = true;
        for(int col = 0; col < n; ++col) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM) {
                all_zeros_c = false;
            }
        }
        if(all_zeros_c && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }

    for(int col = 0; col < n; ++col) {
        bool all_zeros_r = true;
        for(int row = 0; row < n; ++row) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM) {
                all_zeros_r = false;
            }
        }
        if(all_zeros_r && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }
}

void cut_low_vals(double ** mat, int n) {

    for(int row = 0; row < n; ++row) {
        for(int col = 0; col < n; ++col) {
            if(fabs(mat[row][col]) <= SMALL_NUM) {
                mat[row][col] = 0.0;
            }
        }
    }
}

void sort_matrix(double * order_arr, int n, double ** mat) {

    double ** mat_ordered = mat2D(n);

    sort_mat(mat, n, order_arr, mat_ordered);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat[row][c] = mat_ordered[row][c];
        }
    }

    free_mat2D(mat_ordered, n);
}

void init_mat_inv(double ** mat_inv, int n) {
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
}

void check_leading_zeros(double ** mat, int n, bool & is_singular) {
    // Check if matrix is singular
    for(int row = 0; row < n; ++row) {
        int num_lead_zeros = count_leading_zeros(mat, n, row);

        if(num_lead_zeros >= row + 1 && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }
}

void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double * order_arr = new double[n];

    // Initialize matrix inverse
    init_mat_inv(mat_inv, n);

    // Initialize singularity flag
    bool is_singular = false;

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat[c][c]) <= SMALL_NUM) {
            get_order(mat, n, order_arr);

            sort_matrix(order_arr, n, mat);

            sort_matrix(order_arr, n, mat_inv);

            check_leading_zeros(mat, n, is_singular);
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat[c][col] = mat[c][col] / (mat[c][c] + SMALL_NUM);
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = mat_inv[c][col] / (mat[c][c] + SMALL_NUM);
        }

        mat[c][c] = 1.0;

        // Delete elements in rows below
        for(int row = c + 1; row < n; ++row) {
            if(mat[row][c] != 0) {
                for(int col = c + 1; col < n; ++col) {
                    mat[row][col] = -1.0 * mat[row][c] * mat[c][col] + mat[row][col];
                }
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat[row][c] = 0;
            }
        }
    }

    // Backtrace to convert to reduced row echelon form
    for(int c = n - 1; c > 0; --c) {
        for(int row = c - 1; row > -1; --row) {
            if(mat[row][c] != 0) {
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat[row][c] = 0;
            }
        }
    }

    // Free allocated space
    delete [] order_arr;
}


