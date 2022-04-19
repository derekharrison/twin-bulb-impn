/*
 * lib_cofactor.cpp
 *
 *  Created on: Apr 20, 2022
 *      Author: d-w-h
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "mem_ops.hpp"

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

