/*
 * mat_ops.cpp
 *
 *  Created on: Apr 19, 2022
 *      Author: d-w-h
 */

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
