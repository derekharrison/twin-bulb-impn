//
//  lib.hpp
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>

#include "user_types.h"

double ** mat2D(int n);

void free_mat2D(double ** mat, int n);

void store_fractions(t_params_t t_params, int n);

void init_diffusivities(p_params_t & p_params, int n);

void compute_bulb_compositions(e_params_t e_params,
                                p_params_t p_params,
                                t_params_t t_params,
                                int ng,
                                int n,
                                b_data_t & bulb_data);
#endif /* lib_hpp */
