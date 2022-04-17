//
//  user_types.h
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#ifndef user_types_h
#define user_types_h

const int MAX_OUT = 100;
const int MAX_IN = 100;
const int MAX_INT = 1215752192;

typedef struct node_data {
    double * x;
    int n;
} node_t;

typedef struct coeff_mat_boundary {
    double ** A;
    double ** A_inv;
    int n;
} c_mat_b_t;

typedef struct physical_params {
    double ct;
    double ** D;
    int n;
    double D12;
    double D13;
    double D23;
} p_params_t;

typedef struct bulb_data {
    node_t mol_fracs_bulb1;
    node_t mol_fracs_bulb2;
} b_data_t;

typedef struct time_parameters {
    double to;
    double tf;
    double dt;
    int nt;
} t_params_t;

typedef struct experiment_params {
    double V;
    double d;
    double len;
    double A;
    double dz;
} e_params_t;

typedef struct computation_data {
    int ng;
    int n;

    c_mat_b_t * coeff_mat_boundary;
    
    p_params_t p_params;
    e_params_t e_params;
    t_params_t t_params;

    b_data_t bulb_data_inter;
    b_data_t bulb_data_old;
    b_data_t bulb_data;

    node_t * tube_fracs;
    node_t * tube_fracs_old;
    node_t * tube_fracs_inter;
} c_data_t;

typedef struct order_array_elem {
    int old_row;
    double val;
} oa_elem_t;

#endif /* user_types_h */
