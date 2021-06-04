struct BasicInput {
    float *coefficients;
    unsigned short coefficients_length;
    float **constraint_matrix;
    unsigned short constraint_matrix_length;
    float **right_side_matrix;
    unsigned short right_side_matrix_length;
    float **bounds;
    unsigned short bounds_length;
};

struct RemasteredInput {
    float *coefficients;
    int coefficients_length;
    float **constraint_matrix;
    int constraint_matrix_length;
    float *right_side_array;
    int right_side_array_length;
};

struct FloatVector {
    float *values;
};

struct IntVector {
    int *values;
};

struct AnotherOneOutput {
    int vars_length;
    int constrs_length;
    float *obj;
    int obj_length;
    float **newA;
    int newA_length;
    float *sense;
    int sense_length;
    float *rhs;
    int rhs_length;
    float *lb;
    int lb_length;
};

struct CanonicalForm {
    float *obj_fun_coefficients;
    int obj_fun_coefficients_length;
    float **constrs_matrix;
    int constrs_matrix_rows;
    int constrs_matrix_columns;
    float *right_hand_side;
    int right_hand_side_length;
    float *low_bounds;
    int low_bounds_length;
};

struct ToPartialProblemInput {
    float **matrix;
    int matrix_rows;
    int matrix_columns;
    float *right_side_array;
    int right_side_array_length;
    float *vector;
    int vector_length;
};

struct ToPartialProblemOutput {
    float **constraint_matrix;
    int constraint_matrix_rows;
    int constraint_matrix_columns;
    float **right_side;
    int right_side_rows;
    int right_side_columns;
    int nz;
    int nx;
};

struct MPEC_solver_input {
    float *vector_x0;
    int vector_x0_length;
    float **matrixA;
    int matrixA_rows;
    int matrixA_columns;
    float **matrixB;
    int matrixB_rows;
    int matrixB_columns;
    float **matrixC;
    int matrixC_rows;
    int matrixC_columns;
    float *vector_b;
    int vector_b_length;
    float *vector_c;
    int vector_c_length;
};
struct MPEC_solver_output {
    float *objfun_coeffs;
    int objfun_coeffs_length;
    float **constrs_matrix;
    int constrs_matrix_rows;
    int constrs_matrix_columns;
    float *sense_array;
    int sense_array_length;
    float *right_hand_side;
    int right_hand_side_length;
    float **bounds;
    int bounds_rows;
    int bounds_columns;
    int nx;
    int nbin;
};

struct P_matrix {
    float **matrix;
    int matrix_side;
};

struct RemasteredInput remaster_basic_input(struct BasicInput input);
struct AnotherOneOutput function2(struct RemasteredInput input, struct FloatVector vector);
struct CanonicalForm to_canonical_form(struct RemasteredInput input);
struct ToPartialProblemOutput to_partial_problem(struct ToPartialProblemInput input);
struct MPEC_solver_output solve_inverse_via_MPEC(struct MPEC_solver_input input);
struct P_matrix get_P_matrix(int n, int m);
#ifndef INVERSE_PROGRAMMING_MAIN_H
#define INVERSE_PROGRAMMING_MAIN_H

#endif //INVERSE_PROGRAMMING_MAIN_H
