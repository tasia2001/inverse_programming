#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "main.h"

const float inf = (float) 99999999.0;
const float undefined = (float) 67108864.0;
const float BIG_DIGIT = (float) 88888888.0;

typedef struct {
    float *array;
    size_t used;
    size_t size;
} Array;

typedef struct {
    float **array;
    size_t used;
    size_t size;
} PointerArray;

void initPointerArray(PointerArray *a, size_t initialSize) {
    a->array = malloc(initialSize * sizeof(float *));
    a->used = 0;
    a->size = initialSize;
}

void initArray(Array *a, size_t initialSize) {
    a->array = malloc(initialSize * sizeof(float));
    a->used = 0;
    a->size = initialSize;
}

void insertPointerArray(PointerArray *a, float *element) {
    if (a->used == a->size) {
        a->size *= 2;
        a->array = realloc(a->array, a->size * sizeof(float *));
    }
    a->array[a->used++] = element;
}

void insertArray(Array *a, float element) {
    if (a->used == a->size) {
        a->size *= 2;
        a->array = realloc(a->array, a->size * sizeof(float));
    }
    a->array[a->used++] = element;
}

void printPointerArray(PointerArray *a, unsigned short inner_array_length) {
    for (int i = 0; i < a->used; ++i) {
        for (int j = 0; j < inner_array_length; ++j) {
            printf("%f ", a->array[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printArray(Array *a) {
    for (int i = 0; i < a->used; ++i) {
        printf("%f ", a->array[i]);
    }
    printf("\n");
}

void freePointerArray(PointerArray *a) {
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

void freeArray(Array *a) {
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

void insert_positive_value(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array, int i) {
    insertPointerArray(constraint_matrix, input->constraint_matrix[i]);
    insertArray(right_side_array, input->right_side_matrix[i][0]);
}

void insert_negative_value(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array, int i) {
    float *constraint_matrix_row = malloc(input->coefficients_length * sizeof(float));
    for (int j = 0; j < input->coefficients_length; ++j) {
        constraint_matrix_row[j] = -input->constraint_matrix[i][j];
    }
    insertPointerArray(constraint_matrix, constraint_matrix_row);
    insertArray(right_side_array, -input->right_side_matrix[i][0]);
}

void remake_right_side_matrix(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array) {
    for (int i = 0; i < input->right_side_matrix_length; ++i) {
        if (input->right_side_matrix[i][1] == -1) { // <=
            insert_negative_value(input, constraint_matrix, right_side_array, i);
        } else if (input->right_side_matrix[i][1] == 0) { // ==
            insert_negative_value(input, constraint_matrix, right_side_array, i);
            insert_positive_value(input, constraint_matrix, right_side_array, i);
        } else if (input->right_side_matrix[i][1] == 1) { // >=
            insert_positive_value(input, constraint_matrix, right_side_array, i);
        }
    }
}

void insert_positive_bound(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array, int i) {
    float *constraint_matrix_row = malloc(input->coefficients_length * sizeof(float));
    for (int j = 0; j < input->coefficients_length; ++j) {
        if (i != j) {
            constraint_matrix_row[j] = 0;
        } else {
            constraint_matrix_row[j] = 1;
        }
    }
    insertPointerArray(constraint_matrix, constraint_matrix_row);
    insertArray(right_side_array, input->bounds[i][0]);
}

void insert_negative_bound(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array, int i) {
    float *constraint_matrix_row = malloc(input->coefficients_length * sizeof(float));
    for (int j = 0; j < input->coefficients_length; ++j) {
        if (i != j) {
            constraint_matrix_row[j] = 0;
        } else {
            constraint_matrix_row[j] = -1;
        }
    }
    insertPointerArray(constraint_matrix, constraint_matrix_row);
    insertArray(right_side_array, -input->bounds[i][1]);
}

void remake_bounds(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array) {
    for (int i = 0; i < input->bounds_length; ++i) {
        if (input->bounds[i][0] > -inf && input->bounds[i][1] < +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] <= -inf && input->bounds[i][1] < +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] > -inf && input->bounds[i][1] >= +inf) {
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        }
    }
}
struct RemasteredInput remaster_basic_input(struct BasicInput input) {
    PointerArray new_constraint_matrix;
    initPointerArray(&new_constraint_matrix, 1);

    Array new_right_side_array;
    initArray(&new_right_side_array, 1);
    remake_right_side_matrix(&input, &new_constraint_matrix, &new_right_side_array);
    remake_bounds(&input, &new_constraint_matrix, &new_right_side_array);

    struct RemasteredInput output;
    output.coefficients = input.coefficients;
    output.coefficients_length = input.coefficients_length;

    output.constraint_matrix = new_constraint_matrix.array;
    output.constraint_matrix_length = (int) new_constraint_matrix.used;

    output.right_side_array = new_right_side_array.array;
    output.right_side_array_length = (int) new_right_side_array.used;

    return output;
}

Array get_obj_vector(int units, int zeroes) {
    Array obj;
    initArray(&obj, 1);
    for (int i = 0; i < units; ++i) {
        insertArray(&obj, 1);
    }
    for (int i = 0; i < zeroes; ++i) {
        insertArray(&obj, 0);
    }
    return obj;
}

float *multiply_matrix_and_vector(float **matrix, int matrix_length, int coefficients_length,
                                  float *vector, int vector_length, float *return_vector) {
    for (int i = 0; i < vector_length; ++i) {
        return_vector[i] = 0;
    }
    for (int i = 0; i < matrix_length; ++i) {
        for (int j = 0; j < coefficients_length; ++j) {
            return_vector[i] += matrix[i][j] * vector[j];
        }
    }
    return return_vector;
}

PointerArray get_newA_matrix(float **constraint_matrix, int constraint_matrix_length,
                             float *right_side_array, int right_side_array_length,
                             int vars_length, Array B) {
    PointerArray newA_matrix;
    initPointerArray(&newA_matrix, 1);


    for (int i = 0; i < vars_length; ++i) {
        Array tmp_value;
        initArray(&tmp_value, 1);

        for (int j = 0; j < vars_length; ++j) {
            insertArray(&tmp_value, 0);
        }
        for (int j = 0; j < vars_length; ++j) {
            if (i == j) {
                insertArray(&tmp_value, -1);
            } else {
                insertArray(&tmp_value, 0);
            }
        }

        for (int j = 0; j < constraint_matrix_length; ++j) {
            bool inside_of_B = false;
            for (int k = 0; k < B.used; ++k) {
                if ((int) B.array[k] == j) {
                    insertArray(&tmp_value, constraint_matrix[j][i]);
                    inside_of_B = true;
                    break;
                }
            }
            if (!inside_of_B) {
                inside_of_B = false;
                insertArray(&tmp_value, 0);
            }
        }
        insertPointerArray(&newA_matrix, tmp_value.array);
    }


    for (int i = 0; i < vars_length; ++i) {
        Array tmp_value1, tmp_value2;
        initArray(&tmp_value1, 1);
        initArray(&tmp_value2, 1);
        for (int j = 0; j < vars_length; ++j) {
            if (i == j) {
                insertArray(&tmp_value1, -1);
            } else {
                insertArray(&tmp_value1, 0);
            }
        }

        for (int j = 0; j < vars_length; ++j) {
            if (i == j) {
                insertArray(&tmp_value1, 1);
            } else {
                insertArray(&tmp_value1, 0);
            }
        }

        for (int j = 0; j < constraint_matrix_length; ++j) {
            insertArray(&tmp_value1, 0);
        }
        insertPointerArray(&newA_matrix, tmp_value1.array);

        for (int j = 0; j < vars_length; ++j) {
            if (i == j) {
                insertArray(&tmp_value2, 1);
            } else {
                insertArray(&tmp_value2, 0);
            }
        }

        for (int j = 0; j < vars_length; ++j) {
            if (i == j) {
                insertArray(&tmp_value2, 1);
            } else {
                insertArray(&tmp_value2, 0);
            }
        }

        for (int j = 0; j < constraint_matrix_length; ++j) {
            insertArray(&tmp_value2, 0);
        }

        insertPointerArray(&newA_matrix, tmp_value2.array);

    }
    return newA_matrix;
}

Array get_sense_vector(int vars_length) {
    Array sense;
    initArray(&sense, 1);
    for (int i = 0; i < vars_length; ++i) {
        insertArray(&sense, 0);
    }
    for (int i = 0; i < vars_length * 2; i += 2) {
        insertArray(&sense, -1);
        insertArray(&sense, 1);
    }
    return sense;
}

Array get_rhs_vector(float *coefficients, int vars_length) {
    Array rhs;
    initArray(&rhs, 1);
    for (int i = 0; i < vars_length; ++i) {
        insertArray(&rhs, 0);
    }
    for (int i = 0; i < vars_length; ++i) {
        insertArray(&rhs, coefficients[i]);
        insertArray(&rhs, coefficients[i]);
    }
    return rhs;
}

Array get_B(float **constraint_matrix, int constraint_matrix_length, int coefficients_length,
            float *right_side_array, int right_side_array_length, float *vector) {
    Array B;
    initArray(&B, 1);

    float multiply_result[constraint_matrix_length];
    multiply_matrix_and_vector(constraint_matrix, constraint_matrix_length, coefficients_length,
                               vector, right_side_array_length, multiply_result);
    for (int i = 0; i < constraint_matrix_length; ++i) {
        if (multiply_result[i] == right_side_array[i]) {
            insertArray(&B, (float) i);
        }
    }
    return B;
}

Array get_lb_vector(int vars_length, Array B, int matrix_length) {
    Array lb;
    initArray(&lb, 1);
    for (int i = 0; i < vars_length; ++i) {
        insertArray(&lb, 0);
    }

    for (int i = 0; i < vars_length; ++i) {
        insertArray(&lb, -inf);
    }

    for (int i = 0; i < matrix_length; ++i) {
        bool inside_of_B = false;
        for (int j = 0; j < B.used; ++j) {
            if (i == (int) B.array[j]) {
                insertArray(&lb, 0);
                inside_of_B = true;
                break;
            }
        }
        if (!inside_of_B) {
            insertArray(&lb, -inf);
            inside_of_B = false;
        }
    }
    return lb;
}

struct AnotherOneOutput function2(struct RemasteredInput input, struct FloatVector vector) {
    struct AnotherOneOutput output;
    output.vars_length = input.coefficients_length;
    output.constrs_length = input.constraint_matrix_length;

    Array obj = get_obj_vector(input.coefficients_length, input.constraint_matrix_length + input.coefficients_length);
    output.obj = obj.array;
    output.obj_length = 2 * input.coefficients_length + input.constraint_matrix_length;

    Array B = get_B(input.constraint_matrix, input.constraint_matrix_length, input.coefficients_length,
                    input.right_side_array, input.right_side_array_length, vector.values);
    PointerArray newA = get_newA_matrix(input.constraint_matrix, input.constraint_matrix_length, input.right_side_array,
                                        input.right_side_array_length, input.coefficients_length, B);
    output.newA = newA.array;
    output.newA_length = 2 * input.coefficients_length + input.constraint_matrix_length;

    Array sense = get_sense_vector(input.coefficients_length);
    output.sense = sense.array;
    output.sense_length = 3 * input.coefficients_length;

    Array rhs = get_rhs_vector(input.coefficients, input.coefficients_length);
    output.rhs = rhs.array;
    output.rhs_length = input.coefficients_length * 3;

    Array lb = get_lb_vector(input.coefficients_length, B, input.constraint_matrix_length);
    output.lb = lb.array;
    output.lb_length = 2 * input.coefficients_length + input.constraint_matrix_length;
    return output;
}

Array get_obj_fun_coefficients(float *coefficients, int coefficients_length, int constrs_matrix_length) {
    Array obj_fun_coefficients;
    initArray(&obj_fun_coefficients, 1);
    for (int i = 0; i < coefficients_length; ++i) {
        insertArray(&obj_fun_coefficients, coefficients[i]);
        insertArray(&obj_fun_coefficients, -coefficients[i]);
    }
    for (int i = 0; i < constrs_matrix_length; ++i) {
        insertArray(&obj_fun_coefficients, 0);
    }
    return obj_fun_coefficients;
}

PointerArray get_constrs_matrix(float **constraint_matrix, int constraint_matrix_rows, int constraint_matrix_columns) {
    PointerArray constrs_matrix;
    initPointerArray(&constrs_matrix, 1);

    for (int i = 0; i < constraint_matrix_rows; ++i) {
        Array constrs_matrix_row;
        initArray(&constrs_matrix_row, 1);
        for (int j = 0; j < constraint_matrix_columns; ++j) {
            insertArray(&constrs_matrix_row, constraint_matrix[i][j]);
            insertArray(&constrs_matrix_row, -constraint_matrix[i][j]);
        }

        for (int j = 0; j < constraint_matrix_rows; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row, -1);
            } else {
                insertArray(&constrs_matrix_row, 0);
            }
        }
        insertPointerArray(&constrs_matrix, constrs_matrix_row.array);
    }

    return constrs_matrix;
}

Array get_low_bounds(int vars, int constrs) {
    Array low_bounds;
    initArray(&low_bounds, 1);
    for (int i = 0; i < 2 * vars; ++i) {
        insertArray(&low_bounds, 0);
    }
    for (int i = 0; i < constrs; ++i) {
        insertArray(&low_bounds, 0);
    }
    return low_bounds;
}

struct CanonicalForm to_canonical_form(struct RemasteredInput input) {
    struct CanonicalForm output;

    Array obj_fun_coefficients = get_obj_fun_coefficients(input.coefficients, input.coefficients_length,
                                                          input.constraint_matrix_length);
    output.obj_fun_coefficients = obj_fun_coefficients.array;
    output.obj_fun_coefficients_length = (int) obj_fun_coefficients.used;

    PointerArray constrs_matrix = get_constrs_matrix(input.constraint_matrix, input.constraint_matrix_length,
                                                     input.coefficients_length);
    output.constrs_matrix = constrs_matrix.array;
    output.constrs_matrix_rows = (int) constrs_matrix.used;
    output.constrs_matrix_columns = (int) input.coefficients_length * 2 + input.constraint_matrix_length;

    output.right_hand_side = input.right_side_array;
    output.right_hand_side_length = input.right_side_array_length;

    Array low_bounds = get_low_bounds(input.coefficients_length, input.constraint_matrix_length);
    output.low_bounds = low_bounds.array;
    output.low_bounds_length = (int) low_bounds.used;

    return output;
}

float *get_zeroes_line(int length, float *line) {
    for (int i = 0; i < length; ++i) {
        line[i] = 0;
    }
    return line;
}

float *get_undefined_positions(float *vector, int vector_length) {
    Array positions;
    initArray(&positions, 1);
    for (int i = 0; i < vector_length; ++i) {
        if (vector[i] == undefined) {
            insertArray(&positions, (float) i);
        }
    }
    return positions.array;
}

int get_undefined_count(float *vector, int vector_length) {
    int count = 0;
    for (int i = 0; i < vector_length; ++i) {
        if (vector[i] == undefined) {
            ++count;
        }
    }
    return count;
}

float *get_defined_positions(float *vector, int vector_length) {
    Array positions;
    initArray(&positions, 1);
    for (int i = 0; i < vector_length; ++i) {
        if (vector[i] != undefined && vector[i] != 0.0) {
            insertArray(&positions, (float) i);
        }
    }
    return positions.array;
}

int get_defined_count(float *vector, int vector_length) {
    int count = 0;
    for (int i = 0; i < vector_length; ++i) {
        if (vector[i] != undefined && vector[i] != 0.0) {
            ++count;
        }
    }
    return count;
}

PointerArray get_constraint_matrix(float **matrix, int matrix_columns, int matrix_rows,
                                   float *undefined_positions, int undefined_positions_count,
                                   float *defined_positions, int defined_positions_count) {
    PointerArray constraint_matrix;
    initPointerArray(&constraint_matrix, 1);
    int c0_length = 3 * matrix_columns + matrix_rows;
    float *c0 = NULL;
    // 1
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[matrix_columns + (int) undefined_positions[i]] = 1;
        c0[2 * matrix_columns + matrix_rows + (int) undefined_positions[i]] = -1 / BIG_DIGIT;
        insertPointerArray(&constraint_matrix, c0);
    }
    // 2
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[matrix_columns + (int) undefined_positions[i]] = 1;
        c0[2 * matrix_columns + matrix_rows + (int) undefined_positions[i]] = -BIG_DIGIT;
        insertPointerArray(&constraint_matrix, c0);
    }
    // 3
    for (int i = 0; i < defined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[(int) defined_positions[i]] = -1;
        for (int j = 2 * matrix_columns; j < 2 * matrix_columns + matrix_rows; ++j) {
            c0[j] = matrix[j - 2 * matrix_columns][(int) defined_positions[i]];
        }
        insertPointerArray(&constraint_matrix, c0);
    }
    // 4
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[(int) undefined_positions[i]] = -1;
        for (int j = 2 * matrix_columns; j < 2 * matrix_columns + matrix_rows; ++j) {
            c0[j] = matrix[j - 2 * matrix_columns][(int) undefined_positions[i]];
        }
        c0[2 * matrix_columns + matrix_rows + (int) undefined_positions[i]] = BIG_DIGIT;
        insertPointerArray(&constraint_matrix, c0);
    }
    // 5
    for (int i = 0; i < matrix_columns; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[i] = -1;
        for (int j = 2 * matrix_columns; j < 2 * matrix_columns + matrix_rows; ++j) {
            c0[j] = matrix[j - 2 * matrix_columns][i];
        }
        insertPointerArray(&constraint_matrix, c0);
    }
    // 6
    for (int i = 0; i < matrix_rows; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        for (int j = 0; j <= undefined_positions_count; ++j) {
            c0[matrix_columns + (int) undefined_positions[j]] = matrix[i][(int) undefined_positions[j]];
        }
        insertPointerArray(&constraint_matrix, c0);
    }
    // 7
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[matrix_columns + (int) undefined_positions[i]] = 1;
        insertPointerArray(&constraint_matrix, c0);
    }
    // 8
    for (int i = 0; i < defined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[matrix_columns + (int) defined_positions[i]] = 1;
        insertPointerArray(&constraint_matrix, c0);
    }

    return constraint_matrix;
}

PointerArray get_right_side(float **matrix, int matrix_columns, int matrix_rows,
                            float *undefined_positions, int undefined_positions_count,
                            float *defined_positions, int defined_positions_count,
                            float *vector, float *right_side_array) {
    PointerArray right_side;
    initPointerArray(&right_side, 1);
    int c0_length = 2;
    float *c0 = NULL;
    // 1
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = 0;
        c0[1] = 1;
        insertPointerArray(&right_side, c0);
    }
    // 2
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = 0;
        c0[1] = -1;
        insertPointerArray(&right_side, c0);
    }
    // 3
    for (int i = 0; i < defined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = 0;
        c0[1] = 0;
        insertPointerArray(&right_side, c0);
    }
    // 4
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = BIG_DIGIT;
        c0[1] = -1;
        insertPointerArray(&right_side, c0);
    }
    // 5
    for (int i = 0; i < matrix_rows; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = 0;
        c0[1] = 1;
        insertPointerArray(&right_side, c0);
    }
    // 6
    for (int i = 0; i < matrix_columns; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        float k = 0;
        for (int j = 0; j < defined_positions_count; ++j) {
            k += matrix[j][i] * vector[(int) defined_positions[j]];
        }
        c0[0] = right_side_array[i] - k;
        c0[1] = 1;
        insertPointerArray(&right_side, c0);
    }
    // 7
    for (int i = 0; i < undefined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = 0;
        c0[1] = 1;
        insertPointerArray(&right_side, c0);
    }
    // 8
    for (int i = 0; i < defined_positions_count; ++i) {
        c0 = malloc(sizeof(float) * c0_length);
        get_zeroes_line(c0_length, c0);
        c0[0] = vector[(int) defined_positions[i]];
        c0[1] = 0;
        insertPointerArray(&right_side, c0);
    }
    return right_side;
}

struct ToPartialProblemOutput to_partial_problem(struct ToPartialProblemInput input) {
    struct ToPartialProblemOutput output;
    int undefined_positions_count = get_undefined_count(input.vector, input.matrix_columns);
    float *undefined_positions = get_undefined_positions(input.vector, input.matrix_columns);
    int defined_positions_count = get_defined_count(input.vector, input.matrix_columns);
    float *defined_positions = get_defined_positions(input.vector, input.matrix_columns);
    PointerArray constraint_matrix;
    initPointerArray(&constraint_matrix, 1);
    constraint_matrix = get_constraint_matrix(input.matrix, input.matrix_columns, input.matrix_rows,
                                              undefined_positions, undefined_positions_count,
                                              defined_positions, defined_positions_count);

    PointerArray right_side;
    initPointerArray(&right_side, 1);
    right_side = get_right_side(input.matrix, input.matrix_columns, input.matrix_rows,
                                undefined_positions, undefined_positions_count,
                                defined_positions, defined_positions_count,
                                input.vector, input.right_side_array);

    output.constraint_matrix = constraint_matrix.array;
    output.constraint_matrix_rows = (int) constraint_matrix.used;
    output.constraint_matrix_columns = 3 * input.matrix_columns + input.matrix_rows;

    output.right_side = right_side.array;
    output.right_side_rows = (int) right_side.used;
    output.right_side_columns = 2;

    output.nz = input.vector_length;
    output.nx = output.constraint_matrix_columns - output.nz;

    return output;
}

Array get_objfun_coeffs(int vars, int constrA) {
    Array objfun_coeffs;
    initArray(&objfun_coeffs, 4 * vars + 2 * constrA);
    for (int i = 0; i < 2 * vars + 2 * constrA; ++i) {
        insertArray(&objfun_coeffs, 0);
    }
    for (int i = 0; i < vars; ++i) {
        insertArray(&objfun_coeffs, 1);
    }
    for (int i = 0; i < vars; ++i) {
        insertArray(&objfun_coeffs, 0);
    }
    return objfun_coeffs;
}

PointerArray get_constrs_matrix_MPEC(float **matrixA, int matrixA_rows, int matrixA_columns,
                                     float **matrixB, int matrixB_rows, int matrixB_columns,
                                     float **matrixC, int matrixC_rows, int matrixC_columns) {
    PointerArray constrs_matrix;
    initPointerArray(&constrs_matrix, 1);

    // 1
    for (int i = 0; i < matrixA_rows; ++i) {
        Array constrs_matrix_row;
        initArray(&constrs_matrix_row, 1);
        for (int j = 0; j < matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row, matrixA[i][j]);
        }
        for (int j = 0; j < matrixA_rows; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row, -1);
            } else {
                insertArray(&constrs_matrix_row, 0);
            }
        }
        for (int j = 0; j < 3 * matrixA_columns + matrixA_rows; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }
        insertPointerArray(&constrs_matrix, constrs_matrix_row.array);
    }

    // 2
    for (int i = 0; i < matrixA_columns; ++i) {
        Array constrs_matrix_row1;
        Array constrs_matrix_row2;
        initArray(&constrs_matrix_row1, 1);
        initArray(&constrs_matrix_row2, 1);
        for (int j = 0; j < matrixA_rows + matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row1, 0);
            insertArray(&constrs_matrix_row2, 0);
        }

        for (int j = 0; j < matrixA_rows; ++j) {
            insertArray(&constrs_matrix_row1, matrixA[j][i]);
            insertArray(&constrs_matrix_row2, matrixA[j][i]);
        }

        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row1, -1);
                insertArray(&constrs_matrix_row2, -1);
            } else {
                insertArray(&constrs_matrix_row1, 0);
                insertArray(&constrs_matrix_row2, 0);
            }
        }

        // constrs_matrix_row1
        for (int j = 0; j < 2 * matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row1, 0);
        }

        // constrs_matrix_row2
        for (int j = 0; j < matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row2, 0);
        }

        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row2, inf);
            } else {
                insertArray(&constrs_matrix_row2, 0);
            }
        }

        insertPointerArray(&constrs_matrix, constrs_matrix_row1.array);
        insertPointerArray(&constrs_matrix, constrs_matrix_row2.array);
    }

    // 3
    for (int i = 0; i < matrixA_columns; ++i) {
        Array constrs_matrix_row;
        initArray(&constrs_matrix_row, 1);
        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row, 1);
            } else {
                insertArray(&constrs_matrix_row, 0);
            }
        }

        for (int j = 0; j < 2 * matrixA_rows + 2 * matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }

        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row, -inf);
            } else {
                insertArray(&constrs_matrix_row, 0);
            }
        }
        insertPointerArray(&constrs_matrix, constrs_matrix_row.array);
    }

    // 4
    for (int i = 0; i < matrixB_rows; ++i) {
        Array constrs_matrix_row;
        initArray(&constrs_matrix_row, 1);

        for (int j = 0; j < matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }

        for (int j = 0; j < matrixB_columns; ++j) {
            insertArray(&constrs_matrix_row, matrixB[i][j]);
        }

        for (int j = 0; j < 3 * matrixA_columns + matrixA_rows; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }

        insertPointerArray(&constrs_matrix, constrs_matrix_row.array);
    }

    // 5
    for (int i = 0; i < matrixC_rows; ++i) {
        Array constrs_matrix_row;
        initArray(&constrs_matrix_row, 1);

        for (int j = 0; j < 2 * matrixA_rows + matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }

        for (int j = 0; j < matrixC_columns; ++j) {
            insertArray(&constrs_matrix_row, matrixC[i][j]);
        }

        for (int j = 0; j < 2 * matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row, 0);
        }

        insertPointerArray(&constrs_matrix, constrs_matrix_row.array);
    }

    // 6
    for (int i = 0; i < matrixA_columns; ++i) {
        Array constrs_matrix_row1;
        Array constrs_matrix_row2;
        initArray(&constrs_matrix_row1, 1);
        initArray(&constrs_matrix_row2, 1);

        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row1, 1);
                insertArray(&constrs_matrix_row2, 1);
            } else {
                insertArray(&constrs_matrix_row1, 0);
                insertArray(&constrs_matrix_row2, 0);
            }
        }

        for (int j = 0; j < 2 * matrixA_rows + matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row1, 0);
            insertArray(&constrs_matrix_row2, 0);
        }

        for (int j = 0; j < matrixA_columns; ++j) {
            if (i == j) {
                insertArray(&constrs_matrix_row1, -1);
                insertArray(&constrs_matrix_row2, 1);
            } else {
                insertArray(&constrs_matrix_row1, 0);
                insertArray(&constrs_matrix_row2, 0);
            }
        }

        for (int j = 0; j < matrixA_columns; ++j) {
            insertArray(&constrs_matrix_row1, 0);
            insertArray(&constrs_matrix_row2, 0);
        }

        insertPointerArray(&constrs_matrix, constrs_matrix_row1.array);
        insertPointerArray(&constrs_matrix, constrs_matrix_row2.array);
    }

    return constrs_matrix;
}

Array get_sense_array(int matrixA_rows, int matrixA_columns, int matrixB_rows, int matrixC_rows) {
    Array sense_array;
    initArray(&sense_array, 1);
    for (int i = 0; i < matrixA_rows; ++i) {
        insertArray(&sense_array, 0);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        insertArray(&sense_array, 1);
        insertArray(&sense_array, -1);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        insertArray(&sense_array, -1);
    }

    for (int i = 0; i < matrixB_rows; ++i) {
        insertArray(&sense_array, 1);
    }

    for (int i = 0; i < matrixC_rows; ++i) {
        insertArray(&sense_array, 1);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        insertArray(&sense_array, 1);
        insertArray(&sense_array, -1);
    }

    return sense_array;
}

Array get_right_hand_side(int matrixA_rows, int matrixA_columns,
                          float *vector_x0, int vector_x0_length,
                          float *vector_b, int vector_b_length,
                          float *vector_c, int vector_c_length) {
    Array right_hand_side;
    initArray(&right_hand_side, 1);

    for (int i = 0; i < matrixA_rows + matrixA_columns; ++i) {
        insertArray(&right_hand_side, 0);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        insertArray(&right_hand_side, inf);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        insertArray(&right_hand_side, 0);
    }

    for (int i = 0; i < vector_b_length; ++i) {
        insertArray(&right_hand_side, vector_b[i]);
    }

    for (int i = 0; i < vector_c_length; ++i) {
        insertArray(&right_hand_side, vector_c[i]);
    }

    for (int i = 0; i < vector_x0_length; ++i) {
        insertArray(&right_hand_side, vector_x0[i]);
        insertArray(&right_hand_side, vector_x0[i]);
    }

    return right_hand_side;
}

PointerArray get_bounds(int matrixA_rows, int matrixA_columns) {
    PointerArray bounds;
    initPointerArray(&bounds, 1);

    for (int i = 0; i < matrixA_columns; ++i) {
        Array bounds_row;
        initArray(&bounds_row, 1);
        insertArray(&bounds_row, 0);
        insertArray(&bounds_row, inf);
        insertPointerArray(&bounds, bounds_row.array);
    }

    for (int i = 0; i < 2 * matrixA_rows + matrixA_columns; ++i) {
        Array bounds_row;
        initArray(&bounds_row, 1);
        insertArray(&bounds_row, -inf);
        insertArray(&bounds_row, inf);
        insertPointerArray(&bounds, bounds_row.array);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        Array bounds_row;
        initArray(&bounds_row, 1);
        insertArray(&bounds_row, 0);
        insertArray(&bounds_row, inf);
        insertPointerArray(&bounds, bounds_row.array);
    }

    for (int i = 0; i < matrixA_columns; ++i) {
        Array bounds_row;
        initArray(&bounds_row, 1);
        insertArray(&bounds_row, 0);
        insertArray(&bounds_row, 1);
        insertPointerArray(&bounds, bounds_row.array);
    }

    return bounds;
}

struct MPEC_solver_output solve_inverse_via_MPEC(struct MPEC_solver_input input) {
    struct MPEC_solver_output output;
    Array objfun_coeffs;
    objfun_coeffs = get_objfun_coeffs(input.vector_x0_length, input.matrixA_rows);

    output.objfun_coeffs = objfun_coeffs.array;
    output.objfun_coeffs_length = (int) objfun_coeffs.used;

    PointerArray constrs_matrix;
    initPointerArray(&constrs_matrix, 1);
    constrs_matrix = get_constrs_matrix_MPEC(input.matrixA, input.matrixA_rows, input.matrixA_columns,
                                             input.matrixB, input.matrixB_rows, input.matrixB_columns,
                                             input.matrixC, input.matrixC_rows, input.matrixC_columns);

    output.constrs_matrix = constrs_matrix.array;
    output.constrs_matrix_rows = (int) constrs_matrix.used;
    output.constrs_matrix_columns = 2 * input.matrixA_rows + 4 * input.matrixA_columns;

    Array sense_array;
    initArray(&sense_array, 1);
    sense_array = get_sense_array(input.matrixA_rows, input.matrixA_columns, input.matrixB_rows, input.matrixC_rows);

    output.sense_array = sense_array.array;
    output.sense_array_length = (int) sense_array.used;

    Array right_hand_side;
    initArray(&right_hand_side, 1);
    right_hand_side = get_right_hand_side(input.matrixA_rows, input.matrixA_columns,
                                          input.vector_x0, input.vector_x0_length,
                                          input.vector_b, input.vector_b_length,
                                          input.vector_c, input.vector_c_length);
    output.right_hand_side = right_hand_side.array;
    output.right_hand_side_length = (int) right_hand_side.used;

    PointerArray bounds;
    initPointerArray(&bounds, 1);
    bounds = get_bounds(input.matrixA_rows, input.matrixA_columns);
    output.bounds = bounds.array;
    output.bounds_rows = 2 * input.matrixA_rows + 4 * input.matrixA_columns;
    output.bounds_columns = 2;

    output.nbin = input.vector_x0_length;
    output.nx = output.constrs_matrix_columns - output.nbin;

    return output;
}

struct P_matrix get_P_matrix(int n, int m) {
    struct P_matrix p;
    PointerArray matrix;
    initPointerArray(&matrix, 1);
    for (int i = 0; i < 3 * n + m; ++i) {
        Array row;
        initArray(&row, 1);
        for (int j = 0; j < 3 * n + m; ++j) {
            insertArray(&row, 0);
        }
        insertPointerArray(&matrix, row.array);
    }
    for (int i = 0; i < n; ++i) {
        matrix.array[i][i + n] = 1;
    }
    p.matrix = matrix.array;
    p.matrix_side = 3 * n + m;
    return p;
}