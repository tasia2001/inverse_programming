#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

const float inf = 99999999.0;

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
<<<<<<< HEAD
        if (input->bounds[i][0] > -inf && input->bounds[i][1] < +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] <= -inf && input->bounds[i][1] < +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] > -inf && input->bounds[i][1] >= +inf) {
=======
        if (input->bounds[i][0] != -inf && input->bounds[i][1] != +inf) {
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] != -inf) {
>>>>>>> 743fafc (operation signs changed)
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][1] != +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
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
    for (int i = 0; i < vars_length * 2; ++i) {
        insertArray(&lb, 0);
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