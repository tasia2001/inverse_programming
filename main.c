#include <stdlib.h>
#include <stdio.h>

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

struct BasicOutput {
    float *coefficients;
    int coefficients_length;
    float **constraint_matrix;
    int constraint_matrix_length;
    float *right_side_array;
    int right_side_array_length;
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
            insert_positive_value(input, constraint_matrix, right_side_array, i);
        } else if (input->right_side_matrix[i][1] == 0) { // ==
            insert_positive_value(input, constraint_matrix, right_side_array, i);
            insert_negative_value(input, constraint_matrix, right_side_array, i);
        } else if (input->right_side_matrix[i][1] == 1) { // >=
            insert_negative_value(input, constraint_matrix, right_side_array, i);
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
    insertArray(right_side_array, input->bounds[i][1]);
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
    insertArray(right_side_array, -input->bounds[i][0]);
}

void remake_bounds(struct BasicInput *input, PointerArray *constraint_matrix, Array *right_side_array) {
    for (int i = 0; i < input->bounds_length; ++i) {
        if (input->bounds[i][0] != -inf && input->bounds[i][1] != +inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][0] != -inf) {
            insert_negative_bound(input, constraint_matrix, right_side_array, i);
        } else if (input->bounds[i][1] != +inf) {
            insert_positive_bound(input, constraint_matrix, right_side_array, i);
        }
    }
}

struct BasicOutput remake_basic_input(struct BasicInput input) {
    PointerArray new_constraint_matrix;
    initPointerArray(&new_constraint_matrix, 1);

    Array new_right_side_array;
    initArray(&new_right_side_array, 1);

    remake_right_side_matrix(&input, &new_constraint_matrix, &new_right_side_array);
    remake_bounds(&input, &new_constraint_matrix, &new_right_side_array);

    struct BasicOutput output;
    output.coefficients = input.coefficients;
    output.coefficients_length = input.coefficients_length;

    output.constraint_matrix = new_constraint_matrix.array;
    output.constraint_matrix_length = (int) new_constraint_matrix.used;

    output.right_side_array = new_right_side_array.array;
    output.right_side_array_length = (int) new_right_side_array.used;

    return output;
}
