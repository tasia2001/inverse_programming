#include <python3.8/Python.h>
#include "main.h"

float *build_1D_c_array(PyListObject *array, int array_length) {
    float *c_array = (float *) malloc(sizeof(float) * array_length);
    for (int i = 0; i < array_length; ++i) {
        PyObject *item;
        item = PyList_GetItem(array, i);
        if (!PyFloat_Check(item))
            c_array[i] = 0;
        c_array[i] = (float) PyFloat_AsDouble(item);
    }
    return c_array;
}

float **build_2D_c_array(PyListObject *array, int array_length, int inner_array_length) {
    float **c_array = (float **) malloc(sizeof(float *) * array_length);
    for (int i = 0; i < array_length; ++i) {
        c_array[i] = (float *) malloc(sizeof(float) * inner_array_length);
    }
    for (int i = 0; i < array_length; ++i) {
        for (int j = 0; j < inner_array_length; ++j) {
            PyListObject *item_list;
            PyObject *item;
            item_list = PyList_GetItem(array, i);
            item = PyList_GetItem(item_list, j);
            if (!PyFloat_Check(item)) {
                c_array[i][j] = 0;
            }
            c_array[i][j] = (float) PyFloat_AsDouble(item);
        }
    }
    return c_array;
}

PyObject *build_1D_python_list(float *array, int array_length) {
    PyObject *py_array, *item;
    py_array = PyList_New(array_length);
    for (int i = 0; i < array_length; ++i) {
        item = PyFloat_FromDouble(array[i]);
        PyList_SetItem(py_array, i, item);
    }
    return py_array;
}

PyObject *build_2D_python_list(float **array, int array_length, int inner_array_length) {
    PyObject *py_array, *inner_py_array, *item;
    py_array = PyList_New(array_length);
    for (int i = 0; i < array_length; ++i) {
        inner_py_array = PyList_New(inner_array_length);
        for (int j = 0; j < inner_array_length; ++j) {
            item = PyFloat_FromDouble(array[i][j]);
            PyList_SetItem(inner_py_array, j, item);
        }
        PyList_SetItem(py_array, i, inner_py_array);
    }
    return py_array;
}

static PyObject *to_standard_form(PyObject *self, PyObject *args) {
    Py_ssize_t coefficients_length, constraint_matrix_length, constraint_matrix_width, right_side_matrix_length,
            right_side_matrix_width, bounds_length, bounds_width;
    PyListObject *coefficients, *constraint_matrix, *right_side_matrix, *bounds;
    PyObject *temp_value;
    if (!PyArg_ParseTuple(args, "OOOO", &coefficients, &constraint_matrix, &right_side_matrix, &bounds)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return NULL;
    }

    if (!PyList_Check(coefficients)) {
        PyErr_Format(PyExc_TypeError, "The coefficients must be of list or subtype of list");
        return NULL;
    }
    coefficients_length = PyList_Size(coefficients);
    if (coefficients_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_coefficients = build_1D_c_array(coefficients, coefficients_length);

    if (!PyList_Check(constraint_matrix)) {
        PyErr_Format(PyExc_TypeError, "The constraint_matrix must be of list or subtype of list");
        return NULL;
    }
    constraint_matrix_length = PyList_Size(constraint_matrix);
    if (constraint_matrix_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(constraint_matrix, 0);
    constraint_matrix_width = PyList_Size(temp_value);
    if (constraint_matrix_width <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_constraint_matrix = build_2D_c_array(constraint_matrix, constraint_matrix_length,
                                                   constraint_matrix_width);

    if (!PyList_Check(right_side_matrix)) {
        PyErr_Format(PyExc_TypeError, "The right_side_matrix must be of list or subtype of list");
        return NULL;
    }
    right_side_matrix_length = PyList_Size(right_side_matrix);
    if (right_side_matrix_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(right_side_matrix, 0);
    right_side_matrix_width = PyList_Size(temp_value);
    if (right_side_matrix_width <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_right_side_matrix = build_2D_c_array(right_side_matrix, right_side_matrix_length,
                                                   right_side_matrix_width);
    if (!PyList_Check(bounds)) {
        PyErr_Format(PyExc_TypeError, "The bounds must be of list or subtype of list");
        return NULL;
    }
    bounds_length = PyList_Size(bounds);
    if (bounds_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(bounds, 0);
    bounds_width = PyList_Size(temp_value);
    if (bounds_width <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_bounds = build_2D_c_array(bounds, bounds_length, bounds_width);

    struct BasicInput input;
    input.coefficients = c_coefficients;
    input.coefficients_length = (unsigned short) coefficients_length;
    input.constraint_matrix = c_constraint_matrix;
    input.constraint_matrix_length = (unsigned short) constraint_matrix_length;
    input.right_side_matrix = c_right_side_matrix;
    input.right_side_matrix_length = (unsigned short) right_side_matrix_length;
    input.bounds = c_bounds;
    input.bounds_length = (unsigned short) bounds_length;
    struct RemasteredInput output = remaster_basic_input(input);

    PyObject *py_coefficients = build_1D_python_list(output.coefficients, output.coefficients_length);
    PyObject *py_constraint_matrix = build_2D_python_list(output.constraint_matrix, output.constraint_matrix_length,
                                                          output.coefficients_length);
    PyObject *py_right_side_array = build_1D_python_list(output.right_side_array, output.right_side_array_length);

    return Py_BuildValue("(OOO)", py_coefficients, py_constraint_matrix, py_right_side_array);
}

static PyObject *get_inverse_cost_vector_model(PyObject *self, PyObject *args) {
    Py_ssize_t coefficients_length, constraint_matrix_length, constraint_matrix_width, right_side_array_length,
            values_length;
    PyListObject *coefficients, *constraint_matrix, *right_side_array, *values;
    PyObject *temp_value;
    if (!PyArg_ParseTuple(args, "OOOO", &coefficients, &constraint_matrix, &right_side_array, &values)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return NULL;
    }

    if (!PyList_Check(coefficients)) {
        PyErr_Format(PyExc_TypeError, "The coefficients must be of list or subtype of list");
        return NULL;
    }
    coefficients_length = PyList_Size(coefficients);
    if (coefficients_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_coefficients = build_1D_c_array(coefficients, coefficients_length);

    if (!PyList_Check(constraint_matrix)) {
        PyErr_Format(PyExc_TypeError, "The constraint_matrix must be of list or subtype of list");
        return NULL;
    }
    constraint_matrix_length = PyList_Size(constraint_matrix);
    if (constraint_matrix_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(constraint_matrix, 0);
    constraint_matrix_width = PyList_Size(temp_value);
    if (constraint_matrix_width <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_constraint_matrix = build_2D_c_array(constraint_matrix, constraint_matrix_length,
                                                   constraint_matrix_width);

    if (!PyList_Check(right_side_array)) {
        PyErr_Format(PyExc_TypeError, "The right_side_array must be of list or subtype of list");
        return NULL;
    }
    right_side_array_length = PyList_Size(right_side_array);
    if (right_side_array_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_right_side_array = build_1D_c_array(right_side_array, right_side_array_length);

    if (!PyList_Check(values)) {
        PyErr_Format(PyExc_TypeError, "The values must be of list or subtype of list");
        return NULL;
    }
    values_length = PyList_Size(values);
    if (values_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_values = build_1D_c_array(values, values_length);


    struct RemasteredInput input;
    input.coefficients = c_coefficients;
    input.coefficients_length = coefficients_length;
    input.constraint_matrix = c_constraint_matrix;
    input.constraint_matrix_length = constraint_matrix_length;
    input.right_side_array = c_right_side_array;
    input.right_side_array_length = right_side_array_length;
    struct FloatVector vector;
    vector.values = c_values;
    struct AnotherOneOutput output = function2(input, vector);

    PyObject *py_obj = build_1D_python_list(output.obj, output.obj_length);
    PyObject *py_newA = build_2D_python_list(output.newA, 3 * output.vars_length,
                                             output.newA_length);
    PyObject *py_sense = build_1D_python_list(output.sense, output.sense_length);
    PyObject *py_rhs = build_1D_python_list(output.rhs, output.rhs_length);
    PyObject *py_lb = build_1D_python_list(output.lb, output.lb_length);

    return Py_BuildValue("(iiOOOOO)", output.vars_length, output.constrs_length,
                         py_obj, py_newA, py_sense, py_rhs, py_lb);
}

static PyObject *module_to_canonical_form(PyObject *self, PyObject *args) {
    Py_ssize_t coefficients_length, constraint_matrix_length, constraint_matrix_width, right_side_array_length;
    PyListObject *coefficients, *constraint_matrix, *right_side_array;
    PyObject *temp_value;
    if (!PyArg_ParseTuple(args, "OOO", &coefficients, &constraint_matrix, &right_side_array)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return NULL;
    }

    if (!PyList_Check(coefficients)) {
        PyErr_Format(PyExc_TypeError, "The coefficients must be of list or subtype of list");
        return NULL;
    }
    coefficients_length = PyList_Size(coefficients);
    if (coefficients_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_coefficients = build_1D_c_array(coefficients, coefficients_length);

    if (!PyList_Check(constraint_matrix)) {
        PyErr_Format(PyExc_TypeError, "The constraint_matrix must be of list or subtype of list");
        return NULL;
    }
    constraint_matrix_length = PyList_Size(constraint_matrix);
    if (constraint_matrix_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(constraint_matrix, 0);
    constraint_matrix_width = PyList_Size(temp_value);
    if (constraint_matrix_width <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_constraint_matrix = build_2D_c_array(constraint_matrix, constraint_matrix_length,
                                                   constraint_matrix_width);

    if (!PyList_Check(right_side_array)) {
        PyErr_Format(PyExc_TypeError, "The right_side_array must be of list or subtype of list");
        return NULL;
    }
    right_side_array_length = PyList_Size(right_side_array);
    if (right_side_array_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_right_side_array = build_1D_c_array(right_side_array, right_side_array_length);


    struct RemasteredInput input;
    input.coefficients = c_coefficients;
    input.coefficients_length = coefficients_length;
    input.constraint_matrix = c_constraint_matrix;
    input.constraint_matrix_length = constraint_matrix_length;
    input.right_side_array = c_right_side_array;
    input.right_side_array_length = right_side_array_length;
    struct CanonicalForm output = to_canonical_form(input);

    PyObject *py_obj_fun_coefficients = build_1D_python_list(output.obj_fun_coefficients,
                                                             output.obj_fun_coefficients_length);
    PyObject *py_constrs_matrix = build_2D_python_list(output.constrs_matrix, output.constrs_matrix_rows,
                                                       output.constrs_matrix_columns);
    PyObject *py_right_hand_side = build_1D_python_list(output.right_hand_side, output.right_hand_side_length);
    PyObject *py_low_bounds = build_1D_python_list(output.low_bounds, output.low_bounds_length);

    return Py_BuildValue("(OOOO)", py_obj_fun_coefficients, py_constrs_matrix, py_right_hand_side, py_low_bounds);
}

static PyObject *get_partial_problem_model(PyObject *self, PyObject *args) {
    Py_ssize_t matrix_rows, matrix_columns, vector_length, right_side_array_length;
    PyListObject *matrix, *vector, *right_side_array;
    PyObject *temp_value;
    if (!PyArg_ParseTuple(args, "OOO", &matrix, &right_side_array, &vector)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return NULL;
    }

    if (!PyList_Check(matrix)) {
        PyErr_Format(PyExc_TypeError, "The matrix must be of list or subtype of list");
        return NULL;
    }
    matrix_rows = PyList_Size(matrix);
    if (matrix_rows <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(matrix, 0);
    matrix_columns = PyList_Size(temp_value);
    if (matrix_columns <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_matrix = build_2D_c_array(matrix, matrix_rows, matrix_columns);

    if (!PyList_Check(right_side_array)) {
        PyErr_Format(PyExc_TypeError, "The right_side_array must be of list or subtype of list");
        return NULL;
    }
    right_side_array_length = PyList_Size(right_side_array);
    if (right_side_array_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_right_side_array = build_1D_c_array(right_side_array, right_side_array_length);

    if (!PyList_Check(vector)) {
        PyErr_Format(PyExc_TypeError, "The vector must be of list or subtype of list");
        return NULL;
    }
    vector_length = PyList_Size(vector);
    if (vector_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_vector = build_1D_c_array(vector, vector_length);


    struct ToPartialProblemInput input;
    input.matrix = c_matrix;
    input.matrix_rows = matrix_rows;
    input.matrix_columns = matrix_columns;
    input.right_side_array = c_right_side_array;
    input.right_side_array_length = right_side_array_length;
    input.vector = c_vector;
    input.vector_length = vector_length;
    struct ToPartialProblemOutput output = to_partial_problem(input);

    PyObject *py_constraint_matrix = build_2D_python_list(output.constraint_matrix, output.constraint_matrix_rows,
                                                          output.constraint_matrix_columns);
    PyObject *py_right_side = build_2D_python_list(output.right_side, output.right_side_rows,
                                                   output.right_side_columns);

    return Py_BuildValue("(OOii)", py_constraint_matrix, py_right_side, output.nz, output.nx);
}

static PyObject *get_inverse_MPEC_model(PyObject *self, PyObject *args) {
    Py_ssize_t vector_x0_length, matrixA_rows, matrixA_columns, matrixB_rows, matrixB_columns, matrixC_rows,
            matrixC_columns, vector_b_length, vector_c_length;
    PyListObject *vector_x0, *matrixA, *matrixB, *matrixC, *vector_b, *vector_c;
    PyObject *temp_value;
    if (!PyArg_ParseTuple(args, "OOOOOO", &vector_x0, &matrixA, &matrixB, &matrixC, &vector_b, &vector_c)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return NULL;
    }

    if (!PyList_Check(vector_x0)) {
        PyErr_Format(PyExc_TypeError, "The vector_x0 must be of list or subtype of list");
        return NULL;
    }
    vector_x0_length = PyList_Size(vector_x0);
    if (vector_x0_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_vector_x0 = build_1D_c_array(vector_x0, vector_x0_length);

    if (!PyList_Check(matrixA)) {
        PyErr_Format(PyExc_TypeError, "The matrix must be of list or subtype of list");
        return NULL;
    }

    matrixA_rows = PyList_Size(matrixA);
    if (matrixA_rows <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(matrixA, 0);
    matrixA_columns = PyList_Size(temp_value);
    if (matrixA_columns <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_matrixA = build_2D_c_array(matrixA, matrixA_rows, matrixA_columns);

    matrixB_rows = PyList_Size(matrixB);
    if (matrixB_rows <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(matrixB, 0);
    matrixB_columns = PyList_Size(temp_value);
    if (matrixB_columns <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_matrixB = build_2D_c_array(matrixB, matrixB_rows, matrixB_columns);

    matrixC_rows = PyList_Size(matrixC);
    if (matrixC_rows <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    temp_value = PyList_GetItem(matrixC, 0);
    matrixC_columns = PyList_Size(temp_value);
    if (matrixC_columns <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float **c_matrixC = build_2D_c_array(matrixC, matrixC_rows, matrixC_columns);

    if (!PyList_Check(vector_b)) {
        PyErr_Format(PyExc_TypeError, "The vector_b must be of list or subtype of list");
        return NULL;
    }
    vector_b_length = PyList_Size(vector_b);
    if (vector_b_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_vector_b = build_1D_c_array(vector_b, vector_b_length);

    if (!PyList_Check(vector_c)) {
        PyErr_Format(PyExc_TypeError, "The vector_c must be of list or subtype of list");
        return NULL;
    }
    vector_c_length = PyList_Size(vector_c);
    if (vector_c_length <= 0) {
        PyErr_Format(PyExc_TypeError, "One of iterable arguments has zero length");
        return NULL;
    }
    float *c_vector_c = build_1D_c_array(vector_c, vector_c_length);

    struct MPEC_solver_input input;
    input.vector_x0 = c_vector_x0;
    input.vector_x0_length = vector_x0_length;
    input.matrixA = c_matrixA;
    input.matrixA_rows = matrixA_rows;
    input.matrixA_columns = matrixA_columns;
    input.matrixB = c_matrixB;
    input.matrixB_rows = matrixB_rows;
    input.matrixB_columns = matrixB_columns;
    input.matrixC = c_matrixC;
    input.matrixC_rows = matrixC_rows;
    input.matrixC_columns = matrixC_columns;
    input.vector_b = c_vector_b;
    input.vector_b_length = vector_b_length;
    input.vector_c = c_vector_c;
    input.vector_c_length = vector_c_length;
    struct MPEC_solver_output output = solve_inverse_via_MPEC(input);

    PyObject *py_objfun_coeffs = build_1D_python_list(output.objfun_coeffs,
                                                      output.objfun_coeffs_length);
    PyObject *py_constrs_matrix = build_2D_python_list(output.constrs_matrix, output.constrs_matrix_rows,
                                                       output.constrs_matrix_columns);
    PyObject *py_sense_array = build_1D_python_list(output.sense_array,
                                                    output.sense_array_length);
    PyObject *py_right_hand_side = build_1D_python_list(output.right_hand_side,
                                                        output.right_hand_side_length);
    PyObject *py_bounds = build_2D_python_list(output.bounds, output.bounds_rows,
                                               output.bounds_columns);

    return Py_BuildValue("(OOOOOii)", py_objfun_coeffs, py_constrs_matrix, py_sense_array,
                         py_right_hand_side, py_bounds, output.nx, output.nbin);
}

static PyMethodDef invp_methods[] = {
        {"to_standard_form",              to_standard_form,              METH_VARARGS, "default docs"},
        {"get_inverse_cost_vector_model", get_inverse_cost_vector_model, METH_VARARGS, "default docs"},
        {"to_canonical_form",             module_to_canonical_form,      METH_VARARGS, "default docs"},
        {"get_partial_problem_model",     get_partial_problem_model,     METH_VARARGS, "default docs"},
        {"get_inverse_MPEC_model",        get_inverse_MPEC_model,        METH_VARARGS, "default docs"},
        {NULL, NULL, 0, NULL},
};

static struct PyModuleDef invp = {
        PyModuleDef_HEAD_INIT,
        "invp",
        "Inverse Programming Module",
        -1,
        invp_methods
};

PyMODINIT_FUNC PyInit_invp(void) {
    return PyModule_Create(&invp);
}