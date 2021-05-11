from ctypes import *

inf = 99999999.0
undefined = 67108864.0

lib = CDLL('./compiled.so')


class BasicInput(Structure):
    _fields_ = [
        ('coefficients', POINTER(c_float)),
        ('coefficients_length', c_int),
        ('constraint_matrix', POINTER(POINTER(c_float))),
        ('constraint_matrix_length', c_int),
        ('right_side_matrix', POINTER(POINTER(c_float))),
        ('right_side_matrix_length', c_int),
        ('bounds', POINTER(POINTER(c_float))),
        ('bounds_length', c_int),
    ]


class RemasteredInput(Structure):
    _fields_ = [
        ('coefficients', POINTER(c_float)),
        ('coefficients_length', c_int),
        ('constraint_matrix', POINTER(POINTER(c_float))),
        ('constraint_matrix_length', c_int),
        ('right_side_array', POINTER(c_float)),
        ('right_side_array_length', c_int),
    ]


class AnotherOneOutput(Structure):
    _fields_ = [
        ('vars_length', c_int),
        ('constrs_length', c_int),
        ('obj', POINTER(c_float)),
        ('obj_length', c_int),
        ('newA', POINTER(POINTER(c_float))),
        ('newA_length', c_int),
        ('sense', POINTER(c_float)),
        ('sense_length', c_int),
        ('rhs', POINTER(c_float)),
        ('rhs_length', c_int),
        ('lb', POINTER(c_float)),
        ('lb_length', c_int),
    ]


class CanonicalForm(Structure):
    _fields_ = [
        ('obj_fun_coefficients', POINTER(c_float)),
        ('obj_fun_coefficients_length', c_int),
        ('constrs_matrix', POINTER(POINTER(c_float))),
        ('constrs_matrix_rows', c_int),
        ('constrs_matrix_columns', c_int),
        ('right_hand_side', POINTER(c_float)),
        ('right_hand_side_length', c_int),
        ('low_bounds', POINTER(c_float)),
        ('low_bounds_length', c_int),
    ]


class ToPartialProblemInput(Structure):
    _fields_ = [
        ('matrix', POINTER(POINTER(c_float))),
        ('matrix_rows', c_int),
        ('matrix_columns', c_int),
        ('right_side_array', POINTER(c_float)),
        ('right_side_array_length', c_int),
        ('vector', POINTER(c_float)),
    ]


class ToPartialProblemOutput(Structure):
    _fields_ = [
        ('constraint_matrix', POINTER(POINTER(c_float))),
        ('constraint_matrix_rows', c_int),
        ('constraint_matrix_columns', c_int),
        ('right_side', POINTER(POINTER(c_float))),
        ('right_side_rows', c_int),
        ('right_side_columns', c_int),
    ]


class FloatVector(Structure):
    _fields_ = [
        ('values', POINTER(c_float)),
    ]


class MPEC_solver_input(Structure):
    _fields_ = [
        ('vector_x0', POINTER(c_float)),
        ('vector_x0_length', c_int),
        ('matrixA', POINTER(POINTER(c_float))),
        ('matrixA_rows', c_int),
        ('matrixA_columns', c_int),
        ('matrixB', POINTER(POINTER(c_float))),
        ('matrixB_rows', c_int),
        ('matrixB_columns', c_int),
        ('matrixC', POINTER(POINTER(c_float))),
        ('matrixC_rows', c_int),
        ('matrixC_columns', c_int),
        ('vector_b', POINTER(c_float)),
        ('vector_b_length', c_int),
        ('vector_c', POINTER(c_float)),
        ('vector_c_length', c_int),
    ]


class MPEC_solver_output(Structure):
    _fields_ = [
        ('objfun_coeffs', POINTER(c_float)),
        ('objfun_coeffs_length', c_int),
        ('constrs_matrix', POINTER(POINTER(c_float))),
        ('constrs_matrix_rows', c_int),
        ('constrs_matrix_columns', c_int),
        ('sense_array', POINTER(c_float)),
        ('sense_array_length', c_int),
        ('right_hand_side', POINTER(c_float)),
        ('right_hand_side_length', c_int),
        ('bounds', POINTER(POINTER(c_float))),
        ('bounds_rows', c_int),
        ('bounds_columns', c_int),
    ]


class P_matrix(Structure):
    _fields_ = [
        ('matrix', POINTER(POINTER(c_float))),
        ('matrix_side', c_int),
    ]


if __name__ == '__main__':
    lib.remaster_basic_input.restype = RemasteredInput
    lib.remaster_basic_input.argtypes = [BasicInput]

    lib.function2.restype = AnotherOneOutput
    lib.function2.argtypes = [RemasteredInput, FloatVector]

    lib.to_canonical_form.restype = CanonicalForm
    lib.to_canonical_form.argtypes = [RemasteredInput]

    lib.to_partial_problem.restype = ToPartialProblemOutput
    lib.to_partial_problem.argtypes = [ToPartialProblemInput]

    lib.solve_inverse_via_MPEC.restype = MPEC_solver_output
    lib.solve_inverse_via_MPEC.argtypes = [MPEC_solver_input]

    lib.get_P_matrix.restype = P_matrix
    lib.get_P_matrix.argtypes = [c_int, c_int]

    # coefficients = [1.0, 3.0, -1.0]
    # basic_input = BasicInput()
    # coefficients = eval(input())
    # basic_input.coefficients = (c_float * 3)(*coefficients)
    # basic_input.coefficients_length = 3

    # resource = [[1.0, 2.0, 0.0], [3.0, 3.0, 1.0], [-1.0, 2.0, 9.0]]
    # resource = eval(input())
    # constraint_matrix = (POINTER(c_float) * 3)()
    # for i in range(3):
    #     constraint_matrix[i] = (c_float * 3)(*resource[i])
    # basic_input.constraint_matrix = constraint_matrix
    # basic_input.constraint_matrix_length = 3

    # resource = [[4.0, -1.0], [-5.0, 1.0], [7.0, 0.0]]
    # resource = eval(input())
    # right_side_matrix = (POINTER(c_float) * 3)()
    # for i in range(3):
    #     right_side_matrix[i] = (c_float * 2)(*resource[i])
    # basic_input.right_side_matrix = right_side_matrix
    # basic_input.right_side_matrix_length = 3

    # resource = [[-1.0, 3.0], [-inf, 10.0], [-inf, inf]]
    # resource = eval(input())
    # bounds = (POINTER(c_float) * 3)()
    # for i in range(3):
    #     bounds[i] = (c_float * 2)(*resource[i])
    # basic_input.bounds = bounds
    # basic_input.bounds_length = 3

    # result = lib.remaster_basic_input(basic_input)

    # resource = [-1.0, 10.0, -32.0]
    # float_vector = FloatVector()
    # resource = eval(input())
    # float_vector.values = (c_float * 3)(*resource)

    # print("Исходный вектор: ", end='')
    # for i in range(result.coefficients_length):
    #     print(result.coefficients[i], end=' ')
    # print("\nНовая матрица ограничений: ")
    # for i in range(result.constraint_matrix_length):
    #     for j in range(result.coefficients_length):
    #         print(result.constraint_matrix[i][j], end=' ')
    #     print()
    # print("Новый массив правых частей: ")
    # for i in range(result.right_side_array_length):
    #     print(result.right_side_array[i], end=' ')
    # print()

    # new_result = lib.function2(result, float_vector)
    # print('vars_length: ', new_result.vars_length)
    # print('constrs_length: ', new_result.constrs_length)
    # print('obj_length: ', new_result.obj_length)
    # print('obj: ')
    # for i in range(new_result.obj_length):
    #     print(new_result.obj[i], end=' ')
    # print()
    # print('newA_length: ', new_result.newA_length)
    # print('newA: ')
    # for i in range(3 * new_result.vars_length):
    #     for j in range(new_result.newA_length):
    #         print(new_result.newA[i][j], end=' ')
    #     print()
    # print('sense_length: ', new_result.sense_length)
    # print('sense: ')
    # for i in range(new_result.sense_length):
    #     print(new_result.sense[i], end=' ')
    # print()
    # print('rhs_length: ', new_result.rhs_length)
    # for i in range(new_result.rhs_length):
    #     print(new_result.rhs[i], end=' ')
    # print()
    # print('lb_length: ', new_result.lb_length)
    # for i in range(new_result.lb_length):
    #     print(new_result.lb[i], end=' ')
    # print()

    # new_result = lib.to_canonical_form(result)
    # print('obj_fun_coefficients: ')
    # for i in range(new_result.obj_fun_coefficients_length):
    #     print(new_result.obj_fun_coefficients[i], end=' ')
    # print('\nconstrs_matrix: ')
    # for i in range(new_result.constrs_matrix_rows):
    #     for j in range(new_result.constrs_matrix_columns):
    #         print(new_result.constrs_matrix[i][j], end=' ')
    #     print()
    # print('right_hand_side: ')
    # for i in range(new_result.right_hand_side_length):
    #     print(new_result.right_hand_side[i], end=' ')
    # print('\nlow_bounds: ')
    # for i in range(new_result.low_bounds_length):
    #     print(new_result.low_bounds[i], end=' ')
    # print()

    # to_partial_problem_input = ToPartialProblemInput()
    #
    # # resource = [[1, 1], [-1, 1], [1, 0]]
    # resource = eval(input())
    # matrix = (POINTER(c_float) * len(resource))()
    # for i in range(len(resource)):
    #     matrix[i] = (c_float * len(resource[0]))(*resource[i])
    # to_partial_problem_input.matrix = matrix
    # to_partial_problem_input.matrix_rows = len(resource)
    # to_partial_problem_input.matrix_columns = len(resource[0]) if resource != [] else 0
    #
    # # resource = [undefined, 1]
    # resource = eval(input())
    # matrix = (c_float * len(resource))(*resource)
    # to_partial_problem_input.right_side_array = matrix
    # to_partial_problem_input.right_side_array_length = len(resource)
    #
    # # resource = [3, 1, 2]
    # resource = eval(input())
    # to_partial_problem_input.vector = (c_float * len(resource))(*resource)
    #
    # result = lib.to_partial_problem(to_partial_problem_input)
    # print('constraint_matrix: ')
    # for i in range(result.constraint_matrix_rows):
    #     for j in range(result.constraint_matrix_columns):
    #         print(result.constraint_matrix[i][j], end=' ')
    #     print()
    # print('\nright_side: ')
    # for i in range(result.right_side_rows):
    #     for j in range(result.right_side_columns):
    #         print(result.right_side[i][j], end=' ')
    #     print()

    # solver_input = MPEC_solver_input()
    #
    # # resource = [-4.0, 10.0, -6.0]
    # resource = eval(input())
    # vector_x0 = (c_float * len(resource))(*resource)
    # solver_input.vector_x0 = vector_x0
    # solver_input.vector_x0_length = len(resource)
    #
    # # resource = [[-6.0, -6.0, -9.0], [2.0, -4.0, 2.0], [7.0, -6.0, -4.0], [-6.0, 6.0, -7.0]]
    # resource = eval(input())
    # matrixA = (POINTER(c_float) * len(resource))()
    # for i in range(len(resource)):
    #     matrixA[i] = (c_float * len(resource[0]))(*resource[i])
    # solver_input.matrixA = matrixA
    # solver_input.matrixA_rows = len(resource)
    # solver_input.matrixA_columns = len(resource[0]) if resource != [] else 0
    #
    # # resource = [[7.0, -2.0, -5.0, -2.0], [1.0, -5.0, 5.0, 3.0], [0.0, -5.0, -3.0, 5.0]]
    # resource = eval(input())
    # matrixB = (POINTER(c_float) * len(resource))()
    # for i in range(len(resource)):
    #     matrixB[i] = (c_float * len(resource[0]))(*resource[i])
    # solver_input.matrixB = matrixB
    # solver_input.matrixB_rows = len(resource)
    # solver_input.matrixB_columns = len(resource[0]) if resource != [] else 0
    #
    # # resource = [[3.0, -8.0, 5.0], [-3.0, 9.0, -1.0]]
    # resource = eval(input())
    # matrixC = (POINTER(c_float) * len(resource))()
    # for i in range(len(resource)):
    #     matrixC[i] = (c_float * len(resource[0]))(*resource[i])
    # solver_input.matrixC = matrixC
    # solver_input.matrixC_rows = len(resource)
    # solver_input.matrixC_columns = len(resource[0]) if resource != [] else 0
    #
    # # resource = [10.0, 2.0, 7.0]
    # resource = eval(input())
    # vector_b = (c_float * len(resource))(*resource)
    # solver_input.vector_b = vector_b
    # solver_input.vector_b_length = len(resource)
    #
    # # resource = [-9.0, 6.0]
    # resource = eval(input())
    # vector_c = (c_float * len(resource))(*resource)
    # solver_input.vector_c = vector_c
    # solver_input.vector_c_length = len(resource)
    #
    # result = lib.solve_inverse_via_MPEC(solver_input)
    #
    # print('objfun_coeffs: ')
    # for i in range(result.objfun_coeffs_length):
    #     print(result.objfun_coeffs[i], end=' ')
    # print('\nconstrs_matrix: ')
    # for i in range(result.constrs_matrix_rows):
    #     for j in range(result.constrs_matrix_columns):
    #         print(result.constrs_matrix[i][j], end=' ')
    #     print()
    # print('sense_array: ')
    # for i in range(result.sense_array_length):
    #     print(result.sense_array[i], end=' ')
    # print('\nright_hand_side: ')
    # for i in range(result.right_hand_side_length):
    #     print(result.right_hand_side[i], end=' ')
    # print('\nbounds: ')
    # for i in range(result.bounds_rows):
    #     for j in range(result.bounds_columns):
    #         print(result.bounds[i][j], end=' ')
    #     print()

    m = lib.get_P_matrix(5, 3)
    print(m.matrix_side)
    for i in range(m.matrix_side):
        for j in range(m.matrix_side):
            print(m.matrix[i][j], end=' ')
        print()
