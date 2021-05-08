from ctypes import *

inf = 99999999.0

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


class FloatVector(Structure):
    _fields_ = [
        ('values', POINTER(c_float)),
    ]


if __name__ == '__main__':
    lib.remaster_basic_input.restype = RemasteredInput
    lib.remaster_basic_input.argtypes = [BasicInput]

    lib.function2.restype = AnotherOneOutput
    lib.function2.argtypes = [RemasteredInput, FloatVector]

    basic_input = BasicInput()
    # coefficients = [1.0, 3.0, -1.0]
    coefficients = eval(input())
    basic_input.coefficients = (c_float * 3)(*coefficients)
    basic_input.coefficients_length = 3

    # resource = [[1.0, 2.0, 0.0], [3.0, 3.0, 1.0], [-1.0, 2.0, 9.0]]
    resource = eval(input())
    constraint_matrix = (POINTER(c_float) * 3)()
    for i in range(3):
        constraint_matrix[i] = (c_float * 3)(*resource[i])
    basic_input.constraint_matrix = constraint_matrix
    basic_input.constraint_matrix_length = 3

    # resource = [[4.0, -1.0], [-5.0, 1.0], [7.0, 0.0]]
    resource = eval(input())
    right_side_matrix = (POINTER(c_float) * 3)()
    for i in range(3):
        right_side_matrix[i] = (c_float * 2)(*resource[i])
    basic_input.right_side_matrix = right_side_matrix
    basic_input.right_side_matrix_length = 3

    # resource = [[-1.0, 3.0], [-inf, 10.0], [-inf, inf]]
    resource = eval(input())
    bounds = (POINTER(c_float) * 3)()
    for i in range(3):
        bounds[i] = (c_float * 2)(*resource[i])
    basic_input.bounds = bounds
    basic_input.bounds_length = 3

    result = lib.remaster_basic_input(basic_input)

    float_vector = FloatVector()
    # resource = [-1.0, 10.0, -32.0]
    resource = eval(input())
    float_vector.values = (c_float * 3)(*resource)

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

    new_result = lib.function2(result, float_vector)
    print('vars_length: ', new_result.vars_length)
    print('constrs_length: ', new_result.constrs_length)
    print('obj_length: ', new_result.obj_length)
    print('obj: ')
    for i in range(new_result.obj_length):
        print(new_result.obj[i], end=' ')
    print()
    print('newA_length: ', new_result.newA_length)
    print('newA: ')
    for i in range(3 * new_result.vars_length):
        for j in range(new_result.newA_length):
            print(new_result.newA[i][j], end=' ')
        print()
    print('sense_length: ', new_result.sense_length)
    print('sense: ')
    for i in range(new_result.sense_length):
        print(new_result.sense[i], end=' ')
    print()
    print('rhs_length: ', new_result.rhs_length)
    for i in range(new_result.rhs_length):
        print(new_result.rhs[i], end=' ')
    print()
    print('lb_length: ', new_result.lb_length)
    for i in range(new_result.lb_length):
        print(new_result.lb[i], end=' ')
    print()
