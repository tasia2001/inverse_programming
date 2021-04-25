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


class BasicOutput(Structure):
    _fields_ = [
        ('coefficients', POINTER(c_float)),
        ('coefficients_length', c_int),
        ('constraint_matrix', POINTER(POINTER(c_float))),
        ('constraint_matrix_length', c_int),
        ('right_side_array', POINTER(c_float)),
        ('right_side_array_length', c_int),
    ]


if __name__ == '__main__':
    lib.remake_basic_input.restype = BasicOutput
    lib.remake_basic_input.argtypes = [BasicInput]

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

    result = lib.remake_basic_input(basic_input)
    print("Исходный вектор: ", end='')
    for i in range(result.coefficients_length):
        print(result.coefficients[i], end=' ')
    print("\nНовая матрица ограничений: ")
    for i in range(result.constraint_matrix_length):
        for j in range(result.coefficients_length):
            print(result.constraint_matrix[i][j], end=' ')
        print()
    print("Новый массив правых частей: ")
    for i in range(result.right_side_array_length):
        print(result.right_side_array[i], end=' ')
    print()
