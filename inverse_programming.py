from ctypes import *
import traceback
import sys

inf = 99999999.0
undefined = 67108864.0


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


class InverseProgrammingSolver:
    def __init__(self):

        self.inf = inf
        self.undefined = undefined

        # TODO убедиться, что работает во всех системах
        # TODO сделать отдельную функцию для компиляции самого этого файла на случай если его не окажется в файловой системе
        try:
            lib = CDLL('./compiled.so')
        except OSError as e:
            print("The compiled file cannot be used. Traceback and original exception:")
            traceback.print_exc()
            exit(1)

        self._lib = lib

        self._prepare_ctypes_functions()

    def _prepare_ctypes_functions(self):
        self._prepare_remaster_basic_input()
        self._prepare_function2()
        self._prepare_to_canonical_form()
        self._prepare_to_partial_problem()
        self._prepare_solve_inverse_via_MPEC()
        self._prepare_get_P_matrix()

    def _prepare_remaster_basic_input(self):
        self._lib.remaster_basic_input.restype = RemasteredInput
        self._lib.remaster_basic_input.argtypes = [BasicInput]

    def validate_remaster_basic_input_input(self, coefficients, constraint_matrix, right_side_matrix, bounds):
        validated_input = BasicInput()
        prepared_coefficients = [float(i) for i in coefficients]
        validated_input.coefficients = (c_float * len(coefficients))(*prepared_coefficients)
        validated_input.coefficients_length = len(coefficients)

        prepared_constraint_matrix = (POINTER(c_float) * len(constraint_matrix))()
        for i in range(len(constraint_matrix)):
            resource = [float(j) for j in constraint_matrix[i]]
            prepared_constraint_matrix[i] = (c_float * len(resource))(*resource)
        validated_input.constraint_matrix = prepared_constraint_matrix
        validated_input.constraint_matrix_length = len(constraint_matrix)

        prepared_right_side_matrix = (POINTER(c_float) * len(right_side_matrix))()
        for i in range(len(right_side_matrix)):
            resource = [float(j) for j in right_side_matrix[i]]
            prepared_right_side_matrix[i] = (c_float * len(resource))(*resource)
        validated_input.right_side_matrix = prepared_right_side_matrix
        validated_input.right_side_matrix_length = len(right_side_matrix)

        prepared_bounds = (POINTER(c_float) * len(bounds))()
        for i in range(len(bounds)):
            resource = [float(j) for j in bounds[i]]
            prepared_bounds[i] = (c_float * len(resource))(*resource)
        validated_input.bounds = prepared_bounds
        validated_input.bounds_length = len(bounds)

        return validated_input

    def validate_remaster_basic_input_output(self, output):
        coefficients = [float(output.coefficients[i]) for i in range(output.coefficients_length)]
        constraint_matrix = [[float(output.constraint_matrix[j][i]) for i in range(output.coefficients_length)]
                             for j in range(output.constraint_matrix_length)]
        right_side_array = [float(output.right_side_array[i]) for i in range(output.right_side_array_length)]
        return coefficients, constraint_matrix, right_side_array

    def remaster_basic_input(self, coefficients, constraint_matrix, right_side_array, bounds):
        validated_input = self.validate_remaster_basic_input_input(coefficients, constraint_matrix, right_side_array,
                                                                   bounds)

        output = self._lib.remaster_basic_input(validated_input)

        return self.validate_remaster_basic_input_output(output)

    def _prepare_function2(self):
        self._lib.function2.restype = AnotherOneOutput
        self._lib.function2.argtypes = [RemasteredInput, FloatVector]

    def validate_function2_input(self, coefficients, constraint_matrix, right_side_array, x0):
        validated_input = RemasteredInput()
        validated_vector = FloatVector()
        prepared_coefficients = [float(i) for i in coefficients]
        validated_input.coefficients = (c_float * len(coefficients))(*prepared_coefficients)
        validated_input.coefficients_length = len(coefficients)

        prepared_constraint_matrix = (POINTER(c_float) * len(constraint_matrix))()
        for i in range(len(constraint_matrix)):
            resource = [float(j) for j in constraint_matrix[i]]
            prepared_constraint_matrix[i] = (c_float * len(resource))(*resource)
        validated_input.constraint_matrix = prepared_constraint_matrix
        validated_input.constraint_matrix_length = len(constraint_matrix)

        prepared_right_side_array = [float(j) for j in right_side_array]
        validated_input.right_side_array = (c_float * len(right_side_array))(*prepared_right_side_array)
        validated_input.right_side_array_length = len(right_side_array)

        prepared_x0 = [float(i) for i in x0]
        validated_vector.values = (c_float * len(x0))(*prepared_x0)
        return validated_input, validated_vector

    def validate_function2_output(self, output):
        vars_length = output.vars_length
        constrs_length = output.constrs_length
        obj = [float(output.obj[i]) for i in range(output.obj_length)]
        newA = [[float(output.newA[j][i]) for i in range(output.newA_length)]
                for j in range(3 * vars_length)]
        sense = [int(output.sense[i]) for i in range(output.sense_length)]
        rhs = [float(output.rhs[i]) for i in range(output.rhs_length)]
        lb = [int(output.lb[i]) for i in range(output.lb_length)]

        return vars_length, constrs_length, obj, newA, sense, rhs, lb

    def function2(self, coefficients, constraint_matrix, right_side_array, x0):
        validated_input, validated_vector = self.validate_function2_input(coefficients, constraint_matrix,
                                                                          right_side_array, x0)

        output = self._lib.function2(validated_input, validated_vector)

        return self.validate_function2_output(output)

    def _prepare_to_canonical_form(self):
        self._lib.to_canonical_form.restype = CanonicalForm
        self._lib.to_canonical_form.argtypes = [RemasteredInput]

    def validate_to_canonical_form_input(self):
        pass

    def validate_to_canonical_form_output(self):
        pass

    def to_canonical_form(self):
        pass

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

    def _prepare_to_partial_problem(self):
        self._lib.to_partial_problem.restype = ToPartialProblemOutput
        self._lib.to_partial_problem.argtypes = [ToPartialProblemInput]

    def validate_to_partial_problem_input(self):
        pass

    def validate_to_partial_problem_output(self):
        pass

    def to_partial_problem(self):
        pass

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

    def _prepare_solve_inverse_via_MPEC(self):
        self._lib.solve_inverse_via_MPEC.restype = MPEC_solver_output
        self._lib.solve_inverse_via_MPEC.argtypes = [MPEC_solver_input]

    def validate_solve_inverse_via_MPEC_input(self):
        pass

    def validate_solve_inverse_via_MPEC_output(self):
        pass

    def solve_inverse_via_MPEC(self):
        pass

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

    def _prepare_get_P_matrix(self):
        self._lib.get_P_matrix.restype = P_matrix
        self._lib.get_P_matrix.argtypes = [c_int, c_int]

    def validate_get_P_matrix_input(self):
        pass

    def validate_get_P_matrix_output(self):
        pass

    def get_P_matrix(self):
        pass
    # m = lib.get_P_matrix(5, 3)
    # print(m.matrix_side)
    # for i in range(m.matrix_side):
    #     for j in range(m.matrix_side):
    #         print(m.matrix[i][j], end=' ')
    #     print()


sys.modules[__name__] = InverseProgrammingSolver()
