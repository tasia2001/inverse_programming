from __future__ import division

import sys
from ctypes import *
import traceback
import pyomo.environ as pyo
import gurobipy

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
        ('vector_length', c_int),
    ]


class ToPartialProblemOutput(Structure):
    _fields_ = [
        ('constraint_matrix', POINTER(POINTER(c_float))),
        ('constraint_matrix_rows', c_int),
        ('constraint_matrix_columns', c_int),
        ('right_side', POINTER(POINTER(c_float))),
        ('right_side_rows', c_int),
        ('right_side_columns', c_int),
        ('nz', c_int),
        ('nx', c_int),
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
        ('nx', c_int),
        ('nbin', c_int),
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
        self._prepare_pyomo_functions()

    def _prepare_pyomo_functions(self):
        self._prepare_pyomo_function1()
        self._prepare_pyomo_get_inverse_cost_vector_model()
        self._prepare_pyomo_function3()

    def _prepare_pyomo_function1(self):
        model1 = pyo.AbstractModel()

        model1.m = pyo.Param(within=pyo.NonNegativeIntegers)  # количество ограничений
        model1.n = pyo.Param(within=pyo.NonNegativeIntegers)  # количество переменных

        model1.I = pyo.RangeSet(1, model1.m)  # индексы ограничений
        model1.J = pyo.RangeSet(1, model1.n)  # индексы переменных

        model1.a = pyo.Param(model1.I, model1.J)  # объявляем матрицу ограничений
        model1.b = pyo.Param(model1.I)  # правые части
        model1.c = pyo.Param(model1.J)  # коэффициенты цф
        model1.lb = pyo.Param(model1.J)  # нижние границы
        model1.sense = pyo.Param(model1.I)  # знаки в ограничениях

        # the next line declares a variable indexed by the set J
        model1.x = pyo.Var(model1.J)

        def obj_expression(m):
            return pyo.summation(m.c, m.x)

        model1.OBJ = pyo.Objective(rule=obj_expression)

        def ax_constraint_rule(m, i):
            if m.sense[i] == -1:
                return sum(m.a[i, j] * m.x[j] for j in m.J) <= m.b[i]
            if m.sense[i] == 0:
                return sum(m.a[i, j] * m.x[j] for j in m.J) == m.b[i]
            if m.sense[i] == 1:
                return sum(m.a[i, j] * m.x[j] for j in m.J) >= m.b[i]

        # the next line creates one constraint for each member of the set model.I
        model1.AxbConstraint = pyo.Constraint(model1.I, rule=ax_constraint_rule)

        def bounds_rule(m, j):
            return (m.lb[j], m.x[j], None)

        model1.boundx = pyo.Constraint(model1.J, rule=bounds_rule)
        self._inverse_cost_vector_model = model1

    def solve_tool_inverse_cost_vector(self, m, n, v, con, obj, newA, sense, rhs, lb, solver):
        data = {
            None: {
                'm': {None: m},
                'n': {None: n},
                'a': {(i + 1, j + 1): newA[i][j] for i in range(m) for j in range(n)},
                'b': {i + 1: rhs[i] for i in range(m)},
                'c': {i + 1: obj[i] for i in range(n)},
                'lb': {i + 1: lb[i] for i in range(n)},
                'sense': {i + 1: sense[i] for i in range(m)}
            }
        }
        ins = self._inverse_cost_vector_model.create_instance(data)
        opt = pyo.SolverFactory(solver, solver_io='python')
        opt.solve(ins)
        res = [pyo.value(ins.x[i]) for i in range(v + 1, 2 * v + 1)]
        return res

    def _prepare_pyomo_get_inverse_cost_vector_model(self):
        model2 = pyo.AbstractModel()

        model2.m = pyo.Param(within=pyo.NonNegativeIntegers)  # количество ограничений
        model2.n = pyo.Param(within=pyo.NonNegativeIntegers)  # количество переменных
        model2.nx = pyo.Param(within=pyo.NonNegativeIntegers)  # количество обычных переменных
        model2.nbin = pyo.Param(within=pyo.NonNegativeIntegers)  # количество бинарных переменных

        model2.I = pyo.RangeSet(1, model2.m)  # индексы ограничений
        model2.J = pyo.RangeSet(1, model2.n)  # индексы переменных
        model2.Jx = pyo.RangeSet(1, model2.nx)  # индексы обычных переменных
        model2.Jbin = pyo.RangeSet(1, model2.nbin)  # индексы бинарных переменных

        model2.Ax = pyo.Param(model2.I, model2.Jx)  # объявляем матрицу ограничений для x
        model2.Abin = pyo.Param(model2.I, model2.Jbin)  # объявляем матрицу ограничений для bin
        model2.b = pyo.Param(model2.I)  # правые части
        model2.c = pyo.Param(model2.Jx)  # коэффициенты цф
        model2.bounds = pyo.Param(model2.J)  # границы
        model2.sense = pyo.Param(model2.I)  # знаки в ограничениях

        # the next line declares a variable indexed by the set J
        model2.x = pyo.Var(model2.Jx, domain=pyo.Reals)
        model2.bin = pyo.Var(model2.Jbin, domain=pyo.Binary)

        def obj_expression(m):
            return pyo.summation(m.c, m.x)

        model2.OBJ = pyo.Objective(rule=obj_expression)

        def ax_constraint_rule(m, i):
            if m.sense[i] == -1:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Abin[i, j] * m.bin[j] for j in m.Jbin) <= m.b[i]
            if m.sense[i] == 0:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Abin[i, j] * m.bin[j] for j in m.Jbin) == m.b[i]
            if m.sense[i] == 1:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Abin[i, j] * m.bin[j] for j in m.Jbin) >= m.b[i]

        # the next line creates one constraint for each member of the set model.I
        model2.AxbConstraint = pyo.Constraint(model2.I, rule=ax_constraint_rule)

        def bounds_rule(m, j):
            return (m.bounds[j][0], m.x[j], m.bounds[j][1])

        model2.boundx = pyo.Constraint(model2.Jx, rule=bounds_rule)
        self._solve_tool_inverse_MPEC_model = model2

    def solve_tool_inverse_MPEC(self, m, n, nx, nbin, obj, A, sense, rhs, bounds, solver):
        data = {None: {
            'm': {None: m},
            'n': {None: n},
            'nx': {None: nx},
            'nbin': {None: nbin},
            'Ax': {(i + 1, j + 1): A[i][j] for i in range(m) for j in range(nx)},
            'Abin': {(i + 1, j + 1): A[i][j + nx] for i in range(m) for j in range(nbin)},
            'b': {i + 1: rhs[i] for i in range(m)},
            'c': {i + 1: obj[i] for i in range(nx)},
            'bounds': {i + 1: bounds[i] for i in range(nx)},
            'sense': {i + 1: sense[i] for i in range(m)}
        }}
        ins = self._solve_tool_inverse_MPEC_model.create_instance(data)
        opt = pyo.SolverFactory(solver, solver_io="python")
        opt.solve(ins)
        res = [pyo.value(ins.x[i]) for i in range(1, nx + 1)]
        return res

    def _prepare_pyomo_function3(self):
        model3 = pyo.AbstractModel()

        model3.m = pyo.Param(within=pyo.NonNegativeIntegers)  # количество ограничений
        model3.n = pyo.Param(within=pyo.NonNegativeIntegers)  # количество переменных всего
        model3.nx = pyo.Param(within=pyo.NonNegativeIntegers)  # количество обычных переменных
        model3.nz = pyo.Param(within=pyo.NonNegativeIntegers)  # количество бинарных переменных

        model3.I = pyo.RangeSet(1, model3.m)  # индексы ограничений
        model3.J = pyo.RangeSet(1, model3.n)  # индексы переменных
        model3.Jx = pyo.RangeSet(1, model3.nx)  # индексы обычных переменных
        model3.Jz = pyo.RangeSet(1, model3.nz)  # индексы бинарных переменных
        model3.two = pyo.RangeSet(1, 2)

        model3.Ax = pyo.Param(model3.I, model3.Jx)  # объявляем матрицу ограничений для x
        model3.Az = pyo.Param(model3.I, model3.Jz)  # объявляем матрицу ограничений для z
        model3.rhs = pyo.Param(model3.I, model3.two)  # rhs
        model3.P = pyo.Param(model3.Jx, model3.Jx)  # коэффициенты цф для x

        # variables
        model3.x = pyo.Var(model3.Jx, domain=pyo.Reals)
        model3.z = pyo.Var(model3.Jz, domain=pyo.Binary)

        def obj_expression(m):
            return sum(m.x[i] * m.x[i + m.nz] for i in m.Jz)

        model3.OBJ = pyo.Objective(rule=obj_expression)

        def ax_constraint_rule(m, i):
            if m.rhs[i, 2] == -1:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Az[i, j] * m.z[j] for j in m.Jz) <= m.rhs[i, 1]
            if m.rhs[i, 2] == 0:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Az[i, j] * m.z[j] for j in m.Jz) == m.rhs[i, 1]
            if m.rhs[i, 2] == 1:
                return sum(m.Ax[i, j] * m.x[j] for j in m.Jx) + sum(m.Az[i, j] * m.z[j] for j in m.Jz) >= m.rhs[i, 1]

        # the next line creates one constraint for each member of the set model.I
        model3.AxbConstraint = pyo.Constraint(model3.I, rule=ax_constraint_rule)
        self._solve_tool_partial_inverse_model = model3

    def solve_tool_partial_inverse(self, m, n, nx, nz, A, rhs, P):
        data = {None: {
            'm': {None: m},
            'n': {None: n},
            'nx': {None: nx},
            'nz': {None: nz},
            'Ax': {(i + 1, j + 1): A[i][j] for i in range(m) for j in range(nx)},
            'Az': {(i + 1, j + 1): A[i][j + nx] for i in range(m) for j in range(nz)},
            'rhs': {(i + 1, j + 1): rhs[i][j] for i in range(m) for j in range(2)},
            'P': {(i + 1, j + 1): P[i][j] for i in range(nx) for j in range(nx)}
        }}
        ins = self._solve_tool_partial_inverse_model.create_instance(data)
        optimizer = pyo.SolverFactory('gurobi', solver_io='python')
        optimizer.options["NonConvex"] = 2
        optimizer.solve(ins)
        res = [pyo.value(ins.x[i]) for i in range(1, 2 * nz + 1)]
        return res

    def _prepare_ctypes_functions(self):
        self._prepare_to_standard_form()
        self._prepare_get_inverse_cost_vector_model()
        self._prepare_to_canonical_form()
        self._prepare_get_partial_problem_model()
        self._prepare_get_inverse_MPEC_model()
        self._prepare_get_QPmatrix_for_partial()

    def _prepare_to_standard_form(self):
        self._lib.remaster_basic_input.restype = RemasteredInput
        self._lib.remaster_basic_input.argtypes = [BasicInput]

    def validate_to_standard_form_input(self, coefficients, constraint_matrix, right_side_matrix, bounds):
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

    def validate_to_standard_form_output(self, output):
        coefficients = [float(output.coefficients[i]) for i in range(output.coefficients_length)]
        constraint_matrix = [[float(output.constraint_matrix[j][i]) for i in range(output.coefficients_length)]
                             for j in range(output.constraint_matrix_length)]
        right_side_array = [float(output.right_side_array[i]) for i in range(output.right_side_array_length)]
        return coefficients, constraint_matrix, right_side_array

    def to_standard_form(self, coefficients, constraint_matrix, right_side_array, bounds):
        validated_input = self.validate_to_standard_form_input(coefficients, constraint_matrix, right_side_array,
                                                               bounds)

        output = self._lib.remaster_basic_input(validated_input)

        return self.validate_to_standard_form_output(output)

    def _prepare_get_inverse_cost_vector_model(self):
        self._lib.function2.restype = AnotherOneOutput
        self._lib.function2.argtypes = [RemasteredInput, FloatVector]

    def validate_get_inverse_cost_vector_model_input(self, coefficients, constraint_matrix, right_side_array, x0):
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

    def validate_get_inverse_cost_vector_model_output(self, output):
        vars_length = output.vars_length
        constrs_length = output.constrs_length
        obj = [float(output.obj[i]) for i in range(output.obj_length)]
        newA = [[float(output.newA[j][i]) for i in range(output.newA_length)]
                for j in range(3 * vars_length)]
        sense = [int(output.sense[i]) for i in range(output.sense_length)]
        rhs = [float(output.rhs[i]) for i in range(output.rhs_length)]
        lb = [int(output.lb[i]) for i in range(output.lb_length)]
        m = len(newA)
        n = len(obj)
        return m, n, vars_length, constrs_length, obj, newA, sense, rhs, lb

    def get_inverse_cost_vector_model(self, coefficients, constraint_matrix, right_side_array, x0):
        validated_input, validated_vector = self.validate_get_inverse_cost_vector_model_input(coefficients,
                                                                                              constraint_matrix,
                                                                                              right_side_array, x0)

        output = self._lib.function2(validated_input, validated_vector)

        return self.validate_get_inverse_cost_vector_model_output(output)

    def solve_inverse_cost_vector_problem(self, *args, solver='gurobi'):
        values = self.get_inverse_cost_vector_model(*args)
        return self.solve_tool_inverse_cost_vector(*values, solver)

    def _prepare_to_canonical_form(self):
        self._lib.to_canonical_form.restype = CanonicalForm
        self._lib.to_canonical_form.argtypes = [RemasteredInput]

    def validate_to_canonical_form_input(self):
        pass

    def validate_to_canonical_form_output(self):
        pass

    def to_canonical_form(self):
        pass
    # def solve_inverse_cost_vector_problem(self, *args, solver='gurobi'):
    #     values = self.get_inverse_MPEC_model(*args)
    #     return self.solve_tool_inverse_cost_vector(*values, solver)

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

    def _prepare_get_partial_problem_model(self):
        self._lib.to_partial_problem.restype = ToPartialProblemOutput
        self._lib.to_partial_problem.argtypes = [ToPartialProblemInput]

    def validate_get_partial_problem_model_input(self):
        pass

    def validate_get_partial_problem_model_output(self):
        pass

    def get_partial_problem_model(self):
        pass
    # def solve_inverse_cost_vector_problem(self, *args, solver='gurobi'):
    #     values = self.get_inverse_MPEC_model(*args)
    #     return self.solve_tool_inverse_cost_vector(*values, solver)

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

    def _prepare_get_inverse_MPEC_model(self):
        self._lib.solve_inverse_via_MPEC.restype = MPEC_solver_output
        self._lib.solve_inverse_via_MPEC.argtypes = [MPEC_solver_input]

    def validate_get_inverse_MPEC_model_input(self):
        pass

    def validate_get_inverse_MPEC_model_output(self):
        pass

    def get_inverse_MPEC_model(self):
        pass

    # def solve_inverse_cost_vector_problem(self, *args, solver='gurobi'):
    #     values = self.get_inverse_MPEC_model(*args)
    #     return self.solve_tool_inverse_cost_vector(*values, solver)

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

    def _prepare_get_QPmatrix_for_partial(self):
        self._lib.get_P_matrix.restype = P_matrix
        self._lib.get_P_matrix.argtypes = [c_int, c_int]

    def validate_get_QPmatrix_for_partial_input(self):
        pass

    def validate_get_QPmatrix_for_partial_output(self):
        pass

    def get_QPmatrix_for_partial(self):
        m = self._lib.get_P_matrix(5, 3)
        print(m.matrix_side)
        for i in range(m.matrix_side):
            for j in range(m.matrix_side):
                print(m.matrix[i][j], end=' ')
            print()
    # def solve_inverse_cost_vector_problem(self, *args, solver='gurobi'):
    #     values = self.get_inverse_MPEC_model(*args)
    #     return self.solve_tool_inverse_cost_vector(*values, solver)


    # m = lib.get_P_matrix(5, 3)
    # print(m.matrix_side)
    # for i in range(m.matrix_side):
    #     for j in range(m.matrix_side):
    #         print(m.matrix[i][j], end=' ')
    #     print()


sys.modules[__name__] = InverseProgrammingSolver()
