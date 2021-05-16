import inverse_programming as invp

if __name__ == '__main__':
    inf = invp.inf
    undefined = invp.undefined
    c = [1, 3, -1]
    A = [[1, 2, 0], [3, 3, 1], [-1, 2, 9]]
    RS = [[4, -1], [-5, 1], [7, 0]]
    bounds = [[-1, 3], [-inf, 10], [-inf, +inf]]
    coefficients, constraint_matrix, right_side_array = invp.remaster_basic_input(c, A, RS, bounds)
    print(coefficients)
    print(constraint_matrix)
    print(right_side_array)

    x0 = [-1, 10, -32]
    vars_length, constrs_length, obj, newA, sense, rhs, lb = invp.function2(coefficients, constraint_matrix,
                                                                            right_side_array, x0)
    print(vars_length)
    print(constrs_length)
    print(obj)
    print(newA)
    print(sense)
    print(rhs)
    print(lb)
