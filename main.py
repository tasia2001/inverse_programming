import inverse_programming as invp

if __name__ == '__main__':
    inf = invp.inf
    undefined = invp.undefined

    c = [-2, 3, -6, -1]
    A = [[2, 1, -2, 1], [-2, -1, 2, -1], [-1, -2, 0, -4], [1, -1, 2, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
         [0, 0, 0, 1]]
    RS = [24, -24, -22, 0, 0, 0, 0, 0]
    x0 = [26 / 3, 20 / 3, 0, 0]

    vars_length, constrs_length, obj, newA, sense, rhs, lb = invp.function2(c, A, RS, x0)
    print(vars_length)
    print(constrs_length)
    print(obj)
    print(newA)
    print(sense)
    print(rhs)
    print(lb)
