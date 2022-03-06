import numpy

epsilon = 10.0 ** -4


def initial():
    numpy.set_printoptions(suppress=True, precision=4, floatmode="fixed")
    b = numpy.array([[1.2], [2.2], [4.0], [0.0], [-1.2]])
    C = numpy.array([
        [0.01, 0.0, -0.02, 0.0, 0.0],
        [0.01, 0.01, -0.02, 0.0, 0.0],
        [0.0, 0.01, 0.01, 0.0, -0.02],
        [0.0, 0.0, 0.01, 0.01, 0.0],
        [0.0, 0.0, 0.0, 0.01, 0.01]
    ])
    D = numpy.array([
        [1.33, 0.21, 0.17, 0.12, -0.13],
        [-0.13, -1.33, 0.11, 0.17, 0.12],
        [0.12, -0.13, -1.33, 0.11, 0.17],
        [0.17, 0.12, -0.13, -1.33, 0.11],
        [0.11, 0.67, 0.12, -0.13, -1.33]
    ])
    A = 7 * C + D

    # # test example 1
    # A = numpy.array([
    #     [2.0, 1.0],
    #     [1.0, -2.0]
    # ])
    # b = numpy.array([[3.0], [1.0]])
    # # test example 1

    # # test example 2 with error
    # A = numpy.array([
    #     [-1.0, 0.5, 0.6],
    #     [0.0, -1.0, 0.5],
    #     [0.5, 0.0, -1.0]
    # ])
    # b = numpy.array([[1.0], [1.0], [1.0]])
    # # test example 2 with error

    return A, b


def output_initial(a_to_work, b_to_work):
    print(f'Входные данные:\n\nA:\n {a_to_work}')
    print(f'b:\n {b_to_work}\n\n')


def element_is_zero():
    print("Element == 0!")
    quit()


def get_matrix_b(a_to_work):
    n = a_to_work.shape[0]
    alpha = numpy.zeros((n, n))
    for i in range(n):
        for j in range(n):
            alpha[i, j] = - a_to_work[i, j] / a_to_work[i, i]
        alpha[i, i] = 0
    return alpha


def get_all_norms(a_to_work):
    n = a_to_work.shape[0]
    first_norm = max(numpy.absolute(a_to_work[i]).sum() for i in range(n))
    second_norm = max(numpy.absolute(a_to_work.T[j]).sum() for j in range(n))
    third_norm = ((a_to_work ** 2).sum()) ** (1 / 2)
    return first_norm, second_norm, third_norm


def method_of_simple_iterations(a_to_work, b_to_work):
    n = a_to_work.shape[0]

    for i in range(n):
        if a_to_work[i, i] == 0.0:
            element_is_zero()

    matrix_b = get_matrix_b(a_to_work)
    if not (min(get_all_norms(matrix_b)) < 1):
        print("Одно из достаточных условий сходимости не выполняется!")
        # quit(0)

    c_column = numpy.zeros((n, 1))
    for i in range(n):
        c_column[i] = b_to_work[i] / a_to_work[i, i]

    print(f'B:\n {matrix_b}')

    current_iteration_x = numpy.zeros((n, 1))

    count_of_iterations = 0
    delta_one = epsilon
    delta_second = epsilon

    while delta_one + delta_second > epsilon:
        previous_iteration_x = current_iteration_x.copy()
        current_iteration_x = c_column + matrix_b.dot(current_iteration_x)  # .dot - скалярное произведение

        delta_one = numpy.absolute((current_iteration_x - previous_iteration_x)).max()
        delta_second = numpy.absolute((a_to_work.dot(current_iteration_x) - b_to_work)).max()

        count_of_iterations += 1

    print(f"\nКоличество итераций, при расчете методом простых итераций: {count_of_iterations}")
    print(f'x: {current_iteration_x.T}\n\n')


def check_is_sufficient_condition_is_satisfied(a_to_work):  # ???
    n = a_to_work.shape[0]

    sum_of_row = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                sum_of_row += a_to_work[i][j]
        if abs(a_to_work[i][i]) < sum_of_row:
            print('Достаточное условие сходимости метода Зейделя не выполняется!')
            quit(0)
        else:
            sum_of_row = 0


def method_of_seidel(a_to_work, b_to_work):
    n = a_to_work.shape[0]

    # check_is_sufficient_condition_is_satisfied(a_to_work.copy())  # ???

    for i in range(n):
        if a_to_work[i, i] == 0.0:
            element_is_zero()

    matrix_b = get_matrix_b(a_to_work)

    print(f'B:\n {matrix_b}')

    if not (min(get_all_norms(matrix_b)[:2]) < 1):
        print("Одно из достаточных условий сходимости не выполняется!")
        # quit(0)

    current_iteration_x = numpy.zeros((n, 1))
    count_of_iterations = 0
    delta_one = epsilon
    delta_second = epsilon

    while delta_one + delta_second > epsilon:
        previous_iteration_x = current_iteration_x.copy()
        for i in range(n):
            s = 0
            for j in range(n):
                s += a_to_work[i, j] * current_iteration_x[j]
            s -= b_to_work[i]
            current_iteration_x[i] = current_iteration_x[i] - s / a_to_work[i, i]

        delta_one = numpy.absolute((current_iteration_x - previous_iteration_x)).max()
        delta_second = numpy.absolute((a_to_work.dot(current_iteration_x) - b_to_work)).max()

        count_of_iterations += 1

    print(f"\nКоличество итераций, при расчете методом Зейделя: {count_of_iterations}")
    print(f'x: {current_iteration_x.T}\n\n')


def main():
    a_to_work, b_to_work = initial()
    output_initial(a_to_work.copy(), b_to_work.copy())
    method_of_simple_iterations(a_to_work.copy(), b_to_work.copy())
    method_of_seidel(a_to_work.copy(), b_to_work.copy())


if __name__ == '__main__':
    main()
