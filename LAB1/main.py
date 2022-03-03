import numpy


def initial():
    print("\n")
    print("Решение систем линейных алгебраических уравнений (СЛАУ) методом Гаусса и с помощью его модификаций")
    print("Вариант 7 \n")


def inputValues():
    numpy.set_printoptions(suppress=True, precision=4, floatmode="fixed")
    k = 7
    b = numpy.array([[4.2], [4.2], [4.2], [4.2], [4.2]])
    C = numpy.array([
        [0.2, 0.0, 0.2, 0.0, 0.0],
        [0.0, 0.2, 0.0, 0.2, 0.0],
        [0.2, 0.0, 0.2, 0.0, 0.2],
        [0.0, 0.2, 0.0, 0.2, 0.0],
        [0.0, 0.0, 0.2, 0.0, 0.2]
    ])
    D = numpy.array([
        [2.33, 0.81, 0.67, 0.92, -0.53],
        [-0.53, 2.33, 0.81, 0.67, 0.92],
        [0.92, -0.53, 2.33, 0.81, 0.67],
        [0.67, 0.92, -0.53, 2.33, 0.81],
        [0.81, 0.67, 0.92, -0.53, 2.33]
    ])
    A = k * C + D

    # # test example 1 normal
    # A = numpy.array([
    #     [2.0, 1.0],
    #     [1.0, -2.0]
    # ])
    # b = numpy.array([[3.0], [1.0]])
    # # test example 1 normal

    # # test example 2 normal
    # A = numpy.array([
    #     [4.0, 1.0, 1.0],
    #     [1.0, 4.0, 1.0],
    #     [1.0, 1.0, 4.0]
    # ])
    # b = numpy.array([[1.0], [1.0], [1.0]])
    # # test example 2 normal

    # # test example 3 with unlimited solutions but 0 on main matrix isn't square
    # A = numpy.array([
    #     [1.0, 2.0, 1.0, 1.0, 3.0, 1.0],
    #     [1.0, 2.0, 1.0, 2.0, 1.0, -1.0],
    #     [1.0, 2.0, 1.0, -1.0, 5.0, -1.0],
    #     [1.0, 2.0, 1.0, -2.0, -4.0, 4.0]
    # ])
    # b = numpy.array([[7.0], [1.0], [2.0], [-1.0]])
    # # test example 3 with unlimited solutions but 0 on main matrix isn't square

    # # test example 4 with unlimited solutions but 0 on main matrix isn't square
    # A = numpy.array([
    #     [1.0, 2.0, -3.0, 5.0],
    #     [1.0, 3.0, -13.0, 22.0],
    #     [3.0, 5.0, 1.0, -2.0],
    #     [2.0, 3.0, 4.0, -7.0]
    # ])
    # b = numpy.array([[1.0], [-1.0], [5.0], [4.0]])
    # # test example 4 with unlimited solutions but 0 on main matrix isn't square

    # # test example 5 no solutions but 0 on main matrix isn't square
    # A = numpy.array([
    #     [2.0, -1.0, 3.0],
    #     [2.0, -1.0, -1.0],
    #     [4.0, -2.0, 6.0],
    #     [6.0, 8.0, -7.0]
    # ])
    # b = numpy.array([[1.0], [-2.0], [0.0], [2.0]])
    # # test example 5 no solutions but 0 on main matrix isn't square

    # # test example 6
    # A = numpy.array([
    #     [3.0, 2.0, -5.0],
    #     [2.0, -1.0, 3.0],
    #     [1.0, 2.0, -1.0]
    # ])
    # b = numpy.array([[-1.0], [13.0], [9.0]])
    # # test example 6

    # # test example 6 unlimited of solutions
    # A = numpy.array([
    #     [1.0, 2.0, -3.0, 5.0],
    #     [1.0, 3.0, -13.0, 22.0],
    #     [3.0, 5.0, 1.0, -2.0],
    #     [2.0, 3.0, 4.0, -7.0]
    # ])
    # b = numpy.array([[1.0], [-1.0], [5.0], [4.0]])
    # # test example 6 unlimited of solutions

    A[1][2] = 10.8100
    A[4][4] = 5.7300

    print("b:\n", b)
    print("\nC:\n", C)
    print("\nD:\n", D)
    print("\nA:\n", A)

    return A, b


def elementIsZero():
    print("Element == 0!")
    quit()


def checkIsUnlimitedOfSolutions(A, n):
    for i in range(n - 1, 0, -1):
        counter = 0
        for j in range(A.shape[1]):
            if numpy.around(A[i][j], decimals=4) != 0:
                counter += 1
        if counter == 1:
            return
        else:
            print('Бесконечное множество решений!')
            quit()


def checkIsNoSolutions(A, n, b):
    for i in range(n - 1, 0, -1):
        counter = 0
        for j in range(A.shape[1]):
            if numpy.around(A[i][j], decimals=4) != 0:
                counter += 1
        if counter == 0 and b[i] != 0:
            print('Нет решений!')
            quit()
        else:
            return


def baseMethod(A, b):
    if A.shape[0] != A.shape[1]:
        print("Матрица не квадратная!")
        quit()
    n = min(A.shape[0], A.shape[1])  # length of array

    for k in range(n):  # n is stop
        for i in range(k + 1, n):  # (start, stop, step)
            # if k == n - 2:
            #     checkIsUnlimitedOfSolutions(A, n)

            if A[k, k] == 0.0:
                checkIsNoSolutions(A, n, b)
                checkIsUnlimitedOfSolutions(A, n)
                elementIsZero()
            q = A[i, k] / A[k, k]

            for j in range(n):
                A[i, j] -= A[k, j] * q
            b[i] -= b[k] * q

    print("\nПреобразованный вектор b, схемой единственного деления:\n", b)
    print("\nПреобразованная матрица A, схемой единственного деления:\n", A)
    checkIsNoSolutions(A, n, b)
    checkIsUnlimitedOfSolutions(A, n)
    x = numpy.zeros((n, 1))

    for i in range(n - 1, -1, -1):  # (start, stop, step)
        for j in range(i + 1, n):
            b[i] -= A[i, j] * x[j]
        if A[i, i] == 0.0:
            elementIsZero()
        x[i] = b[i] / A[i, i]

    output("\nСхема единственного деления", x)


def partialChoiceMethod(A, b):
    if A.shape[0] != A.shape[1]:
        print("Матрица не квадратная!")
        quit()
    n = min(A.shape[0], A.shape[1])  # length of array

    for k in range(n):  # n is stop

        maxEl = k

        for temp in range(k, n):
            if abs(A[temp, k]) > abs(A[maxEl, k]):
                maxEl = temp
        A[[k, maxEl]] = A[[maxEl, k]]
        b[[k, maxEl]] = b[[maxEl, k]]

        for i in range(k + 1, n):  # (start, stop, step)

            if A[k, k] == 0.0:
                checkIsUnlimitedOfSolutions(A, n)
                elementIsZero()
            q = A[i, k] / A[k, k]

            for j in range(n):
                A[i, j] -= A[k, j] * q
            b[i] -= b[k] * q

    print("\nПреобразованный вектор b, схемой частичного выбора:\n", b)
    print("\nПреобразованная матрица A, схемой частичного выбора:\n", A)
    checkIsNoSolutions(A, n, b)
    checkIsUnlimitedOfSolutions(A, n)
    x = numpy.zeros((n, 1))

    for i in range(n - 1, -1, -1):  # (start, stop, step)
        for j in range(i + 1, n):
            b[i] -= A[i, j] * x[j]
        if A[i, i] == 0.0:
            # checkIsUnlimitedOfSolutions(A, n)
            elementIsZero()
        x[i] = b[i] / A[i, i]

    output("\nСхема частичного выбора ", x)


def fullChoiceMethod(A, b):
    if A.shape[0] != A.shape[1]:
        print("Матрица не квадратная!")
        quit()
    n = min(A.shape[0], A.shape[1])  # length of array

    for i in range(n):

        maxEl = (i, 0)

        for k in range(i, n):
            for j in range(n):
                if abs(A[k, j]) > abs(A[maxEl[0], maxEl[1]]):
                    maxEl = (k, j)  # position of max element
        A[[i, maxEl[0]]] = A[[maxEl[0], i]]  # swapping rows
        b[[i, maxEl[0]]] = b[[maxEl[0], i]]

        for k in range(n):
            if k != i:  # if == - main element
                curElForQ02 = A[i, i]
                if curElForQ02 == 0.0:
                    checkIsUnlimitedOfSolutions(A, n)
                    elementIsZero()
                curElForQ01 = A[k, i]
                q = curElForQ01 / curElForQ02
                for j in range(n):
                    A[k, j] -= A[i, j] * q
                b[k] -= b[i] * q

    print("\nПреобразованный вектор b, схемой полного выбора:\n", b)
    print("\nПреобразованная матрица A, схемой полного выбора:\n", A)
    checkIsNoSolutions(A, n, b)
    checkIsUnlimitedOfSolutions(A, n)
    x = numpy.zeros((n, 1))

    for j in range(n):

        maxEl = (j, 0)

        for i in range(n):
            if abs(A[i, j]) > abs(A[maxEl[0], maxEl[1]]):
                maxEl = (i, j)
        if A[maxEl[0], maxEl[1]] == 0.0:
            elementIsZero()
        x[j] = b[maxEl[0]] / A[maxEl[0], maxEl[1]]

    output("\nСхема полного выбора ", x)


def output(whichMethod, x):
    print(whichMethod)
    print("x:")
    print(x.T)

    numpy.set_printoptions(suppress=True, precision=10, floatmode="fixed")
    print("x:")
    print(x.T)
    numpy.set_printoptions(suppress=True, precision=4, floatmode="fixed")


def main():
    initial()
    AToWork, bToWork = inputValues()
    baseMethod(AToWork.copy(), bToWork.copy())
    partialChoiceMethod(AToWork.copy(), bToWork.copy())
    fullChoiceMethod(AToWork.copy(), bToWork.copy())

    # print('\n\n\n test')
    # matrix = numpy.array([
    #     [3.0, 2.0, -5.0],
    #     [0.0, -1.0, 3.0],
    #     [0.0, 2.0, -1.0]
    # ])
    # checkIsUnlimitedOfSolutions(matrix, matrix.shape[0])
    # matrix = numpy.array([
    #     [3.0, 2.0, -5.0],
    #     [0.0, -1.0, 3.0],
    #     [0.0, 0.0, 0.0]
    # ])
    # b = numpy.array([[4.2], [4.2], [4.2]])
    # checkIsNoSolutions(matrix, matrix.shape[0], b)


if __name__ == '__main__':
    main()
