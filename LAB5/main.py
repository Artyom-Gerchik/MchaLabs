import numpy

print("Собственные вектора и собственные значения\n")

EPS = 10.0 ** -4

numpy.set_printoptions(suppress=True, precision=4, floatmode="fixed")


def input_values():
    C = numpy.array([
        [0.2, 0.0, 0.2, 0.0, 0.0],
        [0.0, 0.2, 0.0, 0.2, 0.0],
        [0.2, 0.0, 0.2, 0.0, 0.2],
        [0.0, 0.2, 0.0, 0.2, 0.0],
        [0.0, 0.0, 0.2, 0.0, 0.2]
    ])
    D = numpy.array([
        [2.33, 0.81, 0.67, 0.92, -0.53],
        [0.81, 2.33, 0.81, 0.67, 0.92],
        [0.67, 0.81, 2.33, 0.81, 0.92],
        [0.92, 0.67, 0.81, 2.33, -0.53],
        [-0.53, 0.92, 0.92, -0.53, 2.33]
    ])
    matrix_A = 7 * C + D

    return matrix_A


matrix_A = input_values()
n = len(matrix_A)
print(f'Матрица А:\n {matrix_A}\n')

if abs((matrix_A - matrix_A.T) ** 2).sum() > EPS:
    raise ValueError("Ошибка, матрица А неправильная")

count_of_iterations = 0

ans_matrix_V = numpy.eye(n)
while True:
    count_of_iterations += 1
    max_elem = (0, 1)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix_A[i][j]) > abs(matrix_A[max_elem]):
                max_elem = (i, j)
    (i, j) = max_elem
    if matrix_A[i][i] == matrix_A[j][j]:
        p = numpy.pi / 4
    else:
        p = 2 * matrix_A[i][j] / (matrix_A[i][i] - matrix_A[j][j])
    cos = numpy.cos(1 / 2 * numpy.arctan(p))
    sin = numpy.sin(1 / 2 * numpy.arctan(p))
    V = numpy.eye(n)

    V[i][i] = cos
    V[i][j] = -sin
    V[j][i] = sin
    V[j][j] = cos

    matrix_A = V.T @ matrix_A @ V  # matrix multiply
    ans_matrix_V = ans_matrix_V @ V  # matrix multiply

    if abs(matrix_A - numpy.diag(numpy.diag(matrix_A))).sum() < EPS:
        ansW = numpy.diag(matrix_A)
        break


def normalization(W, V):
    V = numpy.array([(-i if i[0] < 0 else i) for i in V.T]).T
    (W, V) = list(zip(*(sorted(list(zip(W, V.T)), key=lambda t: t[0]))))
    W = numpy.array(W)
    V = numpy.array(V).T
    return (W, V)


matrix_A = input_values()

print('Метод вращений Якоби: ')
(W, V) = normalization(ansW, ans_matrix_V)
print(f"Собственные значения = {W}")
print(f"Собственные векторы = \n {V}")
print("Потребовалось итераций: =", count_of_iterations)

print("\nСтепенной метод:")
matrix_A = input_values()
r = numpy.ones(len(matrix_A))
count_of_iterations = 0
while True:
    count_of_iterations += 1
    old_u = (r.T @ matrix_A @ r) / (r.T @ r)
    r = (matrix_A @ r) / numpy.sqrt(sum((matrix_A @ r) ** 2))
    u = (r.T @ matrix_A @ r) / (r.T @ r)
    if abs(u - old_u) < EPS:
        break
print(f'Максимальное, по модулю, собственное значение = {u:.4f}')
print(f'Максимальный собственный вектор = {r}')
print("Количество итераций = ", count_of_iterations)
