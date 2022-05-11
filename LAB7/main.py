import numpy


def input_values():
    ##### EXAMPLES #####

    # def function(x): return numpy.tan(x)
    #
    # left_border, dots_count, right_border = 0, 6, 1.5

    # def function(x): return numpy.sinh(x)
    #
    # left_border, dots_count, right_border = 0, 6, 2

    ##### EXAMPLES #####

    def function(x): return numpy.cosh(x)  # func of my var

    left_border, dots_count, right_border = 0, 6, 2

    dots = []
    for i in range(dots_count):
        x = left_border + (right_border - left_border) * i / (dots_count - 1)
        y = function(x)
        dots += [(x, y)]

    return dots, function


print("Интерполяция сплайнами\n")

dots, f = input_values()
(x, y) = map(list, zip(*dots))
print("(x,y) =", dots, '\n')


def triple_diagonal_solve(A, b):
    A = A.copy()
    b = b.copy()
    n = len(A)

    A[0][1] /= A[0][0]
    for i in range(1, n - 1):
        A[i][i + 1] /= (A[i][i] - A[i][i - 1] * A[i - 1][i])

    b[0] /= A[0][0]
    for i in range(1, n):
        b[i] = (b[i] - A[i][i - 1] * b[i - 1]) / (A[i][i] - A[i][i - 1] * A[i - 1][i])

    x = numpy.zeros(n)
    x[-1] = b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = b[i] - A[i][i + 1] * x[i + 1]

    return x


def spline_method(dots):
    n = len(dots) - 1
    (x, y) = map(list, zip(*dots))

    h = [None]
    for i in range(1, n + 1):
        h += [x[i] - x[i - 1]]

    A = [[None] * (n) for i in range(n)]
    for i in range(1, n):
        for j in range(1, n):
            A[i][j] = 0.0
    for i in range(1, n - 1):
        A[i + 1][i] = h[i + 1]
    for i in range(1, n):
        A[i][i] = 2 * (h[i] + h[i + 1])
    for i in range(1, n - 1):
        A[i][i + 1] = h[i + 1]

    F = []
    for i in range(1, n):
        F += [3 * ((y[i + 1] - y[i]) / h[i + 1] - (y[i] - y[i - 1]) / h[i])]

    A = [A[i][1:] for i in range(len(A)) if i]

    c = triple_diagonal_solve(A, F)
    c = [0.0] + list(c) + [0.0]

    def evaluate(x_dot):
        for i in range(1, len(x)):
            if x[i - 1] <= x_dot <= x[i]:
                val = 0
                val += y[i]
                b = (y[i] - y[i - 1]) / h[i] + (2 * c[i] + c[i - 1]) * h[i] / 3
                val += b * (x_dot - x[i])
                val += c[i] * ((x_dot - x[i]) ** 2)
                d = (c[i] - c[i - 1]) / (3 * h[i])
                val += d * ((x_dot - x[i]) ** 3)
                return val
        return None

    def output():
        print("Кубический сплайн: ", '\n')
        for i in range(1, len(x)):
            val = 0
            b = (y[i] - y[i - 1]) / h[i] + (2 * c[i] + c[i - 1]) * h[i] / 3
            d = (c[i] - c[i - 1]) / (3 * h[i])
            print(x[i - 1], x[i], )
            print(numpy.poly1d([d, c[i], b, y[i]]), '\n')

    return evaluate, output


spl, cout = spline_method(dots)

cout()

x_dot = 1.0  # my variant & second example
# x_dot = 0.75  # first example
print(f"Исходная функция от ({x_dot}) = ", f(x_dot))
print(f"Кубический сплайн от ({x_dot}) = ", spl(x_dot))
print(f"Дельта от ({x_dot}) = ", abs(f(x_dot) - spl(x_dot)))
