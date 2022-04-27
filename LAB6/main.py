# import matplotlib.pyplot as plt
import numpy as np


def input_values():
    x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    p = [0.0, 0.41, 0.79, 1.13, 1.46, 1.76, 2.04, 2.3, 2.55, 2.79, 3.01]
    k = 7
    m = 3.5
    y = [(p[i] + ((-1) ** k) * m) for i in range(len(x))]
    dots = list(zip(x, y))

    # dots = [(-1, 0), (0, 1), (1, 0)]
    # dots = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 6)]

    return dots


dots = input_values()
(x, y) = map(list, zip(*dots))


def calculate_lagrange_polynom(dots):
    n = len(dots)
    (x, y) = map(list, zip(*dots))
    polynom = np.poly1d([0])
    for i in range(n):
        p = np.poly1d([1])
        for j in range(n):
            if j != i:
                p *= np.poly1d([1, -x[j]]) / (x[i] - x[j])
        polynom += y[i] * p
    return polynom


def calculate_divided_differences(xs):
    n = len(xs)
    diffs = [[None for j in range(n - i)] for i in range(n)]
    for i in range(n):
        diffs[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            diffs[i][j] = ((diffs[i][j - 1] - diffs[i + 1][j - 1]) / (xs[i] - xs[i + j]))
    return diffs


def calculate_inaccurany(xs, xdot):
    n = len(xs)
    diffs = calculate_divided_differences(xs)
    maxdiff = 0.0
    for i in range(len(diffs)):
        for j in range(len(diffs[i])):
            maxdiff = max(maxdiff, abs(diffs[i][j]))
    w = 1
    for i in range(n):
        w *= xdot - xs[i]
    f = 1
    for i in range(1, n + 1 + 1):
        f *= i
    R = maxdiff * w / f
    return R


def calculate_newton_polynom(dots):
    n = len(dots)
    (x, y) = map(list, zip(*dots))

    diffs = calculate_divided_differences(x)

    polynom = np.poly1d([0])
    for i in range(n):
        p = np.poly1d([1])
        for j in range(i):
            p *= np.poly1d([1, -x[j]])
        polynom += p * diffs[0][i]

    return polynom


def calculate_least_square_polynom(dots, m=None):
    n = len(dots) - 1
    if m is None:
        m = n
    assert 0 <= m <= n
    return np.poly1d(np.polyfit(*map(list, zip(*dots)), m))


def main():
    print("Интерполяционные многочлены \n")

    print("Исходные (x,y) = ", dots, '\n')

    lagrange = calculate_lagrange_polynom(dots)
    print("Многочлен Лагранджа = ")
    print(lagrange, '\n')

    newton = calculate_newton_polynom(dots)
    print("Многочлен Ньютона = ")
    print(newton, '\n')

    squares = calculate_least_square_polynom(dots)
    print("Многочлен наименьших квадратов = ")
    print(squares, '\n')

    dot_to_calculate = 0.47
    print(f"Многочлен Лагранджа от ({dot_to_calculate}) = ", lagrange(dot_to_calculate))
    print(f"Многочлен Ньютона от ({dot_to_calculate}) = ", newton(dot_to_calculate))
    print(f"Многочлен наименьших квадратов от ({dot_to_calculate}) = ", squares(dot_to_calculate))
    print(f"|Лаграндж({dot_to_calculate}) - Ньютон({dot_to_calculate})| =",
          abs(lagrange(dot_to_calculate) - newton(dot_to_calculate)))
    print(f"|Лаграндж({dot_to_calculate}) - Наименьшие квадраты({dot_to_calculate})| =",
          abs(lagrange(dot_to_calculate) - squares(dot_to_calculate)))
    print(f"|Ньютон({dot_to_calculate}) - Наименьшие квадраты({dot_to_calculate})| =",
          abs(newton(dot_to_calculate) - squares(dot_to_calculate)))

    print(f"Погрешность от ({dot_to_calculate}) = ", calculate_inaccurany(x, dot_to_calculate))

    # plotdots = 10 ** 4
    #
    # plt.plot(x, y, 'og')
    #
    # xplot = np.linspace(min(x), max(x), plotdots)
    #
    # yplot = [squares(xdot) for xdot in xplot]
    # plt.plot(xplot, yplot, 'y')
    #
    # yplot = [lagrange(xdot) for xdot in xplot]
    # plt.plot(xplot, yplot, 'b')
    #
    # yplot = [newton(xdot) for xdot in xplot]
    # plt.plot(xplot, yplot, 'r--')
    #
    # plt.show()


if __name__ == '__main__':
    main()
