import math

import numpy as np

np.random.seed(42)

# MAIN TASK #

L, R, DerXdot = 0, 2, 1


def f(x): return np.arctan(np.sqrt(x))


def F(x): return x * np.arctan(np.sqrt(x)) - np.sqrt(x) + np.arctan(np.sqrt(x))


def fd1(x): return 1 / (2 * (x ** (1 / 2)) * (1 + x))


def fd2(x): return -1 / (4 * x ** (3 / 2) * (1 + x)) - 1 / (2 * (x ** (1 / 2)) * (1 + x) ** 2)


def fd3(x): return 3 / (8 * x ** (5 / 2) * (1 + x)) + 1 / (2 * x ** (3 / 2) * (1 + x) ** 2) + 1 / (
        (x ** (1 / 2)) * (1 + x) ** 3)


def fd4(x): return -15 / (16 * x ** (7 / 2) * (1 + x)) - 9 / (8 * x ** (5 / 2) * (1 + x) ** 2) - 3 / (
        2 * x ** (3 / 2) * (1 + x) ** 3) - 3 / ((x ** (1 / 2)) * (1 + x) ** 4)


M2deLR, M4deLR = 0, 0

# MAIN TASK #

# # First Example #
#
#
# L, R, DerXdot = -1, 1, 0
#
#
# def f(x): return np.exp(x)
#
#
# def F(x): return np.exp(x)
#
#
# def fd1(x): return np.exp(x)
#
#
# def fd2(x): return np.exp(x)
#
#
# def fd3(x): return np.exp(x)
#
#
# def fd4(x): return np.exp(x)
#
#
# M2deLR, M4deLR = 0, 0
#
# # First Example #

# # Second Example #
#
#
# L, R, DerXdot = -1, 1, 1 / 2
#
#
# def f(x): return np.sqrt(1 - x ** 2)
#
#
# def F(x): return (1 / 2) * x * np.sqrt(1 - x ** 2) + (1 / 2) * np.arcsin(x)
#
#
# def fd1(x): return -x / np.sqrt(1 - x ** 2)
#
#
# def fd2(x): return -x ** 2 / (-x ** 2 + 1) ** (3 / 2) - 1 / np.sqrt(-x ** 2 + 1)
#
#
# def fd3(x): return -3 * x ** 3 / (-x ** 2 + 1) ** (5 / 2) - 3 * x / (-x ** 2 + 1) ** (3 / 2)
#
#
# def fd4(x): return -15 * x ** 4 / (-x ** 2 + 1) ** (7 / 2) - 18 * x ** 2 / (-x ** 2 + 1) ** (5 / 2) - 3 / (
#             -x ** 2 + 1) ** (3 / 2)
#
#
# M2deLR, M4deLR = 0, 0
#
# # Second Example #


IntEps = 0.000001
IntFormatString = "{:.7f}"

DerEps = 0.01
DerFormatString = "{:.3f}"


def diff_first(f, x, d):
    return (f(x + d) - f(x - d)) / (2 * d)


def diff_first_via_estimation(f, x):
    M2 = abs(fd2(x))
    df = 2 * DerEps / M2
    M3 = abs(fd3(x))
    ds = (6 * DerEps / M3) ** (1 / 2)
    return diff_first(f, x, min(df, ds))


def diff_first_via_ten_in_minus_5(f, x):
    d = 10.0 ** -5
    return diff_first(f, x, d)


def derivative_second(f, x, d):
    return (f(x + d) - 2 * f(x) + f(x - d)) / (d ** 2)


def derivative_second_via_estimation(f, x):
    M4 = abs(fd4(x))
    d = (12 * DerEps / M4) ** (1 / 2)
    return derivative_second(f, x, d)


def derivative_second_via_ten_in_minus_4(f, x):
    d = 10.0 ** -4
    return derivative_second(f, x, d)


print("Численное дифференцирование и интегрирование функций \n")

print()

print("Первая производная = " + DerFormatString.format(fd1(DerXdot)))
print("Приближенное    = " + DerFormatString.format(diff_first_via_estimation(f, DerXdot)))
print("Округленное   = " + DerFormatString.format(diff_first_via_ten_in_minus_5(f, DerXdot)))
print()

print("Вторая Производная = " + DerFormatString.format(fd2(DerXdot)))
print("Приближенное     = " + DerFormatString.format(derivative_second_via_estimation(f, DerXdot)))
print("Округленное    = " + DerFormatString.format(derivative_second_via_ten_in_minus_4(f, DerXdot)))
print()


def integral_via_middle_rectangles(f, L, R, N):
    h = (R - L) / N
    x = L + h / 2
    s = 0.0
    while x < R:
        s += f(x) * h
        x += h
    return s


def integral_via_trapezoids(f, L, R, N):
    h = (R - L) / N
    x = L + h / 2
    s = 0.0
    while x < R:
        s += ((f(x - h / 2) + f(x + h / 2)) / 2) * h
        x += h
    return s


def integral_via_simpson(f, L, R, N):
    h = (R - L) / N
    x = L + h / 2
    s = 0.0
    while x < R:
        fa = f(x - h / 2)
        fm = f(x)
        fb = f(x + h / 2)
        s += (fa + 4 * fm + fb) * h / 6
        x += h
    return s


def integral_via_random_segments(method, f, L, R):
    LeftCoeff, RightCoeff = 1 / 3, 1 / 2
    h_prev = R - L
    ans_prev = method(f, L, R, 1)
    while True:
        h_new = h_prev * (LeftCoeff + (RightCoeff - LeftCoeff) * np.random.rand())
        N = math.floor((R - L) / h_new)
        M = L + h_new * N
        ans_new = method(f, L, M, N) + method(f, M, R, 1)
        if abs(ans_new - ans_prev) < IntEps:
            print("\nN =", N)
            return ans_new
        ans_prev = ans_new
        h_prev = h_new


def integral_via_middle_rectangles_via_estimation(f, L, R):
    if M2deLR > 0.0:
        M2 = M2deLR
        h = (24 * IntEps / (R - L) / M2) ** (1 / 2)
        N = np.ceil((R - L) / h)
        return integral_via_middle_rectangles(f, L, R, N)
    else:
        return integral_via_random_segments(integral_via_middle_rectangles, f, L, R)


def integral_via_trapezoids_via_estimation(f, L, R):
    if M2deLR > 0.0:
        M2 = M2deLR
        h = (12 * IntEps / (R - L) / M2) ** (1 / 2)
        N = np.ceil((R - L) / h)
        return integral_via_trapezoids(f, L, R, N)
    else:
        return integral_via_random_segments(integral_via_trapezoids, f, L, R)


def integral_via_simpson_via_estimation(f, L, R):
    if M4deLR > 0.0:
        M4 = M4deLR
        h = (180 * IntEps / (R - L) / M4) ** (1 / 4)
        N = np.ceil((R - L) / h)
        return integral_via_simpson(f, L, R, N)
    else:
        return integral_via_random_segments(integral_via_simpson, f, L, R)


print()
intprecised = F(R) - F(L)
print("Интеграл =            " + IntFormatString.format(intprecised))


def delta(intappr): return np.ceil(abs(intappr - intprecised) * (1 / (IntEps / 10))) * (IntEps / 10)


intappr = integral_via_middle_rectangles_via_estimation(f, L, R)
print(("Метод Средних Прямоугольников = " + IntFormatString + " | дельта = " + IntFormatString).format(intappr,
                                                                                                       delta(intappr)))
intappr = integral_via_trapezoids_via_estimation(f, L, R)
print(("Метод Трапеций = " + IntFormatString + " | дельта = " + IntFormatString).format(intappr, delta(intappr)))
intappr = integral_via_simpson_via_estimation(f, L, R)
print(("Метод Симпсона = " + IntFormatString + " | дельта = " + IntFormatString).format(intappr, delta(intappr)))
print()
