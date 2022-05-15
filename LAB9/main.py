import numpy as np

count_of_dots = 10 ** 3
eps = 10 ** -3


# # EXAMPLE 1 #
# def yder(x, y): return x
#
#
# y0 = 0
# left_border, right_border = -1, 0
#
#
# def ans(x): return (1 / 2) * x ** 2
#
#
# xplot = np.linspace(left_border, right_border, count_of_dots)
# yplot = [ans(x_dot) for x_dot in xplot]
#
#
# # EXAMPLE 1 #


# # EXAMPLE 2 #
# def yder(x, y): return x + y
#
#
# y0 = 0
# left_border, right_border = -1, 1
#
#
# def ans(x): return -x - 1 + np.exp(x)
#
#
# xplot = np.linspace(left_border, right_border, count_of_dots)
# yplot = [ans(x_dot) for x_dot in xplot]
#
#
# # EXAMPLE 2 #

# MAINT TASK #
m, a = 1.0, 0.7


def yder(x, y): return (a * (1 - y ** 2)) / ((1 + m) * x ** 2 + y ** 2 + 1)


y0 = 0
left_border, right_border = 0, 1


# MAIN TASK #


def euler_method(x_dot, N):
    y_dots = [y0]
    h = x_dot / N
    for i in range(N):
        x = i * h
        y = y_dots[-1]
        y_dots += [y + h * yder(x, y)]
    return y_dots


def runge_kutta_method(x_dot, N):
    y_dots = [y0]
    h = x_dot / N
    for i in range(N):
        x = i * h
        y = y_dots[-1]
        K1 = h * yder(x, y)
        K2 = h * yder(x + h / 2, y + K1 / 2)
        K3 = h * yder(x + h / 2, y + K2 / 2)
        K4 = h * yder(x + h, y + K3)
        y_dots += [y + 1 / 6 * (K1 + 2 * K2 + 2 * K3 + K4)]
    return y_dots


def get_value_at_chosed_point(method, x):
    n = 1
    while True:
        olddots, newdots = method(x, n), method(x, 2 * n)
        if max(abs(newdots[2 * i] - olddots[i]) for i in range(n + 1)) < eps:
            return newdots[-1], 2 * n

        else:
            n *= 2


def create_dots(method, x_dots):
    maxn = 0
    midn = []
    for x in x_dots:
        y, n = get_value_at_chosed_point(method, x)
        maxn = max(maxn, n)
        midn += [n]
    midn = sum(midn) / len(x_dots)
    return midn, maxn


print("Методы Эйлера и Рунге-Кутта\n")

print("Точки для вычислений: ", count_of_dots)
print("Точность: ", eps)

x_dots = [left_border + (right_border - left_border) / count_of_dots * i for i in range(count_of_dots + 1)]

midn, maxn = create_dots(euler_method, x_dots)
print('\n\nМЕТОД ЭЙЛЕРА:')
print('Количество отрезков:')
print(f'Среднее: {midn}')
print(f'Максимальное: {maxn}')

midn, maxn = create_dots(runge_kutta_method, x_dots)
print('\n\nМЕТОД РУНГЕ-КУТТА:')
print('Количество отрезков:')
print(f'Среднее: {midn}')
print(f'Максимальное: {maxn}')
