# # EXAMPLE 1 #
# def function(x, y):
#     return x
# # EXAMPLE 1 #


# # EXAMPLE 2 #
# def function(x, y):
#     return x + y
# # EXAMPLE 2 #

# # MAIN TASK #
def function(x, y):
    return (0.7 * (1 - y ** 2)) / ((1 + 1.0) * x ** 2 + y ** 2 + 1)


# # MAIN TASK #


def euler_method(func, n, h, x, y):
    for i in range(n):
        y += h * func(x, y)
        x += h
    return x, y


def runge_kutta_method(func, n, h, x, y):
    yn = 0
    for i in range(n):
        k1 = h * (func(x, y))
        k2 = h * (func((x + h / 2), (y + k1 / 2)))
        k3 = h * (func((x + h / 2), (y + k2 / 2)))
        k4 = h * (func((x + h), (y + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yn = y + k
        y = yn
        x = x + h

    return x, yn


print("Метод Эйлера и Рунге-Кутта(4-го порядка)\n")

n = 10
h = 0.1
x = 0
y = 0

for i in range(7):
    print('++++++++++++++++++++++++++++++')
    print()
    print(f'h = {h:.10f}')
    print()

    point, answer = euler_method(function, n, h, x, y)
    print('Метод Эйлера')
    print(f'В точке {round(point)} Значение = ', "%.6f" % answer)

    print()
    print('Метод Рунге-Кутта(4-го порядка)')
    point, answer = runge_kutta_method(function, n, h, x, y)
    print(f'В Точке {round(point)} Значение = ', "%.6f" % answer)

    n *= 10
    h /= 10
    print()
    print('++++++++++++++++++++++++++++++')
