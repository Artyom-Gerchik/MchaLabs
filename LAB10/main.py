import time

import matplotlib.pyplot as plt
import numpy as np


def runge_kutta_method(function, n, h, x, y):
    yn = 0
    for i in range(n):
        k1 = h * (function(x, y))
        k2 = h * (function((x + h / 2), (y + k1 / 2)))
        k3 = h * (function((x + h / 2), (y + k2 / 2)))
        k4 = h * (function((x + h), (y + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yn = y + k
        y = yn
        x = x + h

    return x, yn


def adams_bash(funct, n, t0, tn, y0):
    h = abs(tn - t0) / n
    t = np.linspace(t0, tn, n + 1)
    y = np.zeros(n + 1)
    y[0:3] = runge_kutta_method(funct, 2, t0, t0 + 2 * h, y0)[1]
    K1 = funct(t[1], y[1])
    K2 = funct(t[0], y[0])
    for i in range(2, n):
        K3 = K2
        K2 = K1
        K1 = funct(t[i], y[i])
        y[i + 1] = y[i] + h * (23 * K1 - 16 * K2 + 5 * K3) / 12

    return tn, y


# def function(x, y):
#     return x + y


def function(x, y):
    return 0.7 * (1 - y ** 2) / (2 * x ** 2 + y ** 2 + 1)


print("\nМетод Рунге-Кутта:")
t = time.perf_counter()
runge_kutta = runge_kutta_method(function, 100000, 0.00001, 0, 0)
t_w = time.perf_counter() - t
print(f"В точке ({round(runge_kutta[0])}) имеет значение: {'%.6f' % runge_kutta[1]}")
print(f"Время: {'%.6f' % t_w}")

print("\nМетод Адамса:")
t = time.perf_counter()
adams_bash = adams_bash(function, 100000, 0, 1, 0)
t_w = time.perf_counter() - t
print(f"В точке ({adams_bash[0]}) имеет значение:\n {adams_bash[1]}")
print(f"Время: {'%.6f' % t_w}")

x_str = [0]
for i in range(100000, 0, -1):
    x_str.append(x_str[100000 - i] + 1 / 100000)
plt.plot(x_str, adams_bash[1], "xb")
plt.grid(True)
plt.axis("equal")
plt.show()
