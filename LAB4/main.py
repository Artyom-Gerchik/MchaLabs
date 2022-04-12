import numpy
import sympy

print("Системы нелинейных уравнений \n")

m = 0.3
a = 0.6

EPS = 10.0 ** -4

(x, y) = sympy.symbols("x y")
eq1 = sympy.tan(x * y + m) - x
eq2 = a * (x ** 2) + 2 * (y ** 2) - 1


# Исходное 1
def val1(x, y):
    return numpy.tan(x * y + m) - x


# Исходное 2
def val2(x, y):
    return a * (x ** 2) + 2 * (y ** 2) - 1


# Выражение из Исходное 1
def eqx(x, y):
    return numpy.tan(x * y + m)


# Выражение из Исходное 2
def eqy(x, y):
    return numpy.sqrt((1 - a * (x ** 2)) / 2)


# Матрица Якоби
def jakobi_matrix(x, y):
    return numpy.array([
        [(1 + numpy.tan(x * y + m) ** 2) * y - 1, (1 + numpy.tan(x * y + m) ** 2) * x],
        [2 * a * x, 4 * y]
    ])


count_of_iterations = 0


def method_of_simple_iterations(x0, y0):
    global count_of_iterations
    count_of_iterations = 0
    (x, y) = (x0, y0)
    while True:
        count_of_iterations += 1
        oldx = x
        oldy = y
        x = eqx(x, y)
        y = eqy(x, y)

        if not (numpy.isfinite(x) and numpy.isfinite(y)):
            raise RuntimeError("Расходится")
        if max(abs(x - oldx), abs(y - oldy)) < EPS:
            return x, y


def method_of_newton(x0, y0):
    global count_of_iterations
    count_of_iterations = 0
    (x, y) = (x0, y0)
    while True:
        count_of_iterations += 1
        w = jakobi_matrix(x, y)
        f = numpy.array([[val1(x, y)], [val2(x, y)]])
        deltas = numpy.linalg.solve(w, -f)
        x += deltas[0][0]
        y += deltas[1][0]
        if not (numpy.isfinite(x) and numpy.isfinite(y)):
            raise RuntimeError("Расходится")
        if max(abs(deltas)) < EPS:
            return x, y


def main():
    print("Первое уравнение: ", eq1, "= 0")
    print("Второе уравнение:", eq2, "= 0")
    print()

    x0 = 1.0
    y0 = 0.5
    print("Начальноe приближение: ", (x0, y0))
    print()

    global count_of_iterations
    try:
        (x, y) = method_of_simple_iterations(x0, y0)
        print("(x, y) = ({:.4f}, {:.4f})".format(x, y))
        print(f"Используя {method_of_simple_iterations.__name__}, заняло {count_of_iterations} итераций")
    except Exception as ex:
        print(f"ашыпка в : {method_of_simple_iterations.__name__}")
    print()

    try:
        (x, y) = method_of_newton(x0, y0)
        print("(x, y) = ({:.4f}, {:.4f})".format(x, y))
        print(f"Используя {method_of_newton.__name__}, заняло {count_of_iterations} итераций")
    except Exception as ex:
        print(f"ашыпка в : {method_of_newton.__name__}")
    print()


if __name__ == '__main__':
    main()
