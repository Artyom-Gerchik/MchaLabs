import numpy

numpy.random.seed(42)

EPS = 10.0 ** -4


def input():
    expr = numpy.poly1d([1.0, 38.4621, 364.594, 914.196])
    return expr


(f) = input()


def sturm_method(f):
    arr = []
    arr.append(f)
    arr.append(numpy.polyder(f))

    while True:
        fn = -numpy.polydiv(arr[-2], arr[-1])[1]
        if (fn.order > 0 or abs(fn[0]) > 0.0):
            arr.append(fn)
        else:
            break

    return arr


def n(stseq, x):
    if (abs(f(x)) < EPS):
        raise ValueError("Число в N() является корнем")
    ans = 0
    for i in range(1, len(stseq)):
        if (stseq[i](x) == 0.0):
            raise ValueError("Элемент ряда Штурма равен 0")
        if (stseq[i - 1](x) * stseq[i](x) < 0):
            ans += 1
    return ans


def get_bounds(f, a, b):
    if ((abs(f(a)) < EPS) or (abs(f(b)) < EPS)):
        raise ValueError("Границы содержат данный корень")
    if (n(sturm_system, a) - n(sturm_system, b) == 0):
        return []
    if (n(sturm_system, a) - n(sturm_system, b) > 1):
        while True:
            M = a + (b - a) / (1.5 + numpy.random.random())
            if (abs(f(M)) > EPS):
                break
        return get_bounds(f, a, M) + get_bounds(f, M, b)
    if (b - a < EPS):
        print("Границы маловаты")
    return [(a, b)]


sturm_system = sturm_method(f)

iters = 0


def binary_search(L, R):
    global iters
    iters += 1
    M = (L + R) / 2
    if (R - L < EPS):
        return M
    # if (abs(f(M)) < EPS):
    #    return M
    if (f(L) * f(M) <= 0):
        return binary_search(L, M)
    elif (f(R) * f(M) <= 0):
        return binary_search(M, R)
    else:
        raise RuntimeError("Ошибка в двоичном поиске")


def chord_method(L, R):
    global iters
    fder2 = numpy.polyder(f, 2)
    if (f(R) * fder2(R) > 0):
        (oldx, x) = (R, L)
    elif (f(L) * fder2(L) > 0):
        (oldx, x) = (L, R)
    else:
        raise ValueError("Границы плоховаты в методе Хорд")
    t = oldx
    while (abs(x - oldx) > EPS):
        iters += 1
        oldx = x
        x = x - f(x) * (t - x) / (f(t) - f(x))
        if (not (numpy.isfinite(x))):
            raise RuntimeError("Метод Хорд Ошибка")
    if ((x < L) or (R < x)):
        raise RuntimeError("Метод Хорд Ошибка")
    return x


def newton(L, R):
    global iters
    fder = numpy.polyder(f)
    fder2 = numpy.polyder(f, 2)
    if (f(L) * fder2(L) > 0):
        (oldx, x) = (R, L)
    elif (f(R) * fder2(R) > 0):
        (oldx, x) = (L, R)
    else:
        raise ValueError("Границы плоховаты в методе Ньютона")
    while (abs(x - oldx) > EPS):
        iters += 1
        oldx = x
        x = x - f(x) / fder(x)
        if (not (numpy.isfinite(x))):
            raise RuntimeError("Метод Ньютона Ошибка")
    if ((x < L) or (R < x)):
        raise RuntimeError("Метод Ньютона Ошибка")
    return x


numpy.set_printoptions(suppress=True, precision=4, floatmode="fixed")

bounds = get_bounds(f, -10, 10)


def test_method(method):
    global iters
    for i in range(len(bounds)):
        iters = 0
        try:
            str = method(*bounds[i])
            if (not str is None):
                str = "{:.4f}".format(str)
            print(f"{str} используя {method.__name__} заняло {iters} итераций)")
            break
        except Exception as ex:
            print("ERROR: {} - in {} method (with {} iterations)".format(ex, method.__name__, iters))


def main():
    print("Нелинейные уравнения\n")

    print(f"{f}\n")

    print(f"Количество корней на промежутке [-10, 10]: {n(sturm_system, -10) - n(sturm_system, 10)}")

    bounds = get_bounds(f, -10, 10)
    print("Корни, входящие в промежуток:")
    print(bounds)
    print('\n')

    test_method(binary_search)
    test_method(chord_method)
    test_method(newton)
    print('\n')
    print('Результат, используя встроенные функции')
    print(f.r)


if __name__ == '__main__':
    main()
