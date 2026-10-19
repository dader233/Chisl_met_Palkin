# -*- coding: utf-8 -*-
# Лабораторная работа 1.3. Метод простых итераций и метод Зейделя.
# Запуск: python 1_3_iter.py input_1_3.txt

import sys


def read_file(path):
    with open(path, encoding="utf-8") as f:
        data = f.read().split()

    n = int(data[0])
    pos = 1

    A = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(float(data[pos]))
            pos += 1
        A.append(row)

    b = []
    for i in range(n):
        b.append(float(data[pos]))
        pos += 1

    eps = float(data[pos])

    return n, A, b, eps


def gauss(A, b):
    # эталонное решение, чтобы было с чем сравнивать итерации
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for k in range(n):
        p = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > abs(M[p][k]):
                p = i
        M[k], M[p] = M[p], M[k]

        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]
    return x


def simple_iteration(A, b, eps):
    # x = B*x + c, где B = E - D^(-1)*A, c = D^(-1)*b
    n = len(A)
    B = [[0.0] * n for _ in range(n)]
    c = [0.0] * n

    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = -A[i][j] / A[i][i]
        c[i] = b[i] / A[i][i]

    x = c[:]
    it = 0

    while True:
        x_new = [0.0] * n
        for i in range(n):
            s = c[i]
            for j in range(n):
                s += B[i][j] * x[j]
            x_new[i] = s

        d = 0.0
        for i in range(n):
            if abs(x_new[i] - x[i]) > d:
                d = abs(x_new[i] - x[i])

        x = x_new
        it += 1

        if d <= eps:
            break

    return x, it


def seidel(A, b, eps):
    n = len(A)
    # то же начальное приближение, что и в методе простых итераций
    x = [b[i] / A[i][i] for i in range(n)]
    it = 0

    while True:
        d = 0.0
        for i in range(n):
            s = b[i]
            for j in range(n):
                if j != i:
                    s -= A[i][j] * x[j]
            new = s / A[i][i]

            if abs(new - x[i]) > d:
                d = abs(new - x[i])
            x[i] = new

        it += 1
        if d <= eps:
            break

    return x, it


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 1_3_iter.py input_1_3.txt")
        return

    n, A, b, eps = read_file(sys.argv[1])

    x1, it1 = simple_iteration(A, b, eps)
    x2, it2 = seidel(A, b, eps)
    x_ref = gauss(A, b)

    e1 = max(abs(x1[i] - x_ref[i]) for i in range(n))
    e2 = max(abs(x2[i] - x_ref[i]) for i in range(n))

    print("1. Заданная точность")
    print("   eps = %.3e" % eps)
    print()

    print("2. Количество итераций для каждого метода")
    print("   Метод простых итераций: %d итераций" % it1)
    print("   Решение: x = (" + ", ".join("%.10f" % v for v in x1) + ")")
    print("   Метод Зейделя: %d итераций" % it2)
    print("   Решение: x = (" + ", ".join("%.10f" % v for v in x2) + ")")
    print()

    print("3. Эталонное решение для проверки результата")
    print("   x = (" + ", ".join("%.10f" % v for v in x_ref) + ")")
    print()

    print("4. Сравнение методов по скорости сходимости")
    print("   Простые итерации: %d итераций, ошибка %.3e" % (it1, e1))
    print("   Зейдель:          %d итераций, ошибка %.3e" % (it2, e2))
    if it2 < it1:
        print("   Зейдель сошёлся быстрее в %.2f раза" % (it1 / it2))
    elif it1 < it2:
        print("   Простые итерации сошлись быстрее в %.2f раза" % (it2 / it1))
    else:
        print("   одинаковое число итераций")


main()
