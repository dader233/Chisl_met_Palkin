# -*- coding: utf-8 -*-
# Лабораторная работа 1.5. QR-разложение и QR-алгоритм поиска собственных значений.
# Запуск: python 1_5_qr.py input_1_5.txt

import sys
import math
import cmath


def read_file(path):
    with open(path, encoding="utf-8") as f:
        data = f.read().split()

    n = int(data[0])
    pos = 1

    A = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(complex(float(data[pos]), 0.0))
            pos += 1
        A.append(row)

    eps = float(data[pos])

    return n, A, eps


def fmt(z):
    re = z.real
    im = z.imag
    if abs(re) < 1e-10:
        re = 0.0
    if abs(im) < 1e-10:
        im = 0.0
    if im == 0.0:
        return "%14.8f" % re
    if re == 0.0:
        return "%14.8fi" % im
    return "%12.8f%+.8fi" % (re, im)


def show(name, M):
    print(name)
    for row in M:
        print("  " + " ".join(fmt(x) for x in row))
    print()


def mul(A, B):
    n = len(A)
    m = len(B[0])
    C = [[0j] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0j
            for k in range(len(B)):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


def qr(A):
    # QR-разложение методом Хаусхолдера (отражений);
    # комплексная арифметика нужна для комплексных сдвигов
    n = len(A)
    R = [row[:] for row in A]
    Q = [[0j] * n for _ in range(n)]
    for i in range(n):
        Q[i][i] = 1.0 + 0j

    for k in range(n - 1):
        # норма части столбца ниже диагонали
        norm = 0.0
        for i in range(k, n):
            norm += abs(R[i][k]) ** 2
        norm = math.sqrt(norm)
        if norm == 0.0:
            continue

        # вектор отражения v = x - alpha*e1
        x1 = R[k][k]
        if abs(x1) == 0.0:
            phase = 1.0 + 0j
        else:
            phase = x1 / abs(x1)
        alpha = -phase * norm

        v = [0j] * (n - k)
        v[0] = x1 - alpha
        for i in range(k + 1, n):
            v[i - k] = R[i][k]

        vv = 0.0
        for z in v:
            vv += abs(z) ** 2
        if vv == 0.0:
            continue

        # R = H*R
        for j in range(k, n):
            s = 0j
            for i in range(k, n):
                s += v[i - k].conjugate() * R[i][j]
            s *= 2.0 / vv
            for i in range(k, n):
                R[i][j] -= s * v[i - k]

        # Q = Q*H
        for i in range(n):
            s = 0j
            for j in range(k, n):
                s += Q[i][j] * v[j - k]
            s *= 2.0 / vv
            for j in range(k, n):
                Q[i][j] -= s * v[j - k].conjugate()

        # ниже диагонали должны быть точные нули
        for i in range(k + 1, n):
            R[i][k] = 0j

    return Q, R


def shift(B, m):
    # сдвиг Уилкинсона: берём ближайшее к B[m-1][m-1] собственное значение блока 2x2
    a = B[m - 2][m - 2]
    b = B[m - 2][m - 1]
    c = B[m - 1][m - 2]
    d = B[m - 1][m - 1]

    root = cmath.sqrt(((a - d) / 2.0) ** 2 + b * c)
    mu1 = (a + d) / 2.0 + root
    mu2 = (a + d) / 2.0 - root

    if abs(mu1 - d) < abs(mu2 - d):
        return mu1
    return mu2


def qr_algorithm(A, eps, max_it=10000):
    # QR-алгоритм со сдвигом Уилкинсона и дефляцией
    B = [row[:] for row in A]
    m = len(B)
    values = []
    it = 0

    while m > 0:
        if m == 1:
            values.append(B[0][0])
            break

        # максимальный по модулю элемент последней строки левее диагонали
        s1 = 0.0
        for j in range(m - 1):
            if abs(B[m - 1][j]) > s1:
                s1 = abs(B[m - 1][j])

        if s1 <= eps:
            # последняя строка обнулилась — последний элемент собственный
            values.append(B[m - 1][m - 1])
            m -= 1
            continue

        # проверка блока 2x2 в правом нижнем углу
        s2 = 0.0
        for j in range(m - 2):
            if abs(B[m - 2][j]) > s2:
                s2 = abs(B[m - 2][j])
            if abs(B[m - 1][j]) > s2:
                s2 = abs(B[m - 1][j])

        if s2 <= eps:
            # собственные значения блока 2x2
            a = B[m - 2][m - 2]
            b = B[m - 2][m - 1]
            c = B[m - 1][m - 2]
            d = B[m - 1][m - 1]

            root = cmath.sqrt(((a - d) / 2.0) ** 2 + b * c)
            values.append((a + d) / 2.0 + root)
            values.append((a + d) / 2.0 - root)
            m -= 2
            continue

        if it >= max_it:
            print("QR-алгоритм не сошёлся за %d итераций" % max_it)
            break

        mu = shift(B, m)

        C = [[0j] * m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                C[i][j] = B[i][j]
            C[i][i] -= mu

        Q, R = qr(C)

        D = mul(R, Q)
        for i in range(m):
            for j in range(m):
                B[i][j] = D[i][j]
            B[i][i] += mu

        it += 1

    return values


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 1_5_qr.py input_1_5.txt")
        return

    n, A, eps = read_file(sys.argv[1])

    print("1. Исходная матрица A")
    show("   A =", A)

    Q, R = qr(A)
    print("2. Матрицы Q и R")
    show("   Q =", Q)
    show("   R =", R)

    print("3. Результат умножения Q*R")
    show("   Q*R =", mul(Q, R))

    values = qr_algorithm(A, eps)
    print("4. Найденные собственные значения")
    for i in range(len(values)):
        z = values[i]
        # мнимая часть в пределах заданной точности считается нулевой
        if abs(z.imag) <= eps:
            z = z.real
        print("   lambda[%d] = %s" % (i, fmt(z)))


main()
