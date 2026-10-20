# Лабораторная работа 1.4. Метод вращений (Якоби) для симметрической матрицы.
# Запуск: python 1_4_jacobi.py input_1_4.txt

import sys
import math


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

    eps = float(data[pos])

    return n, A, eps


def show(name, M):
    print(name)
    for row in M:
        print("  " + " ".join("%12.6f" % x for x in row))
    print()


def mul(A, B):
    n = len(A)
    C = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0.0
            for k in range(n):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


def diff(A, B):
    d = 0.0
    for i in range(len(A)):
        for j in range(len(A)):
            if abs(A[i][j] - B[i][j]) > d:
                d = abs(A[i][j] - B[i][j])
    return d


def max_offdiag(A):
    n = len(A)
    best = 0.0
    p = 0
    q = 1
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j]) > best:
                best = abs(A[i][j])
                p = i
                q = j
    return best, p, q


def jacobi(A, eps):
    n = len(A)
    A = [row[:] for row in A]

    V = [[0.0] * n for _ in range(n)]
    for i in range(n):
        V[i][i] = 1.0

    while True:
        off, p, q = max_offdiag(A)

        if off <= eps:
            break

        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]

        # угол поворота: phi = 0.5*arctg(2*apq/(aqq-app))
        phi = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c = math.cos(phi)
        s = math.sin(phi)

        # A' = R^T*A*R, меняются только строки и столбцы p и q
        for k in range(n):
            if k != p and k != q:
                akp = A[k][p]
                akq = A[k][q]
                A[k][p] = c * akp - s * akq
                A[p][k] = A[k][p]
                A[k][q] = s * akp + c * akq
                A[q][k] = A[k][q]

        A[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        A[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        A[p][q] = 0.0
        A[q][p] = 0.0

        # накопление собственных векторов: V = V*R
        for k in range(n):
            vkp = V[k][p]
            vkq = V[k][q]
            V[k][p] = c * vkp - s * vkq
            V[k][q] = s * vkp + c * vkq

    values = [A[i][i] for i in range(n)]
    return values, V


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 1_4_jacobi.py input_1_4.txt")
        return

    n, A, eps = read_file(sys.argv[1])

    values, V = jacobi(A, eps)

    Lambda = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Lambda[i][i] = values[i]

    print("1. Заданная точность")
    print("   eps = %.3e" % eps)
    print()

    print("2. Найденные собственные значения")
    for i in range(n):
        print("   lambda[%d] = %.10f" % (i, values[i]))
    print()

    print("3. Найденные собственные векторы и матрица собственных векторов")
    for j in range(n):
        print("   v[%d] = (%s)" % (j, ", ".join("%.8f" % V[i][j] for i in range(n))))
    print()
    show("   Матрица собственных векторов V:", V)

    print("4. Проверка A*V = V*Lambda")
    show("   A*V:", mul(A, V))
    show("   V*Lambda:", mul(V, Lambda))
    print("   max|A*V - V*Lambda| = %.3e" % diff(mul(A, V), mul(V, Lambda)))


main()
