# Лабораторная работа 1.1. LU-разложение с выбором главного элемента.
# Запуск: python 1_1_lu.py input_1_1.txt

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

    return n, A, b


def show(name, M):
    print(name)
    for row in M:
        print("  " + " ".join("%12.6f" % x for x in row))
    print()


def mul(A, B):
    n = len(A)
    m = len(B[0])
    C = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0.0
            for k in range(len(B)):
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


def lu(A):
    # P*A = L*U. P хранится как список номеров строк.
    n = len(A)
    U = [row[:] for row in A]
    L = [[0.0] * n for _ in range(n)]
    P = list(range(n))
    swaps = 0

    for k in range(n):
        # выбор главного элемента в столбце k
        p = k
        for i in range(k + 1, n):
            if abs(U[i][k]) > abs(U[p][k]):
                p = i

        if abs(U[p][k]) < 1e-15:
            raise ValueError("матрица вырождена")

        if p != k:
            U[k], U[p] = U[p], U[k]
            P[k], P[p] = P[p], P[k]
            for j in range(k):
                L[k][j], L[p][j] = L[p][j], L[k][j]
            swaps += 1

        L[k][k] = 1.0

        for i in range(k + 1, n):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

    return L, U, P, swaps


def solve(L, U, P, b):
    n = len(L)
    pb = [b[P[i]] for i in range(n)]

    # L*y = P*b
    y = [0.0] * n
    for i in range(n):
        s = pb[i]
        for j in range(i):
            s -= L[i][j] * y[j]
        y[i] = s

    # U*x = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
        x[i] = s / U[i][i]

    return x


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 1_1_lu.py input_1_1.txt")
        return

    n, A, b = read_file(sys.argv[1])

    L, U, P, swaps = lu(A)

    Pm = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Pm[i][P[i]] = 1.0

    print("1. Матрицы L и U")
    show("   L =", L)
    show("   U =", U)

    print("2. Результат умножения L*U")
    show("   L*U =", mul(L, U))

    x = solve(L, U, P, b)
    print("3. Решение системы A*x = b")
    print("   x = (" + ", ".join("%.10f" % v for v in x) + ")")
    print()

    Ainv = [[0.0] * n for _ in range(n)]
    for j in range(n):
        e = [0.0] * n
        e[j] = 1.0
        col = solve(L, U, P, e)
        for i in range(n):
            Ainv[i][j] = col[i]

    print("4. Обратная матрица A^(-1), найденная через LU-разложение")
    show("   A^(-1) =", Ainv)

    det = -1.0 if swaps % 2 else 1.0
    for i in range(n):
        det *= U[i][i]
    print("5. Определитель матрицы A, найденный через LU-разложение")
    print("   det(A) = %.10f" % det)
    print()

    E = [[0.0] * n for _ in range(n)]
    for i in range(n):
        E[i][i] = 1.0

    print("6. Проверка A*A^(-1) = E")
    show("   A*A^(-1) =", mul(A, Ainv))
    print("   max|A*A^(-1) - E| = %.3e" % diff(mul(A, Ainv), E))
    print()

    print("7. Проверка L*U = P*A")
    show("   Матрица перестановок P =", Pm)
    show("   P*A =", mul(Pm, A))
    print("   max|L*U - P*A| = %.3e" % diff(mul(L, U), mul(Pm, A)))


main()
