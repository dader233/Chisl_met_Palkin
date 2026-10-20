# Лабораторная работа 1.2. Метод прогонки.
# Хранятся только три диагонали и правая часть, полная матрица не строится.
# Запуск: python 1_2_progonka.py input_1_2.txt

import sys


def read_file(path):
    with open(path, encoding="utf-8") as f:
        data = f.read().split()

    n = int(data[0])
    pos = 1

    a = [0.0] * n          # нижняя диагональ, a[0] не используется
    for i in range(1, n):
        a[i] = float(data[pos])
        pos += 1

    b = [0.0] * n          # главная диагональ
    for i in range(n):
        b[i] = float(data[pos])
        pos += 1

    c = [0.0] * n          # верхняя диагональ, c[n-1] не используется
    for i in range(n - 1):
        c[i] = float(data[pos])
        pos += 1

    d = [0.0] * n          # правая часть
    for i in range(n):
        d[i] = float(data[pos])
        pos += 1

    return n, a, b, c, d


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 1_2_progonka.py input_1_2.txt")
        return

    n, a, b, c, d = read_file(sys.argv[1])

    # прямой ход
    alpha = [0.0] * n
    beta = [0.0] * n

    alpha[0] = c[0] / b[0]
    beta[0] = d[0] / b[0]

    for i in range(1, n):
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1])
        beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1])

    
    alpha[n - 1] = 0.0

    print("1. Прогоночные коэффициенты")
    print("      i        alpha[i]            beta[i]")
    for i in range(n):
        print("   %4d %18.10f %18.10f" % (i, alpha[i], beta[i]))
    print()

    # обратный ход
    x = [0.0] * n
    x[n - 1] = beta[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = beta[i] - alpha[i] * x[i + 1]

    print("2. Итоговое решение системы")
    for i in range(n):
        print("   x[%d] = %.10f" % (i, x[i]))


main()
