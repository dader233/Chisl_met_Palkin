# -*- coding: utf-8 -*-
# Лабораторная работа 2.2. Метод простой итерации и метод Ньютона
# для системы нелинейных уравнений (вариант 19, a = 1):
#   x1^2 - 2*lg(x2) - 1 = 0
#   x1^2 - a*x1*x2 + a = 0
# Запуск: .venv\Scripts\python.exe 2_2_system.py input_2_2.txt

import os
import sys
import math

# В venv под Windows Tcl иногда не находит свои библиотеки;
# указываем пути базового Python.
if sys.prefix != sys.base_prefix:
    tcl_dir = os.path.join(sys.base_prefix, "tcl", "tcl8.6")
    tk_dir = os.path.join(sys.base_prefix, "tcl", "tk8.6")
    if os.path.isdir(tcl_dir):
        os.environ.setdefault("TCL_LIBRARY", tcl_dir)
    if os.path.isdir(tk_dir):
        os.environ.setdefault("TK_LIBRARY", tk_dir)

import matplotlib.pyplot as plt


def read_file(path):
    with open(path, encoding="utf-8") as f:
        data = f.read().split()

    eps = float(data[0])
    a = float(data[1])
    x1 = float(data[2])
    x2 = float(data[3])

    return eps, a, x1, x2


def f1(x1, x2):
    return x1 * x1 - 2.0 * math.log10(x2) - 1.0


def f2(x1, x2, a):
    return x1 * x1 - a * x1 * x2 + a


def phi1(x2):
    # из первого уравнения: x1 = sqrt(2*lg(x2) + 1)
    return math.sqrt(2.0 * math.log10(x2) + 1.0)


def phi2(x1, a):
    # из второго уравнения: x2 = (x1^2 + a)/(a*x1)
    return (x1 * x1 + a) / (a * x1)


def simple_iteration(x1, x2, a, eps):
    history = [(x1, x2)]
    k = 0

    while True:
        x1_new = phi1(x2)
        x2_new = phi2(x1, a)
        k += 1
        history.append((x1_new, x2_new))

        if max(abs(x1_new - x1), abs(x2_new - x2)) <= eps:
            return x1_new, x2_new, k, history

        x1, x2 = x1_new, x2_new


def newton(x1, x2, a, eps):
    history = [(x1, x2)]
    k = 0

    while True:
        F1 = f1(x1, x2)
        F2 = f2(x1, x2, a)

        # матрица Якоби
        a11 = 2.0 * x1
        a12 = -2.0 / (x2 * math.log(10.0))
        a21 = 2.0 * x1 - a * x2
        a22 = -a * x1

        # решение системы J*dx = -F по формулам Крамера
        det = a11 * a22 - a12 * a21
        d1 = (-F1 * a22 + a12 * F2) / det
        d2 = (-a11 * F2 + F1 * a21) / det

        x1_new = x1 + d1
        x2_new = x2 + d2
        k += 1
        history.append((x1_new, x2_new))

        if max(abs(d1), abs(d2)) <= eps:
            return x1_new, x2_new, k, history

        x1, x2 = x1_new, x2_new


def draw_graph(x10, x20, x1root, x2root, a):
    # первая кривая: x2 = 10^((x1^2 - 1)/2), строим по x2
    xs1 = []
    ys1 = []
    for i in range(401):
        y = 1.6 + (2.8 - 1.6) * i / 400.0
        x = math.sqrt(2.0 * math.log10(y) + 1.0)
        xs1.append(x)
        ys1.append(y)

    # вторая кривая: x2 = (x1^2 + a)/(a*x1)
    xs2 = []
    ys2 = []
    for i in range(401):
        x = 0.9 + (1.6 - 0.9) * i / 400.0
        y = (x * x + a) / (a * x)
        xs2.append(x)
        ys2.append(y)

    plt.figure(figsize=(8, 5.5))
    plt.plot(xs1, ys1, "b-", linewidth=2,
             label="x1^2 - 2lg(x2) - 1 = 0")
    plt.plot(xs2, ys2, "g-", linewidth=2,
             label="x1^2 - a*x1*x2 + a = 0")
    plt.plot([x10], [x20], "o", color="orange",
             label="x0 = (%.1f; %.1f)" % (x10, x20))
    plt.plot([x1root], [x2root], "ro",
             label="x* = (%.4f; %.4f)" % (x1root, x2root))

    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.title("Пересечение кривых (a = %.0f)" % a)
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 2_2_system.py input_2_2.txt")
        return

    eps, a, x10, x20 = read_file(sys.argv[1])

    # эталонное решение нужно только для оценки погрешности
    x1_ref, x2_ref, _, _ = newton(x10, x20, a, 1e-14)

    x1_s, x2_s, k_s, hist_s = simple_iteration(x10, x20, a, eps)
    x1_n, x2_n, k_n, hist_n = newton(x10, x20, a, eps)

    print("1. Система уравнений и начальное приближение")
    print("   x1^2 - 2*lg(x2) - 1 = 0")
    print("   x1^2 - a*x1*x2 + a = 0,  a = %.6f" % a)
    print("   Положительное решение; по графику кривые пересекаются")
    print("   в точке около (1.3; 2.1):")
    print("   начальное приближение: x1 = %.6f, x2 = %.6f" % (x10, x20))
    print("   Точность: eps = %.3e" % eps)
    print()

    print("2. Метод простой итерации")
    print("   x1 = sqrt(2*lg(x2) + 1)")
    print("   x2 = (x1^2 + a)/(a*x1)")
    print("   Количество итераций: %d" % k_s)
    print("   Решение: x1 = %.12f, x2 = %.12f" % (x1_s, x2_s))
    e_s = max(abs(x1_s - x1_ref), abs(x2_s - x2_ref))
    print("   Погрешность max(|x1-x1*|, |x2-x2*|) = %.3e" % e_s)
    print()

    print("3. Метод Ньютона")
    print("   J(x)*dx = -F(x),  x(k+1) = x(k) + dx")
    print("   Количество итераций: %d" % k_n)
    print("   Решение: x1 = %.12f, x2 = %.12f" % (x1_n, x2_n))
    e_n = max(abs(x1_n - x1_ref), abs(x2_n - x2_ref))
    print("   Погрешность max(|x1-x1*|, |x2-x2*|) = %.3e" % e_n)
    print()

    print("4. Зависимость погрешности от количества итераций")
    print("   (погрешность относительно x1* = %.12f, x2* = %.12f)" % (x1_ref, x2_ref))
    print("      k     простые итерации        Ньютон")

    kmax = max(k_s, k_n)
    for k in range(kmax + 1):
        if k <= k_s:
            s = "%.6e" % max(abs(hist_s[k][0] - x1_ref),
                             abs(hist_s[k][1] - x2_ref))
        else:
            s = ""

        if k <= k_n:
            n = "%.6e" % max(abs(hist_n[k][0] - x1_ref),
                             abs(hist_n[k][1] - x2_ref))
        else:
            n = ""

        print("   %4d   %-18s %s" % (k, s, n))

    print()
    print("   Простые итерации сходятся линейно, Ньютон — квадратично.")

    draw_graph(x10, x20, x1_ref, x2_ref, a)


main()
