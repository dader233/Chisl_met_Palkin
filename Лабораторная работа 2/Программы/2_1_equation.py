# -*- coding: utf-8 -*-
# Лабораторная работа 2.1. Метод простой итерации и метод Ньютона
# для нелинейного уравнения x^4 - 2x - 1 = 0.
# Запуск: .venv\Scripts\python.exe 2_1_equation.py input_2_1.txt

import os
import sys

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
    x0 = float(data[1])

    return eps, x0


def f(x):
    return x ** 4 - 2.0 * x - 1.0


def df(x):
    return 4.0 * x ** 3 - 2.0


def phi(x):
    # уравнение x^4 - 2x - 1 = 0 приведено к виду x = (2x + 1)^(1/4)
    return (2.0 * x + 1.0) ** 0.25


def simple_iteration(x0, eps):
    x = x0
    history = [x]
    k = 0

    while True:
        x_new = phi(x)
        k += 1
        history.append(x_new)

        if abs(x_new - x) <= eps:
            return x_new, k, history

        x = x_new


def newton(x0, eps):
    x = x0
    history = [x]
    k = 0

    while True:
        x_new = x - f(x) / df(x)
        k += 1
        history.append(x_new)

        if abs(x_new - x) <= eps:
            return x_new, k, history

        x = x_new


def draw_graph(x0, xroot):
    xmin, xmax = 1.0, 1.7
    xs = []
    ys = []

    for i in range(401):
        x = xmin + (xmax - xmin) * i / 400.0
        xs.append(x)
        ys.append(f(x))

    plt.figure(figsize=(8, 5))
    plt.plot(xs, ys, "b-", linewidth=2, label="f(x) = x^4 - 2x - 1")
    plt.axhline(0.0, color="black", linewidth=1)
    plt.axvline(x0, color="green", linestyle="--", label="x0 = %.1f" % x0)
    plt.plot([xroot], [0.0], "ro", label="x* = %.4f" % xroot)

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("График f(x) = x^4 - 2x - 1")
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    if len(sys.argv) < 2:
        print("Укажите файл с данными: python 2_1_equation.py input_2_1.txt")
        return

    eps, x0 = read_file(sys.argv[1])

    # эталонное значение корня нужно только для оценки погрешности
    x_ref, _, _ = newton(x0, 1e-14)

    x_s, k_s, hist_s = simple_iteration(x0, eps)
    x_n, k_n, hist_n = newton(x0, eps)

    print("1. Уравнение и начальное приближение")
    print("   x^4 - 2x - 1 = 0")
    print("   Положительный корень отделён на отрезке [1.3, 1.4]:")
    print("   f(1.3) = %.6f, f(1.4) = %.6f" % (f(1.3), f(1.4)))
    print("   Начальное приближение (по графику): x0 = %.6f" % x0)
    print("   Точность: eps = %.3e" % eps)
    print()

    print("2. Метод простой итерации")
    print("   x = (2x + 1)^(1/4)")
    print("   Количество итераций: %d" % k_s)
    print("   Корень: x = %.12f" % x_s)
    print("   Погрешность |x - x*| = %.3e" % abs(x_s - x_ref))
    print()

    print("3. Метод Ньютона")
    print("   x(k+1) = x(k) - f(x(k))/f'(x(k))")
    print("   Количество итераций: %d" % k_n)
    print("   Корень: x = %.12f" % x_n)
    print("   Погрешность |x - x*| = %.3e" % abs(x_n - x_ref))
    print()

    print("4. Зависимость погрешности от количества итераций")
    print("   (погрешность относительно x* = %.12f)" % x_ref)
    print("      k     простые итерации        Ньютон")

    kmax = max(k_s, k_n)
    for k in range(kmax + 1):
        if k <= k_s:
            s = "%.6e" % abs(hist_s[k] - x_ref)
        else:
            s = ""

        if k <= k_n:
            n = "%.6e" % abs(hist_n[k] - x_ref)
        else:
            n = ""

        print("   %4d   %-18s %s" % (k, s, n))

    print()
    print("   Простые итерации сходятся линейно, Ньютон — квадратично.")

    draw_graph(x0, x_ref)


main()
