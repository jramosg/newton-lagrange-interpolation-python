"""Newton interpolation with divided differences.

This module implements Newton's divided-difference interpolation formula for
arbitrary one-dimensional nodes. It can be imported as a small utility module or
executed directly to run a plotting example.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def divided_difference_table(x_values, y_values):
    """Return the divided-difference table for the given interpolation nodes."""
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x_values and y_values must be one-dimensional.")
    if len(x) != len(y):
        raise ValueError("x_values and y_values must have the same length.")
    if len(x) == 0:
        raise ValueError("At least one interpolation point is required.")
    if len(np.unique(x)) != len(x):
        raise ValueError("x_values must not contain repeated nodes.")

    n = len(x)
    table = np.zeros((n, n), dtype=float)
    table[:, 0] = y

    for order in range(1, n):
        for row in range(n - order):
            numerator = table[row + 1, order - 1] - table[row, order - 1]
            denominator = x[row + order] - x[row]
            table[row, order] = numerator / denominator

    return table


def newton_coefficients(x_values, y_values):
    """Return the Newton polynomial coefficients from divided differences."""
    table = divided_difference_table(x_values, y_values)
    return table[0, :]


def evaluate_newton_polynomial(coefficients, x_nodes, x_eval):
    """Evaluate a Newton-form interpolation polynomial.

    The polynomial is evaluated with nested multiplication:

    P(x) = c0 + (x - x0)(c1 + (x - x1)(c2 + ...)).
    """
    coefficients = np.asarray(coefficients, dtype=float)
    x_nodes = np.asarray(x_nodes, dtype=float)
    x_eval = np.asarray(x_eval, dtype=float)

    if len(coefficients) != len(x_nodes):
        raise ValueError("coefficients and x_nodes must have the same length.")

    result = np.full_like(x_eval, coefficients[-1], dtype=float)

    for i in range(len(coefficients) - 2, -1, -1):
        result = result * (x_eval - x_nodes[i]) + coefficients[i]

    return result


def interpolate_newton(x_values, y_values, x_eval):
    """Build and evaluate the Newton interpolating polynomial."""
    coefficients = newton_coefficients(x_values, y_values)
    return evaluate_newton_polynomial(coefficients, x_values, x_eval)


def main():
    """Run an example interpolation and plot the resulting polynomial."""
    x_nodes = np.array([0, 2, 4, 6, 8, 10], dtype=float)
    y_nodes = np.array([0, 1.5, 2, 2.5, 3, 3], dtype=float)

    point = 2.25
    value = interpolate_newton(x_nodes, y_nodes, point)
    print(f"P({point}) = {float(value):.6f}")

    x_plot = np.linspace(x_nodes[0], x_nodes[-1], 400)
    y_plot = interpolate_newton(x_nodes, y_nodes, x_plot)

    plt.figure(figsize=(8, 5))
    plt.scatter(x_nodes, y_nodes, marker="*", c="purple", s=180, label="Data")
    plt.plot(x_plot, y_plot, "b-", label="Newton interpolation")
    plt.scatter([point], [value], c="red", zorder=3, label=f"P({point})")
    plt.xlabel("x")
    plt.ylabel("P(x)")
    plt.title("Newton Interpolation Method")
    plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.6)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
