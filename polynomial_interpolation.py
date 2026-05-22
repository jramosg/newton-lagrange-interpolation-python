"""Polynomial interpolation experiments in Python.

The module contains a reusable implementation of Newton interpolation with
divided differences and a command-line experiment inspired by a scientific
computing activity. It compares:

- Barycentric polynomial interpolation
- Lagrange polynomial interpolation
- Newton polynomial interpolation
- Cubic spline interpolation

The default experiment evaluates sin(x), Runge's function and exp(-20*x**2) on
[-1, 1] with equispaced and Chebyshev nodes.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.interpolate import (
        BarycentricInterpolator,
        InterpolatedUnivariateSpline,
        lagrange,
    )
except ImportError:  # pragma: no cover - exercised only without SciPy.
    BarycentricInterpolator = None
    InterpolatedUnivariateSpline = None
    lagrange = None


Array = np.ndarray
Function = Callable[[Array], Array]


@dataclass(frozen=True)
class TestFunction:
    """Function used by the interpolation benchmark."""

    key: str
    label: str
    function: Function


@dataclass(frozen=True)
class InterpolationResult:
    """Numerical result for one function/node/method combination."""

    function: str
    n_nodes: int
    node_type: str
    method: str
    max_error: float
    l2_error: float
    build_time: float
    eval_time: float

    @property
    def total_time(self) -> float:
        """Return construction time plus evaluation time."""
        return self.build_time + self.eval_time


def f_sin(x: Array) -> Array:
    """Return sin(x)."""
    return np.sin(x)


def f_runge(x: Array) -> Array:
    """Return Runge's classical interpolation test function."""
    return 1 / (1 + 25 * x**2)


def f_gaussian(x: Array) -> Array:
    """Return exp(-20*x**2)."""
    return np.exp(-20 * x**2)


TEST_FUNCTIONS = (
    TestFunction("sin", "f1(x)=sin(x)", f_sin),
    TestFunction("runge", "f2(x)=1/(1+25x**2)", f_runge),
    TestFunction("gaussian", "f3(x)=exp(-20x**2)", f_gaussian),
)


def equispaced_nodes(a: float, b: float, n: int) -> Array:
    """Return n equispaced nodes on [a, b]."""
    validate_interval(a, b)
    validate_node_count(n)
    return np.linspace(a, b, n)


def chebyshev_nodes(a: float, b: float, n: int) -> Array:
    """Return n Chebyshev roots mapped from [-1, 1] to [a, b]."""
    validate_interval(a, b)
    validate_node_count(n)
    k = np.arange(1, n + 1)
    roots = np.cos((2 * k - 1) * np.pi / (2 * n))
    return 0.5 * (b - a) * roots + 0.5 * (a + b)


def divided_difference_table(x_values: Iterable[float],
                             y_values: Iterable[float]) -> Array:
    """Return the Newton divided-difference table."""
    x, y = validated_xy(x_values, y_values)
    n = len(x)
    table = np.zeros((n, n), dtype=float)
    table[:, 0] = y

    for order in range(1, n):
        for row in range(n - order):
            numerator = table[row + 1, order - 1] - table[row, order - 1]
            denominator = x[row + order] - x[row]
            table[row, order] = numerator / denominator

    return table


def newton_coefficients(x_values: Iterable[float],
                        y_values: Iterable[float]) -> Array:
    """Return Newton-form polynomial coefficients."""
    return divided_difference_table(x_values, y_values)[0, :]


def evaluate_newton_polynomial(coefficients: Iterable[float],
                               x_nodes: Iterable[float],
                               x_eval: Iterable[float] | float) -> Array:
    """Evaluate a Newton-form interpolation polynomial."""
    coefficients = np.asarray(coefficients, dtype=float)
    x_nodes = np.asarray(x_nodes, dtype=float)
    x_eval = np.asarray(x_eval, dtype=float)

    if coefficients.ndim != 1 or x_nodes.ndim != 1:
        raise ValueError("coefficients and x_nodes must be one-dimensional.")
    if len(coefficients) != len(x_nodes):
        raise ValueError("coefficients and x_nodes must have the same length.")
    if len(coefficients) == 0:
        raise ValueError("At least one coefficient is required.")

    result = np.full_like(x_eval, coefficients[-1], dtype=float)

    for i in range(len(coefficients) - 2, -1, -1):
        result = result * (x_eval - x_nodes[i]) + coefficients[i]

    return result


def interpolate_newton(x_values: Iterable[float],
                       y_values: Iterable[float],
                       x_eval: Iterable[float] | float) -> Array:
    """Build and evaluate the Newton interpolating polynomial."""
    coefficients = newton_coefficients(x_values, y_values)
    return evaluate_newton_polynomial(coefficients, x_values, x_eval)


def error_metrics(y_true: Array, y_approx: Array,
                  x_grid: Array) -> tuple[float, float]:
    """Return maximum absolute error and approximate L2 error."""
    error = np.abs(y_approx - y_true)
    return float(np.max(error)), float(np.sqrt(np.trapezoid(error**2, x_grid)))


def measure_time(operation: Callable[[], object],
                 repetitions: int = 20) -> tuple[object, float]:
    """Return operation result and median runtime in seconds."""
    if repetitions <= 0:
        raise ValueError("repetitions must be positive.")

    timings = []
    result = None

    for _ in range(repetitions):
        start = perf_counter()
        result = operation()
        timings.append(perf_counter() - start)

    return result, float(np.median(timings))


def run_experiment(n_values: Iterable[int] = (11, 21),
                   node_types: Iterable[str] = ("equispaced", "chebyshev"),
                   x_range: tuple[float, float] = (-1.0, 1.0),
                   grid_size: int = 1000,
                   repetitions: int = 20,
                   make_plots: bool = False,
                   output_dir: Path | None = None) -> list[InterpolationResult]:
    """Run the full interpolation comparison."""
    a, b = x_range
    validate_interval(a, b)
    if grid_size < 2:
        raise ValueError("grid_size must be at least 2.")

    x_grid = np.linspace(a, b, grid_size)
    results: list[InterpolationResult] = []

    for test_function in TEST_FUNCTIONS:
        y_true = test_function.function(x_grid)

        for n in n_values:
            for node_type in node_types:
                x_nodes = make_nodes(node_type, a, b, n)
                y_nodes = test_function.function(x_nodes)
                x_nodes, y_nodes = sort_nodes(x_nodes, y_nodes)

                method_values = evaluate_methods(
                    x_nodes, y_nodes, x_grid, repetitions
                )

                for method, (y_approx, build_time, eval_time) in method_values:
                    max_error, l2_error = error_metrics(
                        y_true, y_approx, x_grid
                    )
                    results.append(
                        InterpolationResult(
                            function=test_function.label,
                            n_nodes=n,
                            node_type=node_type,
                            method=method,
                            max_error=max_error,
                            l2_error=l2_error,
                            build_time=build_time,
                            eval_time=eval_time,
                        )
                    )

                if make_plots:
                    plot_comparison(
                        test_function.label,
                        node_type,
                        x_grid,
                        y_true,
                        x_nodes,
                        y_nodes,
                        method_values,
                        output_dir,
                    )

    return results


def evaluate_methods(x_nodes: Array, y_nodes: Array, x_grid: Array,
                     repetitions: int) -> list[tuple[str, tuple[Array,
                                                                float,
                                                                float]]]:
    """Build and evaluate every available interpolation method."""
    values = []

    if BarycentricInterpolator is not None:
        interpolator, build_time = measure_time(
            lambda: BarycentricInterpolator(x_nodes, y_nodes), repetitions
        )
        y_approx, eval_time = measure_time(
            lambda: interpolator(x_grid), repetitions
        )
        values.append(("barycentric", (y_approx, build_time, eval_time)))

    if lagrange is not None:
        polynomial, build_time = measure_time(
            lambda: lagrange(x_nodes, y_nodes), repetitions
        )
        y_approx, eval_time = measure_time(lambda: polynomial(x_grid),
                                           repetitions)
        values.append(("lagrange", (y_approx, build_time, eval_time)))

    coefficients, build_time = measure_time(
        lambda: newton_coefficients(x_nodes, y_nodes), repetitions
    )
    y_approx, eval_time = measure_time(
        lambda: evaluate_newton_polynomial(coefficients, x_nodes, x_grid),
        repetitions,
    )
    values.append(("newton", (y_approx, build_time, eval_time)))

    if InterpolatedUnivariateSpline is not None:
        spline, build_time = measure_time(
            lambda: InterpolatedUnivariateSpline(x_nodes, y_nodes, k=3),
            repetitions,
        )
        y_approx, eval_time = measure_time(lambda: spline(x_grid), repetitions)
        values.append(("cubic_spline", (y_approx, build_time, eval_time)))

    return values


def make_nodes(node_type: str, a: float, b: float, n: int) -> Array:
    """Create interpolation nodes by name."""
    if node_type == "equispaced":
        return equispaced_nodes(a, b, n)
    if node_type == "chebyshev":
        return chebyshev_nodes(a, b, n)
    raise ValueError(f"Unknown node type: {node_type}")


def sort_nodes(x_nodes: Array, y_nodes: Array) -> tuple[Array, Array]:
    """Return nodes and values sorted by x."""
    idx = np.argsort(x_nodes)
    return x_nodes[idx], y_nodes[idx]


def write_results_csv(results: list[InterpolationResult],
                      output_file: Path) -> None:
    """Write experiment results to a CSV file."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=result_fieldnames())
        writer.writeheader()
        for row in results:
            writer.writerow(result_to_dict(row))


def print_results_table(results: list[InterpolationResult]) -> None:
    """Print a compact table with all experiment results."""
    headers = result_fieldnames()
    rows = [result_to_dict(result) for result in results]
    widths = {
        header: max(len(header), *(len(str(row[header])) for row in rows))
        for header in headers
    }
    print(" | ".join(header.ljust(widths[header]) for header in headers))
    print("-+-".join("-" * widths[header] for header in headers))
    for row in rows:
        print(" | ".join(str(row[header]).ljust(widths[header])
                         for header in headers))


def plot_comparison(function_label: str,
                    node_type: str,
                    x_grid: Array,
                    y_true: Array,
                    x_nodes: Array,
                    y_nodes: Array,
                    method_values: list[tuple[str, tuple[Array,
                                                         float,
                                                         float]]],
                    output_dir: Path | None = None) -> None:
    """Plot one function/node comparison."""
    plt.figure(figsize=(9, 5.5))
    plt.plot(x_grid, y_true, color="black", label="Real function")

    styles = {
        "barycentric": "r--",
        "lagrange": "b-.",
        "newton": "g:",
        "cubic_spline": "m-",
    }
    for method, (y_approx, _, _) in method_values:
        plt.plot(x_grid, y_approx, styles.get(method, "-"), label=method)

    plt.scatter(x_nodes, y_nodes, color="black", s=25, label="Nodes", zorder=3)
    plt.title(f"{function_label} with {len(x_nodes)} {node_type} nodes")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True, alpha=0.35)
    plt.legend()
    plt.tight_layout()

    if output_dir is None:
        plt.show()
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        filename = safe_filename(f"{function_label}_{len(x_nodes)}_{node_type}")
        plt.savefig(output_dir / f"{filename}.png", dpi=160)
        plt.close()


def result_fieldnames() -> list[str]:
    """Return CSV/table field names."""
    return [
        "function",
        "n_nodes",
        "node_type",
        "method",
        "max_error",
        "l2_error",
        "build_time",
        "eval_time",
        "total_time",
    ]


def result_to_dict(result: InterpolationResult) -> dict[str, str | int]:
    """Convert one result to printable/CSV values."""
    return {
        "function": result.function,
        "n_nodes": result.n_nodes,
        "node_type": result.node_type,
        "method": result.method,
        "max_error": f"{result.max_error:.6e}",
        "l2_error": f"{result.l2_error:.6e}",
        "build_time": f"{result.build_time:.6e}",
        "eval_time": f"{result.eval_time:.6e}",
        "total_time": f"{result.total_time:.6e}",
    }


def safe_filename(value: str) -> str:
    """Return a filesystem-friendly filename fragment."""
    return "".join(char if char.isalnum() else "_" for char in value).strip("_")


def validated_xy(x_values: Iterable[float],
                 y_values: Iterable[float]) -> tuple[Array, Array]:
    """Validate interpolation data and return NumPy arrays."""
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

    return x, y


def validate_interval(a: float, b: float) -> None:
    """Validate an interpolation interval."""
    if not a < b:
        raise ValueError("The interval must satisfy a < b.")


def validate_node_count(n: int) -> None:
    """Validate the number of interpolation nodes."""
    if n < 2:
        raise ValueError("At least two interpolation nodes are required.")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare Newton, Lagrange, barycentric and cubic spline "
                    "interpolation on scientific-computing examples."
    )
    parser.add_argument("--nodes", nargs="+", type=int, default=[11, 21],
                        help="Node counts to test. Default: 11 21.")
    parser.add_argument("--grid-size", type=int, default=1000,
                        help="Number of evaluation points. Default: 1000.")
    parser.add_argument("--repetitions", type=int, default=20,
                        help="Timing repetitions. Default: 20.")
    parser.add_argument("--plot", action="store_true",
                        help="Show comparison plots.")
    parser.add_argument("--save-plots", type=Path,
                        help="Directory where plots will be saved as PNG.")
    parser.add_argument("--csv", type=Path,
                        help="Optional CSV output path for the results.")
    return parser.parse_args()


def main() -> None:
    """Run the command-line experiment."""
    args = parse_args()
    output_dir = args.save_plots
    results = run_experiment(
        n_values=args.nodes,
        grid_size=args.grid_size,
        repetitions=args.repetitions,
        make_plots=args.plot or output_dir is not None,
        output_dir=output_dir,
    )
    print_results_table(results)

    if args.csv is not None:
        write_results_csv(results, args.csv)
        print(f"\nResults saved to {args.csv}")

    if BarycentricInterpolator is None:
        print("\nSciPy was not found. Only the Newton method was executed.")


if __name__ == "__main__":
    main()
