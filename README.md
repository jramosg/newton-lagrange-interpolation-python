# Newton Interpolation Method in Python

This project implements Newton's divided-difference interpolation method in
Python. The implementation works with unequally spaced interpolation nodes and
can evaluate the resulting polynomial at one point or over a full array of
points.

## Method

Given interpolation points

```text
(x0, y0), (x1, y1), ..., (xn, yn)
```

Newton interpolation writes the polynomial as

```text
P(x) = c0
     + c1(x - x0)
     + c2(x - x0)(x - x1)
     + ...
     + cn(x - x0)(x - x1)...(x - x(n-1))
```

where the coefficients `c0, c1, ..., cn` are obtained from the divided-difference
table.

## Files

- `newton_interpolation.py`: implementation and executable plotting example.

## Usage

Install the required libraries:

```shell
pip install numpy matplotlib
```

Run the example:

```shell
python newton_interpolation.py
```

Use the functions from another script:

```python
import numpy as np
from newton_interpolation import (
    newton_coefficients,
    evaluate_newton_polynomial,
)

x_nodes = np.array([0, 2, 4, 6], dtype=float)
y_nodes = np.array([1, 3, 2, 5], dtype=float)

coefficients = newton_coefficients(x_nodes, y_nodes)
value = evaluate_newton_polynomial(coefficients, x_nodes, 2.5)

print(value)
```