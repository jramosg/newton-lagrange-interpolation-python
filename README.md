Learn Newton interpolation and Lagrange interpolation in Python with
step-by-step examples.

# Metodo de interpolacion de Newton en Python

> Newton interpolation, divided differences, Chebyshev nodes, Runge phenomenon,
> Lagrange interpolation, barycentric interpolation and cubic splines in Python.

Repositorio educativo para estudiantes de programacion cientifica, analisis
numerico y computacion cientifica. El proyecto implementa el metodo de
interpolacion de Newton con diferencias divididas y reproduce un experimento
comparativo con interpolacion baricentrica, Lagrange, Newton y splines cubicos.

## Contenido

- Implementacion reutilizable del polinomio interpolador de Newton.
- Tabla de diferencias divididas.
- Nodos equiespaciados y nodos de Chebyshev en un intervalo general.
- Comparacion de errores para:
  - `f1(x)=sin(x)`
  - `f2(x)=1/(1+25x**2)` funcion de Runge
  - `f3(x)=exp(-20x**2)`
- Error maximo absoluto y norma L2 aproximada.
- Medicion separada de tiempo de construccion y tiempo de evaluacion.
- Graficas opcionales y exportacion de resultados a CSV.

## Instalacion

```shell
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

SciPy se usa para los metodos baricentrico, Lagrange y spline cubico. Si no
esta instalado, el script sigue funcionando con el metodo de Newton.

## Uso rapido

Ejecutar el experimento completo del notebook:

```shell
python polynomial_interpolation.py
```

Guardar resultados en CSV:

```shell
python polynomial_interpolation.py --csv results/interpolation_results.csv
```

Guardar graficas como imagenes PNG:

```shell
python polynomial_interpolation.py --save-plots figures
```

Cambiar el numero de nodos o repeticiones para tiempos:

```shell
python polynomial_interpolation.py --nodes 11 21 41 --repetitions 50
```

## Usar Newton desde otro archivo

```python
import numpy as np
from polynomial_interpolation import (
    evaluate_newton_polynomial,
    newton_coefficients,
)

x_nodes = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
y_nodes = 1 / (1 + 25 * x_nodes**2)

coefficients = newton_coefficients(x_nodes, y_nodes)
value = evaluate_newton_polynomial(coefficients, x_nodes, 0.25)

print(float(value))
```

## Fundamento matematico

Dado un conjunto de puntos:

```text
(x0, y0), (x1, y1), ..., (xn, yn)
```

el polinomio de Newton se escribe como:

```text
P(x) = c0
     + c1(x - x0)
     + c2(x - x0)(x - x1)
     + ...
     + cn(x - x0)(x - x1)...(x - x(n-1))
```

Los coeficientes `c0, c1, ..., cn` se obtienen de la primera fila de la tabla de
diferencias divididas.

## Nodos de Chebyshev

El experimento usa las raices de Chebyshev:

```text
x_k = cos(((2k - 1) * pi) / (2n)), k = 1, ..., n
```

y despues las transforma al intervalo `[a, b]`. Estos nodos ayudan a reducir
las oscilaciones en los extremos del intervalo, especialmente para la funcion de
Runge.

## Resultados esperados

En general, con los mismos nodos, Newton, Lagrange y el metodo baricentrico
construyen el mismo polinomio interpolador, por lo que sus errores son muy
parecidos. La diferencia principal esta en estabilidad numerica y coste
computacional.

Los nodos de Chebyshev suelen reducir mucho el error frente a nodos
equiespaciados en interpolacion polinomica global. Los splines cubicos son
robustos porque no usan un unico polinomio global de grado alto, sino polinomios
de bajo grado por tramos.

## English

# Newton Interpolation Method in Python

Educational Python repository for numerical analysis, scientific computing and
engineering students. It implements Newton interpolation with divided
differences and includes a complete comparison against barycentric
interpolation, SciPy Lagrange interpolation and cubic spline interpolation.

## Features

- Newton divided-difference table.
- Newton polynomial coefficients and nested evaluation.
- Equispaced nodes and Chebyshev nodes.
- Benchmark functions:
  - `sin(x)`
  - Runge function: `1/(1+25x**2)`
  - Gaussian-like function: `exp(-20x**2)`
- Maximum absolute error and approximate L2 error.
- Build time, evaluation time and total runtime.
- Optional plots and CSV export.

## Quick Start

```shell
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python polynomial_interpolation.py
```

Save the benchmark table:

```shell
python polynomial_interpolation.py --csv results/interpolation_results.csv
```

Save plots:

```shell
python polynomial_interpolation.py --save-plots figures
```

## Why This Repository Is Useful

This project is designed as a clear student reference for:

- Newton interpolation in Python
- divided differences in Python
- Chebyshev nodes in Python
- Runge phenomenon examples
- polynomial interpolation comparison
- cubic spline interpolation with SciPy
- numerical methods coursework

## File Structure

```text
index.html                  # SEO static tutorial page for VPS deployment
polynomial_interpolation.py # implementation, CLI experiment and plotting
requirements.txt           # Python dependencies
robots.txt                 # crawler access for static deployment
README.md                  # Spanish and English documentation
```

## References

- NumPy documentation: <https://numpy.org/doc/stable/>
- SciPy interpolate documentation:
  <https://docs.scipy.org/doc/scipy/reference/interpolate.html>
- Newton polynomial interpolation:
  <https://pythonnumericalmethods.berkeley.edu/notebooks/chapter17.05-Newtons-Polynomial-Interpolation.html>
