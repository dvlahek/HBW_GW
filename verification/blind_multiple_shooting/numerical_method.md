# Numerical method

## Vacuum sector

Both Schwarzschild vacuum regions are solved on the real radial axis in the
Riccati form of the Zerilli equation using `scipy.integrate.solve_ivp` with the
DOP853 integrator.

The horizon solution starts from a Frobenius expansion at `R = 2 + eps` and is
integrated outward to `R = 3`. The outgoing solution starts from an eighth-order
asymptotic expansion at a finite `Rmax` and is integrated inward to `R = 3`.
No exterior complex scaling or contour rotation is used.

The vacuum implementation was checked against the standard Schwarzschild
`l = 2` fundamental Zerilli mode. The production verification used `Rmax = 40`,
where the relative complex error is `1.19e-5`. The benchmark values for other
outer boundaries are provided in `schwarzschild_benchmark.csv`.

## Einstein--Vlasov bulk

The even-parity linearized Einstein tensor is derived symbolically from the
specified metric ansatz in `bulk_blind.py`. The collisionless stress response
is then added and projected onto seven independent equations, giving a `7 x 4`
homogeneous matrix.

The frequency-dependent coefficients `chi_J(omega)` and `chi_X(omega)` are
evaluated at every trial frequency from the publication-resolution spectral
nodes. The calculation uses the retarded piecewise-linear Cauchy continuation
implemented in `spectral.py`. No fixed constitutive coefficients are used in
the root search.

For every trial frequency, two physical radial exponents are obtained from a
four-equation minor and then checked against the full row-normalized `7 x 4`
matrix by singular-value decomposition.

## Matter propagation and global condition

The two bulk eigenvectors are mapped to the Zerilli state `(Z, Z_X)` at the
smooth vacuum edge. Since the matter equations have constant coefficients in
the ocean coordinate, each bulk mode is propagated exactly by
`exp(lambda_i L)`, with `L = 9 ln(1/q)`.

The inner ingoing state is propagated through the matter region and required
to be parallel to the independently computed outer outgoing state. The
resulting normalized determinant is the root-search objective.

## Blind search

The search window is

```text
0.34 < Re(M omega) < 0.42
-0.05 < Im(M omega) < -0.005
```

Twelve frequency seeds are combined with two bulk-root seed pairs, giving 24
attempts. All attempts, including failures, are stored in
`all_frequency_seeds.csv`. Converged roots are deduplicated at a tolerance of
`1e-4` in complex frequency.

The single stable candidate is also checked at three outer-vacuum radii and by
a 24-point argument-principle contour. The winding number is one.

## Independence of the numerical check

The blind calculation was performed within the PiPi institutional research
project at the Faculty of Organization and Informatics. The target matter-mode
frequency, manuscript, and production solver were withheld until the search
and audits were complete. The vacuum formulation differs from both the
production exterior-complex-scaling calculation and the Chebyshev
boundary-value cross-check.

Anthropic Claude was used as a coding assistant during implementation. The
researchers defined the physical specification, blind protocol, search window,
and validation criteria and inspected the resulting code and numerical audits
before unblinding.
