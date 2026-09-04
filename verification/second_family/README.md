# Second self-consistent Einstein--Vlasov family

This directory documents the independent profile-family robustness check used in the manuscript.

The alternative stationary family is

```text
F_alt(E,L) = C (E0-E)_+^2 (L-L0)_+^2.
```

The high-resolution continuation path was selected using stationary proper-shell bracketing before any spectral or QNM result was evaluated. The cold-orbit audit passes through N=80. The final global pole is

```text
M omega_alt = 0.375287456824 - 0.021057162748 i
```

which is 0.591145% from the reference q=0.10 pole in complex frequency. The root has unit winding, max full 7x4 relative singular-value residual 1.36e-12, and all ten seed combinations collapse to one deduplicated root.

The files in `results/` contain the publication audit tables. This is a profile-robustness result, not a universality theorem over arbitrary collisionless distributions.
