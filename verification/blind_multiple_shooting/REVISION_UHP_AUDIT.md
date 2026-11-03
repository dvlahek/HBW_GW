# CQG revision upper-half-plane audit

This directory contains the physical-sheet upper-half-plane audit added for the CQG revision of **“Cold collisionless matter can retain hidden dynamics.”**

The tested window is

```text
-0.05 <= Re(M omega) <= 1.0
 0.001 <= Im(M omega) <= 0.5
```

so the contour crosses `Re(M omega) = 0` and includes purely imaginary growing roots within the tested range.

The stored publication-resolution `q = 0.10` spectral input is used. Bulk roots are followed with residual-gated predictor continuation and adaptive subdivision. The acceptance gate is a full-bulk relative singular-value residual below `1e-7`.

Results:

```text
384-point contour:
  winding = 0
  min |D| = 0.1530880619696
  max full-bulk relative SV = 8.289e-08
  closure mismatch = 9.367e-13

768-point contour:
  winding = 0
  min |D| = 0.1492964459583
  max full-bulk relative SV = 9.998e-08
  closure mismatch = 1.827e-13
```

No growing `ell = 2` root is detected in this finite window. This is a finite-window linear spectral audit and is not a nonlinear stability theorem.

The reproducible implementation is `upper_half_plane_audit.py`. The GitHub Actions workflow is `.github/workflows/upper-half-plane-audit.yml`. Frozen numerical summaries are stored in `upper_half_plane_audit_summary.txt` and `upper_half_plane_audit_summary.json`.
