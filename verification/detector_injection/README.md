# Idealized detector-level ringdown injection/recovery

This directory contains a synthetic feasibility test for the predicted `q = 0.1` matter-led ringdown component.

The injected signal contains the Schwarzschild `l = 2` vacuum fundamental and the predicted matter-led mode

```text
M omega_kin = 0.377514814460 - 0.020871654486 i
```

in Gaussian noise shaped with an analytic Advanced-LIGO design-sensitivity PSD proxy. The test compares a vacuum-only damped-sinusoid model with a vacuum + fixed-frequency matter-mode model using

```text
DeltaBIC = BIC(vacuum) - BIC(vacuum + kinetic).
```

Strong evidence is defined as `DeltaBIC > 10`.

The scan uses masses `30, 60, 100, 180 Msun`, full-window ringdown SNRs `12, 20, 30, 40`, injected amplitude ratios `Akin/Avac = 0, 0.02, 0.05, 0.10, 0.20, 0.30`, and analysis starts `0M, 10M, 20M`. Each grid point uses 100 Monte Carlo injections with fixed seed `1701`.

Main result: no zero-kinetic control produced `DeltaBIC > 10`. At full-window SNR 40, an injected ratio `Akin/Avac = 0.20` reaches at least 90% strong-evidence recovery for all tested masses and start times. At SNR 30, the corresponding threshold is `0.20` to `0.30` depending on mass and start time. The tested ratios do not reach the 90% criterion at SNR 12 or 20.

This is an idealized detector-level separability test, not an analysis of measured gravitational-wave data and not an astrophysical detectability forecast. The matter-mode frequency and damping time are fixed to the theoretical prediction. Astrophysical excitation amplitudes are not modeled.

Run
---

```bash
python hbh_detector_injection_recovery.py
```

The script regenerates trial-level, summary, and threshold CSV files. Compact publication-facing outputs are also stored in `results/detector_injection_summary.csv` and `results/detector_thresholds.csv`.
