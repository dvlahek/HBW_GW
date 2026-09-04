# External and analytic benchmarks

These checks are independent of the new Einstein--Vlasov matter-led pole itself.

1. **Published Schwarzschild QNM.** The manuscript Zerilli/ECS vacuum solver gives
   `M omega = 0.3736716843524924 - 0.0889623156416965 i`.
   The independently published Schwarzschild `l=2,n=0` value in
   B. Carneiro da Cunha and J. P. Cavalcante, JHEP 08 (2024) 110 is
   `0.3736716844180418 - 0.0889623156889357 i`.
   Relative complex difference: `2.10e-10`.

2. **Closed-form retarded continuation.** For
   `rho(nu)=1+nu` on `0<=nu<=1`,
   the exact Cauchy transform is

       C(z) = (z+1)[log(z)-log(z-1)] - 1.

   On the retarded second sheet crossed through `0<Re z<1`,

       C_II(z) = C(z) - 2 pi i (1+z).

   The same piecewise-linear continuation routine used in the publication code
   reproduces the exact continuation with maximum absolute error `2.31e-14`
   over the supplied test points.

3. **Published Einstein-cluster endpoint.** The Einstein-cluster equations in
   K. Jusufi, Eur. Phys. J. C 83, 103 (2023) give

       m/r = 2 omega/(1+4 omega),
       rho = omega/[2 pi r^2(1+4 omega)],
       P_t = omega rho,
       f proportional to r^(4 omega).

   Setting `omega=1/2` gives exactly the cold endpoint used in the manuscript:

       m/r = 1/3,
       rho = 1/(12 pi r^2),
       P_r = 0,
       P_t = rho/2,
       f proportional to r^2.

Run:

    python reproduce_external_analytic_benchmarks.py
