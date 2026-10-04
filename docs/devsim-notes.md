# DEVSIM implementation notes

The reference model (`devsim/tcad.py`) is a from-scratch 2D drift-diffusion MOSFET on DEVSIM 2.11. These are the numerical problems found while building it and how each was fixed. They're useful when debugging any TCAD deck, Sentaurus included.

| # | Symptom | Cause | Fix |
|---|---|---|---|
| 1 | Contacts matched 0–1 nodes | DEVSIM 2D contacts must lie on a boundary between a simulated region and another region | Surround the device with `gas` regions; contacts on Si/gas and oxide/gas edges |
| 2 | Doping overflowed (`exp` of +10¹⁰) | The parser binds unary minus tighter than `^`: `exp(-(y/s)^2)` means `exp((-y/s)^2)` | Write `exp(0 - (y/s)^2)` |
| 3 | `log(0)` at n+ contacts | `0.5·(−N + √(N²+4nᵢ²))` cancels to 0 in double precision for N = 5×10¹⁹ | Branch-free exact form n = nᵢ·exp(asinh(N/2nᵢ)), p = nᵢ·exp(−asinh(N/2nᵢ)) |
| 4 | Singular Jacobian at triple points | Si/oxide/nitride corner nodes got redundant continuity equations; S/D contact end node also sat on the Si/nitride interface | Drop the oxide/spacer side-wall interfaces; stop S/D contacts 2 nm short of the spacers |
| 5 | Newton diverged at equilibrium | Wrong sign: DEVSIM's Poisson node term is **+q(n − p − N)** | Match the `simple_physics` conventions (documented in `tcad.py`) |
| 6 | Linear (not quadratic) Newton convergence | DEVSIM's symbolic `diff()` does not differentiate `B()` (Bernoulli) inline | Define `Bern01`, `Bern10`, `vdiff` as edge models with explicit derivatives |
| 7 | Stall at ~10⁻⁵ relative error | Cancellation in n₁B(−v) − n₀B(v) at 5×10¹⁹ cm⁻³ | Compensated summation `kahan3(...)` |
| 8 | Source contact solved to −0.58 V instead of +0.58 V | An `ifelse` model first evaluated *inside* contact assembly took the wrong branch and was cached (forcing evaluation beforehand gave the right value) | Removed `ifelse` from contact models (see #3) and force-evaluate helper models |
| 9 | Mixed-mode LU "divide by zero" | MKL PARDISO pivoting failed on matrices mixing O(1) circuit entries with O(10⁻¹⁴) device entries | Switch to UMFPACK (`import devsim.umfpack.umfshim`) |
| 10 | Empty KCL row for `out` at equilibrium | No current flows before drift-diffusion is enabled | 100 TΩ bleed resistor (18 fA, about 500× below Ioff) |
| 11 | Circuit relative error stuck ~10⁻⁷ at 0 V | Relative error of node voltages near zero is noise | rel. tolerance 10⁻⁶ in mixed mode (≈1 µV on a 1.8 V swing) |

A useful debugging habit throughout: **bisect against a known-good reference.** DEVSIM's own `simple_physics` solved on the same mesh, so swapping one piece at a time (parameters, potential equation, contacts, currents) isolated each bug in a few runs.
