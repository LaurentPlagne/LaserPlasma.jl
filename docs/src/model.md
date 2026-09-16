# The model as implemented

A precise statement of what is simulated, for reference rather than motivation. The
[algorithms](methods.md) that evaluate it are the subject of the next page.

## System

`N` ions of charge `Z = 1`, **fixed** for the whole run, and `N` electrons, in a cubic
box of side

```math
L = \left(\frac{4\pi N}{3}\right)^{1/3}r_{\mathrm{ws}},
```

with periodic boundary conditions in all three directions. Electron dynamics are
relativistic, integrated in the momentum `u = γv`; `d_relat = false` selects the
non-relativistic path.

**Units are atomic** (`mₑ = |e| = 1`, `c = 137`, lengths in Bohr radii). Note that the
runs take `r_ws = 1.442 ≈ 3^{1/3}`, at which density `ω_p = 1` to four digits — which
is why every frequency in the study and in the figure titles is quoted as a bare
multiple of `ω_p`.

## Interactions

Electron–ion and electron–electron pairs both use a soft-core potential
`1/√(r² + ε²)`, with **two independent softening lengths**: `r_pot` for electron–ion
and `eps_tree` for electron–electron. The periodic images are handled by an Ewald
correction (α = 2/L, real and reciprocal sums truncated at 5, tabulated on a 17³ grid
and trilinearly interpolated).

The split between the two is worth stating, because it shapes the code: **the table
holds only `erfc + reciprocal − 1/r`**. The central `1/r` term is excluded from it and
supplied separately by the softened pair potential. Softening therefore never touches
the table, and changing `r_pot` or `eps_tree` does not invalidate it.

Ions are not mutually interacting (they are fixed); in the tree path they are
nonetheless inserted as bodies of the octree, so the electron–ion sum is absorbed into
the tree walk rather than done separately.

## Laser

Uniform in space (dipole approximation), along `x`:

```math
E_x(t) = A(t)\,F_0\cos(\omega t),
```

with `A(t) = 0` for the first `n_relax` laser periods, rising linearly over the next
`n_increase`, and 1 thereafter. The strength parameter used throughout is the quiver
velocity `v₀ = F₀/ω`, reported against the thermal speed as `v₀/v_th`.

## Initial conditions

| Quantity | Draw |
|---|---|
| Ion positions | Sobol sequence, 3 draws per ion |
| Electron positions | Uniform in the box |
| Electron momenta | Gaussian per component, σ = √kT, via Box–Muller (`gasdev`) on `dran2` |

The generators `dran2`, `gasdev` and `sobol_` are reproduced **bit for bit** from the
C original. This is not fidelity for its own sake: it is what makes a particle-by-particle
comparison against the reference program possible at all, and hence what the
verification in [Algorithms](methods.md) rests on. The `sample` keyword selects a
statistically independent realization by consuming the preceding draws exactly as the
C skip loop did.

## Diagnostics

| File | Columns | Content |
|---|---|---|
| `energy.dat` | 4 | `t`, kinetic, potential, total |
| `power.dat` | 12 | `t`, `E_ion` (3), integrated current `j` (3), `j_elec` (3), then `P₁ = j_elec·E` (column 11) and `P₂` |
| `deltaU` | 3 | `t`, `ΔU` over the cycle, running total; first row is the reference, `ΔU = 0` |
| `pos_0001`, `mome_0001`, `ion_pos.dat` | 3 | Final electron positions and momenta; ion positions |

The collision frequency follows from the mean slope of the absorbed energy:

```math
\nu_{ei}=\frac{2}{v_0^2}\frac{d(U_{\mathrm{abs}}/N)}{dt},\qquad v_0=\frac{F_0}{\omega},
```

implemented in `collision_frequency` (`src/analysis.jl`), which also returns the
per-cycle standard deviation and the standard error over the number of cycles.

## Thermostat

Optional, and inherited from the original protocol. During relaxation and ramp-up,
velocities are rescaled each step towards `1.5 N kT`. Afterwards, at each full cycle
(`2·nb_step` steps), `ΔU = E_tot − E_tot_ref` is measured and then removed from the
momenta by the factor `α = (E_kin − ΔU)/E_kin` applied to `u`. Temperature is held
fixed, so cycles are measured under identical conditions and `ΔU` is additive.

The rescaling factor is computed from the RMS velocity spread,
`U_th = 0.5 N (σx² + σy² + σz²)` over `v = u/γ` at the start of the step — despite the
C routine being named `fit_gauss`, nothing is fitted.

## Conventions that bite

- **Output instants differ by integrator.** `run_multistep` records at the **start** of
  an outer step; `run_verlet` at the **end**. This reproduces the C and is deliberate;
  any file comparison must account for the one-step offset.
- **Kinetic energy is only meaningful at cycle boundaries** (`counter−1 ≡ 0 mod
  2·nb_step`). Mid-cycle it is dominated by the quiver: in the `F₀ = 79` run it reaches
  ~7·10⁵ against a cycle-boundary value of ~2·10⁴. Sampling it off-boundary
  manufactures a drift that is not there.
- **`n_pulse` does not set the run length** in this port; `RunConfig.xfinal` does. The
  field is retained for fidelity to the C input format.
- **`xfinal` is truncated to an integer before the step count is derived**, exactly as
  the C did (`maxnbstep = (int) xfinal/delta_x`), so `xfinal < 1` yields zero steps.

## Statistics

Trajectories are chaotic, so only cycle-averaged quantities — `P`, `ΔU`, `ν_ei` — are
comparable between runs or between codes beyond a short initial stretch. The per-cycle
`ΔU` fluctuates considerably more than it drifts, so a comparison between two
configurations without a scatter estimate over several `sample` values establishes
nothing. This is why the demonstration run in [Getting started](getting-started.md),
at `N = 24` with a single measured cycle, can report a negative absorption.
