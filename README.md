# LaserPlasma.jl

[![CI](https://github.com/LaurentPlagne/LaserPlasma.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/LaurentPlagne/LaserPlasma.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/LaurentPlagne/LaserPlasma.jl/actions/workflows/Documentation.yml/badge.svg)](https://laurentplagne.github.io/LaserPlasma.jl/)
[![codecov](https://codecov.io/gh/LaurentPlagne/LaserPlasma.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LaurentPlagne/LaserPlasma.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

**How much energy does a laser leave behind in a hot plasma, and why?**

This is a Julia port of a molecular-dynamics simulation written in C at RWTH Aachen
in 1999–2000. It follows individual electrons moving among fixed ions in a periodic
box driven by a laser field, and measures the energy they absorb — the process known
as collisional absorption, or inverse bremsstrahlung, which governs how a laser heats
a fusion target.

The port keeps the original's numbers (it is validated against the C program digit by
digit) while organizing the code for a reader. It includes Ewald summation for the
periodic forces, an optional Barnes–Hut tree, relativistic dynamics, and two time
integrators — a plain fixed-step one and the adaptive individual-step scheme that was
the original work's contribution.

## Documentation

**[Read the guide](https://laurentplagne.github.io/LaserPlasma.jl/).** The physics is stated as a
[specification](https://laurentplagne.github.io/LaserPlasma.jl/model/) rather than motivated; the explanatory weight is on
the computational side — [the algorithms](https://laurentplagne.github.io/LaserPlasma.jl/methods/) (tabulated Ewald
correction, Barnes–Hut, the recursive individual-step integrator and its coupling to
the tree, and the performance engineering under all of it) and
[the code](https://laurentplagne.github.io/LaserPlasma.jl/code-guide/), which assumes no Julia.

The [figure gallery](https://laurentplagne.github.io/LaserPlasma.jl/gallery/) reproduces **every result and method figure
cited by the original 2000 presentation**, recording for each whether the curve is
recomputed by this port or redrawn from archived output. The 2000 results ship as small
data files, so redrawing a plot does not rerun a 2,000-electron simulation.

## Installation

Requires Julia 1.12. To use the package:

```julia
using Pkg; Pkg.add(url = "https://github.com/LaurentPlagne/LaserPlasma.jl")
```

That installs `src/` — the plasma types, forces, both integrators, diagnostics and
dielectric theory. The `examples/`, `figures/` and `docs/` directories, and the archived
2000 data, are not part of the installed package; for those, clone the repository:

```sh
git clone https://github.com/LaurentPlagne/LaserPlasma.jl.git
cd LaserPlasma.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'

julia --project=. examples/first_run.jl        # a 24-electron run, a few seconds
julia --project=. figures/build_all.jl         # rebuild every figure
julia --project=. -e 'using Pkg; Pkg.test()'   # check against the frozen reference values
```

The example writes to `runs/first_run/`; the figure command writes to `figures/out/` and
refreshes the images used by the guide. To build the guide locally:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl     # then open docs/build/index.html
```
