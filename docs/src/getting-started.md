# Getting started

This page assumes you have never used Julia. It covers the toolchain, a first run, and
the output format; the model those numbers describe is on the [model page](model.md),
and the Julia mechanics of the code itself are in the [code guide](code-guide.md).

## Step 1: install Julia

Download Julia **version 1.12** from [julialang.org](https://julialang.org/downloads/)
and install it as you would any other application. Confirm it worked:

```sh
julia --version
```

If the terminal reports that the command is not found, the installer did not add Julia
to your path; the download page explains how for each system.

## Step 2: choose how you want the code

There are two ways in, and they are not interchangeable — pick by what you intend to do.

### (a) Install the package, to run your own simulations

`Pkg.add` fetches the package straight from GitHub. Nothing is cloned, and you can
`using LaserPlasma` from anywhere:

```sh
julia -e 'using Pkg; Pkg.add(url = "https://github.com/LaurentPlagne/LaserPlasma.jl")'
```

or, equivalently, from the Julia prompt — press `]` to enter package mode:

```julia
pkg> add https://github.com/LaurentPlagne/LaserPlasma.jl
```

To pin a particular commit or tag rather than tracking the default branch, add
`#v0.1.0` or `#<sha>` to the URL.

This gives you everything under `src/`: the plasma types, the forces, both
integrators, the diagnostics and the dielectric theory. It does **not** give you the
`examples/`, `figures/` or `docs/` directories, because those are not part of the
installed package — they live only in the repository.

A complete first script, which you can save anywhere and run with
`julia first_run.jl`:

```julia
using LaserPlasma
const PC = LaserPlasma

p = PC.PlasmaParams(24, 1, 1.442, 6.935, 0.144, 0.144, 1.57, 3.0, 10, 2, 2, true)
cfg = PC.RunConfig(12.0, 8, joinpath(pwd(), "first_run"))
PC.run_multistep(p, cfg, 0.01; sample = 2, thermostat = true)
```

### (b) Clone the repository, to reproduce the figures and examples

The archived 2000 data, the figure scripts, the worked examples and the source of this
guide exist only in the repository:

```sh
git clone https://github.com/LaurentPlagne/LaserPlasma.jl.git
cd LaserPlasma.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` reads `Project.toml` and installs the dependencies at known-good
versions. It takes a minute or two the first time and is instant afterwards.

Two pieces of that command recur throughout this guide:

- `--project=.` means *use the environment described by the directory I am in*. Without
  it, Julia uses your personal global environment, which does not contain this package
  — the cause of essentially every `Package LaserPlasma not found`.
- `-e '…'` evaluates a fragment of Julia and exits, rather than opening the prompt.

**The rest of this page follows route (b)**, from the repository root, because the
examples and figures are the quickest way to see the code do something.

## Step 3: run your first simulation

```sh
julia --project=. examples/first_run.jl
```

This takes **a few seconds** and prints one line: the folder where it wrote its
results, `runs/first_run/`.

`N = 24`, `xfinal = 12`, direct forces (no tree), adaptive individual steps,
thermostat on, `F₀ = 1.57`, `ω = 3`, `n_relax = n_increase = 2` — so the field is off
until `t ≈ 4.2`, ramping until `t ≈ 8.4`, at full amplitude after that.

It is a toy at `N = 24`, but it is the same code path as the 2,000-electron runs and
small enough to experiment with freely.

## Step 4: look at what came out

`runs/first_run/` holds six files: tab-separated text, one row per recorded instant,
atomic units throughout. The full column specification is on the
[model page](model.md); in brief:

| File | Columns |
|---|---|
| `energy.dat` | `t`, kinetic, potential, total |
| `power.dat` | 12 columns; `P = j_elec·E` is column 11 |
| `deltaU` | `t`, `ΔU` over the cycle, running total; first row is the reference (`ΔU = 0`) |
| `pos_0001`, `mome_0001` | final electron positions and momenta, one particle per line |
| `ion_pos.dat` | ion positions |

```sh
head -3 runs/first_run/energy.dat
cat runs/first_run/deltaU
```

`E_tot` jumps between the first and second rows — the thermostat setting the
temperature on the first step — then fluctuates by a few units per step.

`deltaU` has **two** lines and its single measured `ΔU` will very likely be negative.
Full amplitude is only reached at `t ≈ 8.4`, so `xfinal = 12` contains one complete
cycle, and at `N = 24` the per-cycle fluctuation swamps the drift. Nothing is wrong;
it is the expected consequence of the statistics noted on the [model page](model.md).

Positions are wrapped into the box, so a coordinate can jump by `L` between recorded
instants — the same effect that produces the long straight segments in the trajectory
figures of the [gallery](gallery.md).

## Step 5: change something

`examples/first_run.jl` is about twenty lines, one commented parameter per line. The
`PlasmaParams(` block lists the physical settings positionally, in the order given by
its docstring (`?LaserPlasma.PlasmaParams`).

Two lines are worth decoding because their arguments are easy to misread:

- `RunConfig(12.0, 8, outdir)` — final time, then **outer steps per half laser
  period** (the coarsest resolution before adaptive refinement), then the output
  directory.
- `run_multistep(p, cfg, 0.01; …)` — the `0.01` is the integrator's **error
  tolerance**, not a softening length. See the [glossary](glossary.md), where the
  three confusable small numbers are set side by side.

Raising `12.0` to `120.0` turns one measured cycle into fifty-three at no noticeable
cost: the few seconds you waited were Julia compiling, not the simulation running.
Change `F0` too if you want a second configuration to compare — but **change the
output directory as well**, since a run overwrites its folder.

`n_pulse` sets nothing here; run length comes from `RunConfig.xfinal`. The field is
retained for fidelity to the C input format.

## Step 6: rebuild the historical figures

```sh
julia --project=. figures/build_all.jl
```

This regenerates every plot and diagram in this guide: PNG and PDF versions into
`figures/out/`, and copies for these web pages into `docs/src/assets/figures/`.

It is fast, because it mostly redraws the *archived* numbers from 2000 that ship
with the repository, and recomputes the theoretical curves in Julia. It does **not**
launch the expensive 2,000-electron simulations. Two of the scripts — the ones
evaluating dielectric theory — take noticeably longer than the rest.

To produce a single figure instead of all of them, run just its script; the
[gallery](gallery.md) names the script for each figure.

## Step 7 (optional): rebuild this guide

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html` in a web browser. The images are already in the
repository, so this works even if you skipped step 6.

## Checking that the installation is sound

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

This runs the automated test suite: it exercises the random number generation, the
Ewald tables, the forces, the tree, both integrators, and the diagnostics, and
compares each against reference values frozen from the original C program. If every
line reports a pass, your copy is producing the same numbers as the 2000 code.

## If something goes wrong

| Symptom | Likely cause and fix |
|---|---|
| `julia: command not found` | Julia is installed but not on your path; see the installer's instructions. |
| `Package LaserPlasma not found` | Following route (b): you omitted `--project=.`, or the terminal is not in the repository root. Following route (a): the `Pkg.add` did not complete — rerun it. |
| `could not open file examples/first_run.jl` | `examples/` ships only with the cloned repository, not with the installed package. Either clone (route b) or use the standalone script in step 2(a). |
| The first command is very slow | Expected: Julia compiles each method on first call and caches the result. Subsequent runs are fast. |
| A figure script complains about a missing data file | Run it from the repository root, not from inside `figures/`. Each script finds its data relative to its own location, but expects the project environment of the root. |
| Results changed and you do not know why | Check whether you edited a parameter and reused the same output folder. Each run overwrites its folder. |

## Where to go next

[The model page](model.md) gives the exact specification of what was just computed —
parameters, conventions, output columns. [Algorithms](methods.md) covers how it is
computed and how the port is verified, and the [code guide](code-guide.md) covers the
Julia itself.
