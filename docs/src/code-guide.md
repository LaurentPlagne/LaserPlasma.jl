# Code guide

This page is Julia for someone who is fluent in the physics and the mathematics but
does not write this kind of code. It covers the language mechanics this package
actually uses, the design decisions behind it, and where to find things.

## Why Julia, concretely

The original is C. The usual reason to leave C is that one ends up writing two
languages — a fast kernel and a slow driver around it — and the usual reason not to is
that the alternatives are slower. Julia is compiled ahead of each first call, to native
code, specialized on the actual argument types; the result here is a port that is
*faster* than the C at `-O2` while being readable, which is the only justification that
matters.

The cost is that the first call to anything is slow, because it is being compiled. This
is why `julia --project=. examples/first_run.jl` spends a few seconds starting before it
computes anything, and why the simulation itself is almost free by comparison. It is
compilation, not slowness, and it is cached between runs of the same code.

## The project environment

Julia's dependency management is per-project, and this matters more than it sounds.

- `Project.toml` lists the direct dependencies and their allowed versions.
- `Manifest.toml` records the exact resolved versions of everything, transitively.
- `--project=.` tells Julia to use the environment described by the current directory.

Omit `--project=.` and you get your personal global environment, which does not contain
this package; that is the cause of essentially every `Package LaserPlasma not found`.
`Pkg.instantiate()` reads those files and installs what they specify. `docs/` has its
own separate environment, which is why building the documentation uses
`--project=docs`.

## Reading the code

**Structure.** `src/LaserPlasma.jl` is the entry point and does nothing but `include`
the other files in dependency order — definitions must precede their use, so the order
is meaningful. Everything lives in a module named `LaserPlasma`; scripts reach it as
`LaserPlasma.run_multistep`, usually via a short local alias.

**Docstrings.** Every significant function carries one, and reading them from the REPL
is faster than searching files:

```julia
julia> using LaserPlasma          # started with: julia --project=.
julia> ?LaserPlasma.run_multistep
```

The docstrings in this package are unusually load-bearing: several record *why* a
routine is written the way it is, and one (`PerParticleScheme`) documents that its
results are deliberately not bit-identical to the reference. Read them before changing
a numerical path.

**`!` is a convention, not syntax.** A trailing `!` warns that the function mutates one
of its arguments — by convention the first. `init_plasma!(state, params)` fills the
`state` you passed in. Julia does not enforce this; the package keeps the promise.

**Arguments after a semicolon are optional and named:**

```julia
run_multistep(p, cfg, 0.01; thermostat = true, use_tree = true, theta = 0.5)
```

The first three are positional and required; everything after the `;` is a keyword
argument with a default, order-independent and omissible. All the run switches work
this way.

**Unicode identifiers are ordinary.** `p.ω` is the laser frequency; `γ`, `ν`, `Γ`, `ε`
appear where the physics uses them. In the REPL they are typed as `\omega` followed by
Tab.

## Types and dispatch: the one idea that shapes this package

This is the part worth understanding, because it explains why the Julia code looks so
unlike the C.

In Julia, a function can have many **methods**, and which one runs is chosen from the
types of *all* the arguments. That mechanism is used here to replace the C's global
integer flags with types:

| C | Here |
|---|---|
| `do_tree` (0/1) | `DirectScheme` / `TreeScheme{Q}` |
| `usequad` (0/1/2) | the parameter `Q` of `Octree{Q}` |
| `stepper` (integer) | `run_verlet` / `run_multistep` |
| implicit potential choice | `SoftCorePot`, carried by `Forces{E,I}` |

Concretely, `prepare_level!` is defined twice: once for `TreeScheme`, where it rebuilds
the tree and the interaction lists, and once for `DirectScheme`, where it returns
`nothing`. The integrator calls `prepare_level!(scheme, stable)` unconditionally and
never asks which scheme it has. There is no `if do_tree` anywhere in the integrator.

The practical consequences:

- **The branch is gone from the inner loop.** The compiler generates one specialized
  version of the integrator per scheme type, and the no-op hook is inlined away to
  literally nothing.
- **Adding a variant means adding methods**, not editing a switch in every routine that
  consults the flag.
- **Nonsense combinations fail to compile** rather than running and producing a number.

**Type parameters** take this one step further. `Octree{Q}` puts the multipole order
into the *type*, so `Octree{2}` and `Octree{0}` are different types and the code that
walks them is compiled separately, with the "how many multipole terms" branches already
resolved. This is why the constructor contains what looks like a redundant
`if usequad == 2 ... Val(2)`: it converts a runtime integer into a compile-time type
parameter, once, at construction. Similarly, `Forces{E,I}` carries the electron–electron
and electron–ion potential types, so the potential evaluation in the pair loop is a
direct call rather than a dispatch through a pointer.

If you take one thing from this page: **in this code, choosing an algorithm means
choosing a type.**

## Memory layout

**Particles are `3 × N` matrices.** Rows are `x, y, z`; columns are particles;
`pos[2, 17]` is the `y` of particle 17. Julia is column-major, so a particle's three
coordinates are adjacent in memory — which is what the pair loops want, since they read
all three together.

**The octree is a structure of arrays.** Not a `Vector{Node}` of linked objects but
parallel arrays `pos`, `charge`, `dip`, `quad`, `next`, `more`, with a node's identity
being its column index and `0` meaning "none". A tree walk reads one or two fields of
very many nodes; with node objects each visit would drag a whole node through cache to
read a fraction of it. [Algorithms](methods.md) covers the reasoning.

**Hot-path buffers are preallocated.** `MultistepWork` holds every array the recursion
needs, allocated once per run. The tree's node arrays are recycled across rebuilds.
Allocating inside the recursion costs a measured factor of 190, so if you add code
there, add its buffer to the workspace rather than creating arrays locally.

You will see `@inbounds` (skip bounds checking) and `@inline` in these loops. They are
safe only where the indices are already known to be in range; do not scatter them.

## Where things are

| To understand… | Start at |
|---|---|
| Parameters and particle state | `src/plasma.jl`: `PlasmaParams`, `PlasmaState`, `init_plasma!` |
| Reproducible random draws | `src/rng.jl`, `src/sobol.jl` |
| Laser envelope and field | `src/laser.jl` |
| Periodic correction and its table | `src/ewald.jl` |
| Potentials and pair forces | `src/force_types.jl`, `src/forces.jl` |
| The octree | `src/tree.jl` |
| Force-method selection and the tree/integrator hooks | `src/schemes.jl` |
| Fixed-step integrator | `src/integrators.jl` |
| Adaptive individual-step integrator | `src/multistep.jl` |
| Run loop and checkpointing | `src/simulation.jl`: `RunConfig`, `run_multistep`, `run_verlet` |
| Energy, power, collision frequency | `src/diagnostics.jl`, `src/analysis.jl` |
| Dielectric theory | `src/dielectric.jl` |
| Figures | `figures/*.jl`, shared helpers in `figures/common.jl` |

A workable first pass: `plasma.jl` for the state, `laser.jl` because it is short and
self-contained, `force_types.jl` then `forces.jl` for one interaction, `schemes.jl` for
how the two force paths are selected, then `simulation.jl` for how a run is assembled.
Leave `multistep.jl` and `tree.jl` until after [Algorithms](methods.md).

## Three small numbers that are not the same thing

| Name | Role | Effect of changing it |
|---|---|---|
| `r_pot` | Softening length, electron–ion | Physics |
| `eps_tree` | Softening length, electron–electron | Physics |
| the `eps` argument of `run_multistep` | Integration error tolerance | Arithmetic accuracy only |

All three are small positive floats around `0.01`–`0.15`, and confusing the third with
the first two produces a perfectly plausible answer to a different question.

## Run switches

| Keyword | Effect |
|---|---|
| `use_tree` | Barnes–Hut instead of exact pair sums. Required in practice for large `N`. |
| `theta`, `usequad` | Opening criterion, and multipole order (`0` monopole, `1` +dipole, `2` +quadrupole). Historical runs: `0.5`, `2`. |
| `thermostat` | The constant-temperature protocol. Needed for a clean absorption measurement. |
| `sample` | Selects an independent realization of the initial conditions — this is how you get an error bar. |
| `level_max` | Cap on recursion depth. Historical runs reach 12. |
| `interp` | Interpolation stencil for frozen particles; `:quad` is the C scheme and the default. |
| `checkpoint_every`, `restart` | Periodic state saves and resumption. Essential for long runs. |

## Testing and changing things

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite checks the random generators, Ewald table, forces, tree, both integrators and
the diagnostics against values frozen from the C reference. It also asserts type
stability with `@inferred` on a dozen core entry points — the Ewald table builders, the
tree construction and force walk, `Forces`, `PlasmaState`, `energies`. A `@inferred`
failure means the compiler could not determine a concrete return type, usually because
a variable can hold more than one, and that is a performance bug even when the numbers
come out right.

What the suite does and does not establish: it pins down agreement with the reference on
short, controlled cases. It says nothing about whether a long run's `ν_ei` is
trustworthy — that needs several `sample` values and an honest look at the scatter, for
the reasons in [the model page](model.md).

If you change a numerical path and a reference value moves, the question to answer is
*which* of the two is wrong. The C is reconstructible and can be instrumented and
recompiled; speculation about it is not a substitute.

## Beyond this guide

Detailed validation records and research notes on the integrator are kept as internal
working documents and are not part of this repository. The measured agreements they
establish are summarized at the end of [Algorithms](methods.md).
