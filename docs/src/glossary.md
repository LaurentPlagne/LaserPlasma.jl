# Glossary

Computational vocabulary, the parameters of this code that are easy to confuse, and a
short list of conventions inherited from the C original. Physics terms are not defined
here; the [model page](model.md) states the model, and terms specific to *this
implementation* appear below.

## Parameters that are easy to confuse

**`r_pot` / `eps_tree` / `eps`.** Three small positive numbers, three different roles.
`r_pot` is the electron–ion soft-core length and `eps_tree` the electron–electron one —
both physics. The `eps` positional argument of `run_multistep` is the **integration
error tolerance** of the step-doubling test — arithmetic accuracy, not a model
parameter. Mistaking the third for the first two yields a plausible answer to a
different question.

**`nb_step` / `level_max` / `n_pulse`.** `nb_step` is the number of *outer* steps per
half laser period, i.e. the coarsest resolution before any adaptive refinement.
`level_max` caps how deep the adaptive recursion may then go (historical runs reach
12). `n_pulse` sets nothing in this port — run length comes from `RunConfig.xfinal` —
and is retained only for fidelity to the C input format.

**`theta` / `usequad`.** Both control the tree approximation but independently. `theta`
is the opening criterion: a cell is accepted when its size-to-distance ratio falls
below it, so *smaller `theta` is more accurate and slower*. `usequad` is how many terms
of the multipole expansion each accepted cell contributes: `0` monopole, `1` +dipole,
`2` +quadrupole. Historical runs: `theta = 0.5`, `usequad = 2`.

**`sample`.** Selects which realization of the random initial conditions to use, by
consuming the preceding draws exactly as the C skip loop did. Different `sample` values
give statistically independent runs — this is the mechanism for obtaining an error bar,
not a seed in the usual sense.

## Computational terms

**Molecular dynamics (MD).** Integrate every particle's equations of motion
individually, with no fluid or kinetic closure. What this code does.

**Structure of arrays (SoA).** Storing a collection of records as one array per field
(`pos`, `charge`, `quad`, …) rather than one array of record objects. A loop reading a
few fields of many records then reads contiguous memory instead of dragging whole
records through cache. Used for the octree and, in effect, for the `3 × N` particle
matrices. Contrast *array of structures* (AoS), the `Vector{Node}` this code avoids.

**Column-major.** Julia stores matrices column by column, so in a `3 × N` particle
array a single particle's `x, y, z` are adjacent in memory. The layout is chosen for
this reason.

**Step doubling.** Estimating the error of a step by taking it twice: once whole, once
as two halves, and comparing. The basis of the adaptive integrator here, with the
comparison made per particle rather than on a global norm.

**Freezing / interpolation.** In the adaptive scheme, a particle whose step-doubling
test passes is *frozen* for the remainder of the outer interval and is not integrated
again at deeper levels; where its position is needed there, it is interpolated
(quadratic Lagrange through the three RK4 nodes already computed). This is what makes
the scheme cheap rather than merely selective.

**Interaction list.** The set of tree nodes — bodies and accepted cells — that
contribute to the force on one particle, produced by the tree walk. Here they are
packed into flat per-thread arenas with `starts[i]`/`lens[i]` delimiting each particle,
rather than as a vector of vectors. A frozen particle has `lens[i] = 0`.

**Frozen topology.** Within one recursion level, the tree's *structure* is not rebuilt,
though the cell *moments* are refreshed at every RK4 stage. Rebuilding the structure
between stages would make the right-hand side discontinuous within a single step.

**Workspace / preallocation.** Allocating every buffer a hot routine needs once, up
front, and reusing it. `MultistepWork` is this package's; allocating inside the
recursion instead costs a measured factor of 190.

**Type parameter, compile-time specialization.** Putting a choice into a *type*
(`Octree{Q}`, `Forces{E,I}`) so the compiler generates a separate specialized version
of the code for each case and the corresponding branches vanish from the inner loops.
See [Code guide](code-guide.md).

**Multiple dispatch.** Selecting which method of a function runs from the types of all
its arguments. Used here to replace the C's integer flags: choosing an algorithm means
choosing a type.

**Type stability.** A function is type-stable when the compiler can determine its
return type from its argument types alone. Instability forces boxed values and runtime
checks and is a performance bug even when results are correct; the test suite asserts
it with `@inferred`.

**Oracle.** The original C program, which still compiles and runs, and against which
this port's numbers are checked. Verification here means numerical agreement with a
running reference, not a plausibility argument.

**Bit-identical.** Agreeing to the last bit of the floating-point representation, not
merely to within tolerance. Achievable when the same arithmetic is performed in the same
order — which is why the RNG is reproduced exactly and why threading the direct force
loop is legal (it changes no summation order) while `PerParticleScheme` is not
bit-identical (it changes the pairing).

## Environment and tooling

**Package, `Project.toml`, `Manifest.toml`.** The unit of reusable Julia code, the file
listing its direct dependencies and version bounds, and the file recording the exact
resolved versions of everything transitively. `--project=.` selects the environment
described by the current directory.

**REPL.** Julia's interactive prompt. `?name` prints a docstring; this package's
docstrings carry design rationale, not just signatures.

**Compilation latency.** Julia compiles each method the first time it is called with a
given set of argument types. A script's first seconds are usually compilation, not
computation, and it is cached between runs.

**Documenter.** The Julia package that renders `docs/src/*.md` into the HTML you are
reading.

**Archived data.** The small files under `figures/data/`, produced by the C program in
2000 and bundled with the repository. The result plots in the [gallery](gallery.md) are
drawn from them, so redrawing one does not rerun a simulation.
