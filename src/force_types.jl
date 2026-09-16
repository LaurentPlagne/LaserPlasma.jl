# Force model: potentials are types, with the softening parameter carried by
# each instance. The dynamics currently use only the soft-core potential;
# `figures/cutoff.jl` compares it with other potentials for illustration.

"""
    PairPotential

Abstract pair potential: a behaviour, selected by dispatch.
"""
abstract type PairPotential end

"""
    SoftCorePot(eps2)

Softened potential `1/√(r²+ε²)`: carries `ε²` (precomputed product, as in the C).
"""
struct SoftCorePot <: PairPotential
    eps2::Float64
end

"""
    softcore(eps)

Build a `SoftCorePot` from the softening length `eps`.
"""
softcore(eps::Real) = SoftCorePot(eps * eps)

"""
    softened_sqdist(pot, x, y, z)

Squared distance `x² + y² + z² + ε²` of the softened potential.
"""
@inline softened_sqdist(pot::SoftCorePot, x, y, z) = x * x + y * y + z * z + pot.eps2

"""
    Forces(ee, ion, q_ion)
    Forces(p::PlasmaParams)

Force context: both potentials (e-e and e-ion) and the ion charge.
"""
struct Forces{E<:PairPotential,I<:PairPotential}
    ee::E
    ion::I
    q_ion::Float64
end

Forces(p::PlasmaParams) = Forces(softcore(p.eps_tree), softcore(p.r_pot), Float64(p.q_ion))
