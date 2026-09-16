# Force schemes: how the electron forces are evaluated, as a type rather than an
# integer flag (the C `do_tree`). `DirectScheme` is the pair sum of the oracle,
# `TreeScheme` the Barnes-Hut path with the frozen interaction lists of
# `docs/architecture/barnes-hut.md`. Both integrators dispatch on this.

"""
    ForceScheme

How `multistep_step!` evaluates the electron forces: [`DirectScheme`](@ref)
(pair sums, the default) or [`TreeScheme`](@ref) (Barnes-Hut, C `do_tree = 1`).
"""
abstract type ForceScheme end

"""
    DirectScheme()

Pair-sum forces: e-e + laser + ions per particle, as `build_force_i` with
`do_tree = 0`.
"""
struct DirectScheme <: ForceScheme end

"""
    PerParticleScheme()

Direct pair sums assembled **per particle**, as `build_force_i!` does for
`multistep_rec` — electron `i` sums over every other electron and every ion on its
own. [`DirectScheme`](@ref) computes the same physics but through `acceleration!`,
which exploits Newton's third law (`f[j] -= f[i]`, half the pairs) and is the
summation the C `verlet` froze; that is why it stays the default and is left alone.

Two reasons to want this one for a fixed-step run:

  * it parallelizes, because nothing is written twice — the symmetric loop cannot,
    and at N=500 a `verlet_step!` is 7.8 ms of which 2.7 ms is that serial e-e sum;
  * it makes a fixed-step run and an adaptive run sum **in the same order**, so a
    paired comparison of the two integrators isolates the integrator instead of
    mixing it with the force assembly.

It costs twice the arithmetic of the symmetric loop and is *not* bit-identical to the
oracle's `verlet`: use it for measurements in a new regime, never for reproducing an
archived run.
"""
struct PerParticleScheme <: ForceScheme end

"""
    TreeScheme(tree, arena, starts, lens, mass, charge, eps_tree, L)
    TreeScheme(p::PlasmaParams; theta = 0.5, usequad = 2)

Barnes-Hut force scheme: the shared [`Octree`](@ref), the interaction lists
packed in one flat `arena` buffer (`starts[i]`/`lens[i]` delimit electron `i`,
`lens[i] = 0` for a stable particle), the body masses and charges (electrons
then ions) and the softening length. The lists are rebuilt at every recursion
level with the stable mask; the ions are bodies of the tree, as in the C.
"""
mutable struct TreeScheme{Q} <: ForceScheme
    tree::Octree{Q}
    arenas::Vector{Vector{Int}}
    athread::Vector{Int}
    starts::Vector{Int}
    lens::Vector{Int}
    mass::Vector{Float64}
    charge::Vector{Float64}
    eps_tree::Float64
    L::Float64
end

# Below this population the threading overhead is not worth it. Used for the tree's
# list rebuild, where the per-particle work is O(log N). The force loops use
# `par_worth` instead (`forces.jl`), which weighs the work rather than the population.
const PAR_THRESHOLD = 512


function TreeScheme(p::PlasmaParams; theta::Float64 = 0.5, usequad::Int = 2)
    (usequad == 0 || usequad == 1 || usequad == 2) ||
        error("TreeScheme: usequad must be 0, 1 or 2")
    ne = nb_elec(p)
    if usequad == 2
        return _treescheme(Val(2), p, theta, ne)
    elseif usequad == 1
        return _treescheme(Val(1), p, theta, ne)
    end
    return _treescheme(Val(0), p, theta, ne)
end

function _treescheme(::Val{Q}, p::PlasmaParams, theta::Float64, ne::Int) where {Q}
    L = box_length(p)
    tree = _octree(Val(Q), ne + p.nb_ion, L, theta, true)
    mass = vcat(ones(ne), fill(Float64(p.q_ion), p.nb_ion))
    charge = vcat(fill(Q_ELEC, ne), fill(Float64(p.q_ion), p.nb_ion))
    arenas = [Int[] for _ in 1:Threads.nthreads()]
    return TreeScheme{Q}(tree, arenas, ones(Int, ne), zeros(Int, ne), zeros(Int, ne),
                         mass, charge, p.eps_tree, L)
end

# Preparation hooks mirroring the C tree branches around the shared integrator.
# The direct scheme has nothing to prepare.

"""
    prepare_step!(scheme, deb, ions)

Called once per outer step before the reference RK4 (C `multistep_rec` tree
block): build the tree and the full interaction lists.
"""
prepare_step!(::Union{DirectScheme,PerParticleScheme}, deb::ParticleSet,
              ions::Matrix{Float64}) = nothing

"""
    prepare_stage!(scheme, deb, ions)

Called once per RK4 stage before the forces (C `euler1_rec` tree block): refresh
the body positions and the cell moments, keeping the topology.
"""
prepare_stage!(::Union{DirectScheme,PerParticleScheme}, deb::ParticleSet,
               ions::Matrix{Float64}) = nothing

"""
    prepare_level!(scheme, stable)

Called at the entry of every `step_rec!` level (C tree block): rebuild the tree
and the interaction lists with the stable mask.
"""
prepare_level!(::Union{DirectScheme,PerParticleScheme},
               stable::AbstractVector{Bool}) = nothing

function prepare_step!(s::TreeScheme, deb::ParticleSet, ions::Matrix{Float64})
    set_bodies!(s.tree, s.mass, s.charge, deb.pos, ions)
    maketree!(s.tree)
    fill_lists!(s, nothing)
    return s
end

function prepare_stage!(s::TreeScheme, deb::ParticleSet, ions::Matrix{Float64})
    set_bodies!(s.tree, s.mass, s.charge, deb.pos, ions)
    refresh_moments!(s.tree)
    return s
end

function prepare_level!(s::TreeScheme, stable::AbstractVector{Bool})
    maketree!(s.tree)
    fill_lists!(s, stable)
    return s
end

# Pack the lists of the unstable electrons into the flat arenas with the current
# topology (C `make_interaction_list_stable`); the stable ranges are empty and
# never read. One arena per thread: a rebuild is parallel over electron chunks
# and allocation-free once the high-water marks are reached.
function fill_lists!(s::TreeScheme, stable)
    for a in s.arenas
        empty!(a)
    end
    ne = length(s.lens)
    nt = length(s.arenas)
    if nt == 1 || ne < PAR_THRESHOLD
        _fill_lists_chunk!(s, stable, 1, ne, 1)
    else
        chunk = cld(ne, nt)
        Threads.@threads :static for c in 1:nt
            lo = (c - 1) * chunk + 1
            hi = min(c * chunk, ne)
            lo <= hi && _fill_lists_chunk!(s, stable, lo, hi, c)
        end
    end
    return s
end

function _fill_lists_chunk!(s::TreeScheme, stable, lo::Int, hi::Int, a::Int)
    arena = s.arenas[a]
    n = 0
    @inbounds for i in lo:hi
        s.athread[i] = a
        s.starts[i] = n + 1
        if stable === nothing || !stable[i]
            n = walk_interaction_list!(arena, n, s.tree, i, s.L)
        end
        s.lens[i] = n - s.starts[i] + 1
    end
    return s
end
"""
    build_force_i!(f, i, t, pos, p, fm, tab, L, ions)

Force on electron `i`, in the C order `build_force_i`: e-e + laser + ions.
"""
function build_force_i!(f::Matrix{Float64}, i::Int, t::Float64,
                        pos::Matrix{Float64}, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, ions::Matrix{Float64})
    n = size(pos, 2)
    xi = pos[1, i]
    yi = pos[2, i]
    zi = pos[3, i]
    felx = 0.0
    fely = 0.0
    felz = 0.0
    @inbounds for j in 1:n
        j == i && continue
        fx, fy, fz = pair_force_ee(fm.ee, tab, L, xi, yi, zi, pos[1, j], pos[2, j],
                                   pos[3, j])
        felx += fx
        fely += fy
        felz += fz
    end
    fpulsex = Q_ELEC * laser_electric_field(p, t)
    fionx = 0.0
    fiony = 0.0
    fionz = 0.0
    @inbounds for ion in eachcol(ions)
        fx, fy, fz = pair_force_eion(fm.ion, tab, L, fm.q_ion * Q_ELEC, xi, yi, zi,
                                     ion[1], ion[2], ion[3])
        fionx += fx
        fiony += fy
        fionz += fz
    end
    f[1, i] = felx + fpulsex + fionx
    f[2, i] = fely + fiony
    f[3, i] = felz + fionz
    return f
end

"""
    build_force_i!(f, i, t, pos, p, fm, tab, L, ions, scheme::DirectScheme)

Direct scheme: the pair-sum implementation above.
"""
function build_force_i!(f::Matrix{Float64}, i::Int, t::Float64,
                        pos::Matrix{Float64}, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, ions::Matrix{Float64},
                        ::Union{DirectScheme,PerParticleScheme})
    return build_force_i!(f, i, t, pos, p, fm, tab, L, ions)
end

"""
    build_force_i!(f, i, t, pos, p, fm, tab, L, ions, scheme::TreeScheme)

Tree scheme: force from the frozen interaction list of electron `i`
(C `build_force_ele_tree_i`, `-Acc`), laser added; the ions are in the tree, so
no separate ion sum.
"""
function build_force_i!(f::Matrix{Float64}, i::Int, t::Float64,
                        pos::Matrix{Float64}, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, ions::Matrix{Float64},
                        s::TreeScheme)
    s0 = s.starts[i]
    a = s.athread[i]
    acc, _ = tree_force_i(s.tree, i, view(s.arenas[a], s0:s0+s.lens[i]-1), tab, s.L,
                          s.eps_tree)
    fpulsex = Q_ELEC * laser_electric_field(p, t)
    f[1, i] = -acc[1] + fpulsex
    f[2, i] = -acc[2]
    f[3, i] = -acc[3]
    return f
end

