# Physical parameters, plasma state and initialization (init_plasma, ewald.c).
#
# The state separates what the C mixed in `y[6(i-1)+k]`: fixed ions on one
# side, positions and momenta `u = γv` of the electrons on the other. No linear
# index is manipulated anywhere.

"""
    PlasmaParams(nb_ion, q_ion, r_ws, kT, r_pot, eps_tree, F0, ω, n_pulse, n_relax,
                 n_increase, d_relat)

Physical parameters of a run, in atomic units.
"""
struct PlasmaParams
    nb_ion::Int
    q_ion::Int
    r_ws::Float64
    kT::Float64
    r_pot::Float64
    eps_tree::Float64
    F0::Float64
    ω::Float64
    n_pulse::Int
    n_relax::Int
    n_increase::Int
    d_relat::Bool
end

"""
    nb_elec(p)

Number of electrons, `Z N_ion`.
"""
nb_elec(p::PlasmaParams) = p.q_ion * p.nb_ion

"""
    box_length(p)

Box side `L = (4π N_ion/3)^{1/3} r_ws`, with the `pow` cube root of the C.
"""
box_length(p::PlasmaParams) = ((4.0 * π / 3.0) * Float64(p.nb_ion))^(1.0 / 3.0) * p.r_ws

"""
    PlasmaState(p)

Named plasma state: electron positions and momenta `u = γv` (3×N), and the fixed
ion positions (3×N_ion). The C flat layout `y[6(i-1)+k]` is not reproduced.
"""
mutable struct PlasmaState
    pos::Matrix{Float64}
    mom::Matrix{Float64}
    ions::Matrix{Float64}
end

PlasmaState(p::PlasmaParams) = PlasmaState(zeros(3, nb_elec(p)), zeros(3, nb_elec(p)),
                                           zeros(3, p.nb_ion))

"""
    init_plasma!(st, p; sample = 1, rng = Ran2(-1), gd = Gasdev(), sobol = Sobol())

Faithful initialization: uniform positions and Maxwellian velocities for the
electrons, ion positions from the Sobol sequence. `sample` reproduces the C skip
loop (the draws of the previous samples are consumed).
"""
function init_plasma!(st::PlasmaState, p::PlasmaParams; sample::Int = 1,
                      rng::Ran2 = Ran2(-1), gd::Gasdev = Gasdev(),
                      sobol::Sobol = Sobol())
    L = box_length(p)
    n = nb_elec(p)
    for _ in 1:sample-1
        for _ in 1:n
            foreach(_ -> dran2!(rng), 1:3)
            foreach(_ -> gasdev!(gd, rng), 1:3)
            # in plasma mode (bckgrd_ion==2) the C consumes three more draws
            foreach(_ -> dran2!(rng), 1:3)
        end
    end
    sigmap = sqrt(p.kT)
    @inbounds for i in 1:n
        for c in 1:3
            st.pos[c, i] = L * (dran2!(rng) - 0.5)
        end
        for c in 1:3
            st.mom[c, i] = sigmap * gasdev!(gd, rng)
        end
    end
    @inbounds for i in 1:p.nb_ion
        x = sobol_next!(sobol)
        for c in 1:3
            st.ions[c, i] = L * (Float64(x[c]) - 0.5)
        end
    end
    return st
end

"""
    ParticleSet
    ParticleSet(pos, mom)

Bundle of electron positions and momenta used as a state at one time.
"""
struct ParticleSet
    pos::Matrix{Float64}
    mom::Matrix{Float64}
end

"""
    particle_set(n)

Allocate a zero `ParticleSet` of `n` electrons.
"""
particle_set(n::Int) = ParticleSet(zeros(3, n), zeros(3, n))

"""
    nparticles(s)

Number of electrons of `s`.
"""
@inline nparticles(s::ParticleSet) = size(s.pos, 2)

"""
    copy_particles!(dst, src)

Copy positions and momenta of all electrons.
"""
function copy_particles!(dst::ParticleSet, src::ParticleSet)
    copyto!(dst.pos, src.pos)
    copyto!(dst.mom, src.mom)
    return dst
end

"""
    copy_stable!(stable, dst, src)

`dst.x[i] = src.x[i]` for stable particles only (C `copy_stable`).
"""
function copy_stable!(stable::AbstractVector{Bool}, dst::ParticleSet, src::ParticleSet)
    @inbounds for i in axes(dst.pos, 2)
        stable[i] || continue
        for c in 1:3
            dst.pos[c, i] = src.pos[c, i]
            dst.mom[c, i] = src.mom[c, i]
        end
    end
    return dst
end

"""
    wrap_in_box!(st, L)

`back_in_cell`: fold the electrons back into the box, with the C bounds
(`> L/2` and `<= -L/2`), applied after every outer step.
"""
function wrap_in_box!(st::PlasmaState, L::Float64)
    half = 0.5 * L
    mhalf = -0.5 * L
    @inbounds for r in eachcol(st.pos), c in 1:3
        while r[c] > half
            r[c] -= L
        end
        while r[c] <= mhalf
            r[c] += L
        end
    end
    return st
end
