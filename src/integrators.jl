# Integrators. For now the relativistic leapfrog from `verlet` (ppbs.c), the
# fixed-step reference scheme: positions in x, velocities stored as momenta
# `u = γv`.

"""
    verlet_forces!(f, st, p, fm, tab, L, t, scheme)

Total force on every electron at time `t`, under the given [`ForceScheme`](@ref).
The [`DirectScheme`](@ref) path is the pair sum of the oracle and is left
bit-identical (symmetric, hence serial); [`PerParticleScheme`](@ref) assembles the
same physics per particle, which threads; [`TreeScheme`](@ref) rebuilds the tree and
the full interaction lists at every step (no particle is frozen in a fixed-step run).
"""
function verlet_forces!(f::Matrix{Float64}, st::PlasmaState, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, t::Float64, ::DirectScheme)
    return acceleration!(f, st, fm, tab, L, Q_ELEC * laser_electric_field(p, t))
end

function verlet_forces!(f::Matrix{Float64}, st::PlasmaState, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, t::Float64, s::PerParticleScheme)
    n = size(st.pos, 2)
    if par_worth(n, n)
        nt = Threads.nthreads()
        chunk = cld(n, nt)
        Threads.@threads :static for c in 1:nt
            lo = (c - 1) * chunk + 1
            hi = min(c * chunk, n)
            for i in lo:hi
                build_force_i!(f, i, t, st.pos, p, fm, tab, L, st.ions, s)
            end
        end
    else
        for i in 1:n
            build_force_i!(f, i, t, st.pos, p, fm, tab, L, st.ions, s)
        end
    end
    return f
end

function verlet_forces!(f::Matrix{Float64}, st::PlasmaState, p::PlasmaParams, fm::Forces,
                        tab::EwaldTable, L::Float64, t::Float64, s::TreeScheme)
    set_bodies!(s.tree, s.mass, s.charge, st.pos, st.ions)
    maketree!(s.tree)
    fill_lists!(s, nothing)
    n = size(st.pos, 2)
    nt = Threads.nthreads()
    # electron `i` writes only its own column of `f`, so this needs no reduction
    if par_worth(n, n)
        chunk = cld(n, nt)
        Threads.@threads :static for c in 1:nt
            lo = (c - 1) * chunk + 1
            hi = min(c * chunk, n)
            for i in lo:hi
                build_force_i!(f, i, t, st.pos, p, fm, tab, L, st.ions, s)
            end
        end
    else
        for i in 1:n
            build_force_i!(f, i, t, st.pos, p, fm, tab, L, st.ions, s)
        end
    end
    return f
end

"""
    verlet_step!(st, f, p, fm, tab, L, t, h, scheme = DirectScheme())

One leapfrog step of `verlet`, with the forces evaluated at the beginning of the
step (time `t`), as in the C. `scheme` selects how those forces are summed; the
default is the direct pair sum of the oracle.
"""
function verlet_step!(st::PlasmaState, f::Matrix{Float64}, p::PlasmaParams, fm::Forces,
                      tab::EwaldTable, L::Float64, t::Float64, h::Float64,
                      scheme::ForceScheme = DirectScheme())
    verlet_forces!(f, st, p, fm, tab, L, t, scheme)
    c2 = C_LIGHT * C_LIGHT
    coef1 = 0.5 * h
    coef2 = Q_ELEC * 0.5 * h / (C_LIGHT * M_ELEC)
    pos = st.pos
    mom = st.mom
    @inbounds for i in axes(pos, 2)
        u = Vec3(mom[1, i], mom[2, i], mom[3, i])
        acc = Vec3(f[1, i], f[2, i], f[3, i])

        um = u + coef1 * acc
        gamman = sqrt(1.0 + squared_norm(um) / c2)

        xn, yn, zn = pos[1, i], pos[2, i], pos[3, i]
        _, _, _, bx, by, bz = laser_fields(p, t, xn, yn, zn)
        bvec = Vec3(bx, by, bz)

        tvec = (coef2 / gamman) * bvec
        svec = (2.0 / (1.0 + squared_norm(tvec))) * tvec

        uprime = um + cross(um, tvec)
        up = um + cross(uprime, svec)

        unph = up + coef1 * acc
        gammanph = sqrt(1.0 + squared_norm(unph) / c2)
        coef5 = h / gammanph

        pos[1, i] = xn + coef5 * unph[1]
        pos[2, i] = yn + coef5 * unph[2]
        pos[3, i] = zn + coef5 * unph[3]

        mom[1, i] = unph[1]
        mom[2, i] = unph[2]
        mom[3, i] = unph[3]
    end
    return st
end
