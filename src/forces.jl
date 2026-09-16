# Force assembly, ported from force_background / force_laser_ele /
# force_manybody_direct_ele (ppbs.c). In plasma mode the pair force is
# `q_i q_j [ r/(r²+ε²)^{3/2} + F_corr(r) ]`, where `F_corr` comes from the
# Ewald table; displacements use the minimum image.

"""
    par_worth(population, n)

Is a force loop over `population` particles worth splitting across threads? The cost
of one particle is O(N) in the direct scheme and O(log N) with the tree, so the
population alone is the wrong criterion: at N=500 every particle is unstable at level
1 and the step is 250 000 pair evaluations, yet a population threshold of 512 keeps it
on one core. Measured (session 5): that is exactly what happened to every N≤500 run of
`benchmarks/`, which used one core out of ten.

The product `population * n` is the pair count, and 4096 is where the `@threads`
overhead (a few microseconds) stops mattering. Splitting changes no arithmetic —
electron `i` writes only its own columns and sums its own pairs in the same order — so
the direct path stays bit-identical to the oracle.
"""
@inline par_worth(population::Int, n::Int) =
    Threads.nthreads() > 1 && population * n >= 4096

"""
    min_image(d, L)

Fold one component into `[-L/2, L/2[`, with the C bounds.
"""
@inline function min_image(d::Float64, L::Float64)
    half = 0.5 * L
    while d >= half
        d -= L
    end
    while d < -half
        d += L
    end
    return d
end

"""
    pair_force_ee(pot, tab, L, xi, yi, zi, xj, yj, zj)

Electron-electron pair force (charges `q_e² = 1`), minimum image.
"""
@inline function pair_force_ee(pot::SoftCorePot, tab::EwaldTable, L::Float64,
                               xi::Float64, yi::Float64, zi::Float64,
                               xj::Float64, yj::Float64, zj::Float64)
    xij = min_image(xi - xj, L)
    yij = min_image(yi - yj, L)
    zij = min_image(zi - zj, L)
    sqradij = softened_sqdist(pot, xij, yij, zij)
    radij = sqrt(sqradij)
    coef = 1.0 / (radij * sqradij)
    fwx, fwy, fwz = ewald_force(tab, xij, yij, zij)
    return xij * coef + fwx, yij * coef + fwy, zij * coef + fwz
end

"""
    pair_force_eion(pot, tab, L, constout, xi, yi, zi, xj, yj, zj)

Electron-ion pair force with `constout = q_e q_ion`, minimum image.
"""
@inline function pair_force_eion(pot::SoftCorePot, tab::EwaldTable, L::Float64,
                                 constout::Float64,
                                 xi::Float64, yi::Float64, zi::Float64,
                                 xj::Float64, yj::Float64, zj::Float64)
    xij = min_image(xi - xj, L)
    yij = min_image(yi - yj, L)
    zij = min_image(zi - zj, L)
    radius2 = softened_sqdist(pot, xij, yij, zij)
    radius = sqrt(radius2)
    radim3 = constout / (radius2 * radius)
    fwx, fwy, fwz = ewald_force(tab, xij, yij, zij)
    return radim3 * xij + constout * fwx,
           radim3 * yij + constout * fwy,
           radim3 * zij + constout * fwz
end

"""
    pair_pot_ee(pot, tab, L, xi, yi, zi, xj, yj, zj)

Electron-electron pair potential (used by the energy diagnostics).
"""
@inline function pair_pot_ee(pot::SoftCorePot, tab::EwaldTable, L::Float64,
                             xi::Float64, yi::Float64, zi::Float64,
                             xj::Float64, yj::Float64, zj::Float64)
    xij = min_image(xi - xj, L)
    yij = min_image(yi - yj, L)
    zij = min_image(zi - zj, L)
    sqradii = softened_sqdist(pot, xij, yij, zij)
    return sqradii^(-0.5) + ewald_pot(tab, xij, yij, zij)
end

"""
    pair_pot_eion(pot, tab, L, constout, xi, yi, zi, xj, yj, zj)

Electron-ion pair potential (used by the energy diagnostics).
"""
@inline function pair_pot_eion(pot::SoftCorePot, tab::EwaldTable, L::Float64,
                               constout::Float64,
                               xi::Float64, yi::Float64, zi::Float64,
                               xj::Float64, yj::Float64, zj::Float64)
    xij = min_image(xi - xj, L)
    yij = min_image(yi - yj, L)
    zij = min_image(zi - zj, L)
    radius2 = softened_sqdist(pot, xij, yij, zij)
    radius = sqrt(radius2)
    return constout * (1.0 / radius + ewald_pot(tab, xij, yij, zij))
end

"""
    background_force!(f, st, fm, tab, L)

Ion background force. Assigns (does not accumulate) `f[:, i]`.

Electron `i` sums over every ion on its own and writes only its own column, so this
splits across threads with no reduction and **no change of summation order**: the
result is bit-identical to the serial loop, and to the oracle. It is worth splitting
because it is the larger half of a fixed-step force — measured at N=500, 5.1 ms of the
7.8 ms of a `verlet_step!`, against 2.7 ms for the (symmetric, hence serial) e-e sum.
"""
function background_force!(f::Matrix{Float64}, st::PlasmaState, fm::Forces,
                           tab::EwaldTable, L::Float64)
    n = size(st.pos, 2)
    if par_worth(n, size(st.ions, 2))
        nt = Threads.nthreads()
        chunk = cld(n, nt)
        Threads.@threads :static for c in 1:nt
            lo = (c - 1) * chunk + 1
            hi = min(c * chunk, n)
            lo <= hi && _background_chunk!(f, st, fm, tab, L, lo, hi)
        end
    else
        _background_chunk!(f, st, fm, tab, L, 1, n)
    end
    return f
end

@inline function _background_chunk!(f::Matrix{Float64}, st::PlasmaState, fm::Forces,
                                    tab::EwaldTable, L::Float64, lo::Int, hi::Int)
    constout = fm.q_ion * Q_ELEC
    pos = st.pos
    @inbounds for i in lo:hi
        xi, yi, zi = pos[1, i], pos[2, i], pos[3, i]
        fix = fiy = fiz = 0.0
        for ion in eachcol(st.ions)
            fx, fy, fz = pair_force_eion(fm.ion, tab, L, constout, xi, yi, zi,
                                         ion[1], ion[2], ion[3])
            fix += fx
            fiy += fy
            fiz += fz
        end
        f[1, i], f[2, i], f[3, i] = fix, fiy, fiz
    end
    return f
end

"""
    laser_force!(f, force)

Laser force: x component only, added to `f[:, i]`.
"""
function laser_force!(f::Matrix{Float64}, force::Float64)
    @views f[1, :] .+= force
    return f
end

"""
    ee_force!(f, st, fm, tab, L)

Direct electron-electron force, `O(N²)`, symmetric accumulation as in the C.
"""
function ee_force!(f::Matrix{Float64}, st::PlasmaState, fm::Forces,
                   tab::EwaldTable, L::Float64)
    pos = st.pos
    n = size(pos, 2)
    @inbounds for i in 1:n-1
        r = view(pos, :, i)
        for j in i+1:n
            fx, fy, fz = pair_force_ee(fm.ee, tab, L, r[1], r[2], r[3],
                                       pos[1, j], pos[2, j], pos[3, j])
            f[1, i] += fx
            f[2, i] += fy
            f[3, i] += fz
            f[1, j] -= fx
            f[2, j] -= fy
            f[3, j] -= fz
        end
    end
    return f
end

"""
    acceleration!(f, st, fm, tab, L, force_laser)

Total acceleration (`m_e = 1`), in the C order: background, laser, e-e.
"""
function acceleration!(f::Matrix{Float64}, st::PlasmaState, fm::Forces,
                       tab::EwaldTable, L::Float64, force_laser::Float64)
    background_force!(f, st, fm, tab, L)
    laser_force!(f, force_laser)
    ee_force!(f, st, fm, tab, L)
    return f
end
