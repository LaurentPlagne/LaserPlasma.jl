# Diagnostics: energies (make_ekin / make_epot_* from ppbs.c) and full-precision
# numerical outputs.

"""
    kinetic_energy(st, p)

Relativistic (or classical) kinetic energy `make_ekin` of the current momenta.
"""
function kinetic_energy(st::PlasmaState, p::PlasmaParams)
    c2 = C_LIGHT * C_LIGHT
    ekin = 0.0
    @inbounds for i in axes(st.mom, 2)
        vx, vy, vz = st.mom[1, i], st.mom[2, i], st.mom[3, i]
        veloc2 = vx * vx + vy * vy + vz * vz
        ekin += p.d_relat ? M_ELEC * c2 * (sqrt((veloc2 + c2) / c2) - 1.0) :
                0.5 * M_ELEC * veloc2
    end
    return ekin
end

"""
    potential_energy_bkg(st, fm, tab, L)

Electron-ion potential energy (C `make_epot_bkg`).
"""
function potential_energy_bkg(st::PlasmaState, fm::Forces, tab::EwaldTable, L::Float64)
    constout = fm.q_ion * Q_ELEC
    epot = 0.0
    @inbounds for r in eachcol(st.pos), ion in eachcol(st.ions)
        epot += pair_pot_eion(fm.ion, tab, L, constout, r[1], r[2], r[3],
                              ion[1], ion[2], ion[3])
    end
    return epot
end

"""
    potential_energy_ee(st, fm, tab, L)

Electron-electron potential energy, direct `O(N²)` sum (C `make_epot_ele`).
"""
function potential_energy_ee(st::PlasmaState, fm::Forces, tab::EwaldTable, L::Float64)
    pos = st.pos
    n = size(pos, 2)
    epot = 0.0
    @inbounds for i in 1:n-1
        r = view(pos, :, i)
        for j in i+1:n
            epot += pair_pot_ee(fm.ee, tab, L, r[1], r[2], r[3],
                                pos[1, j], pos[2, j], pos[3, j])
        end
    end
    return epot
end

"""
    thermal_energy(st, p)

Thermal energy measured by the C `fit_gauss`: root-mean-square velocity
dispersions `v = u/γ` of the current state, then
`Uth = 0.5 N (σx² + σy² + σz²)`. The Gaussian of the C only feeds a dump.
"""
function thermal_energy(st::PlasmaState, p::PlasmaParams)
    c2 = C_LIGHT * C_LIGHT
    n = size(st.mom, 2)
    sx = sy = sz = sx2 = sy2 = sz2 = 0.0
    @inbounds for i in 1:n
        ux = st.mom[1, i]
        uy = st.mom[2, i]
        uz = st.mom[3, i]
        u2 = ux * ux + uy * uy + uz * uz
        gamma = sqrt((u2 + c2) / c2)
        vx = ux / gamma
        vy = uy / gamma
        vz = uz / gamma
        sx += vx
        sy += vy
        sz += vz
        sx2 += vx * vx
        sy2 += vy * vy
        sz2 += vz * vz
    end
    dn = Float64(n)
    ax = sx / dn
    ay = sy / dn
    az = sz / dn
    stdx = sqrt(sx2 / dn - ax * ax)
    stdy = sqrt(sy2 / dn - ay * ay)
    stdz = sqrt(sz2 / dn - az * az)
    return (n * 0.5) * (stdx * stdx + stdy * stdy + stdz * stdz)
end

"""
    energies(st, p, fm, tab, L)

Return `(ekin, epot, etot)`; `epot = epot_bkg + epot_ee`.
"""
function energies(st::PlasmaState, p::PlasmaParams, fm::Forces, tab::EwaldTable,
                  L::Float64)
    ekin = kinetic_energy(st, p)
    epot = potential_energy_bkg(st, fm, tab, L) + potential_energy_ee(st, fm, tab, L)
    return ekin, epot, ekin + epot
end

"""
    CurrentState()

Current accumulator of `make_field_ion_dire`: the C `j_elec`, which integrates
the mean field seen by the ions (polarization term).
"""
mutable struct CurrentState
    jx::Float64
    jy::Float64
    jz::Float64
    fresh::Bool
end

CurrentState() = CurrentState(0.0, 0.0, 0.0, false)

"""
    measure_power!(cs, io, st, fm, tab, L, t, delta_x, p)

One `make_field_ion_dire` step: write time, Eion, `j_elec`, `jelec` and the two
powers `P`, and advance the integrated current `cs`.
"""
function measure_power!(cs::CurrentState, io::IO, st::PlasmaState, fm::Forces,
                        tab::EwaldTable, L::Float64, t::Float64, delta_x::Float64,
                        p::PlasmaParams)
    ftotbx = ftotby = ftotbz = 0.0
    ftotlx = ftotly = ftotlz = 0.0
    fex = fey = fez = 0.0
    @inbounds for ion in eachcol(st.ions)
        xj, yj, zj = ion[1], ion[2], ion[3]
        for r in eachcol(st.pos)
            fx, fy, fz = pair_force_eion(fm.ion, tab, L, fm.q_ion * Q_ELEC, r[1], r[2],
                                         r[3], xj, yj, zj)
            ftotbx -= fx
            ftotby -= fy
            ftotbz -= fz
        end
        fex, fey, fez = -laser_electric_field(p, t), 0.0, 0.0
        ftotlx -= fex
        ftotly -= fey
        ftotlz -= fez
    end
    norma = 1.0 / p.nb_ion
    Eionx = (ftotbx + ftotlx) * norma
    Eiony = (ftotby + ftotly) * norma
    Eionz = (ftotbz + ftotlz) * norma

    jelecx = jelecy = jelecz = 0.0
    @inbounds for r in eachcol(st.mom)
        jelecx -= r[1]
        jelecy -= r[2]
        jelecz -= r[3]
    end
    half = 0.5 * L
    volume = 8.0 * half * half * half
    jelecx /= volume
    jelecy /= volume
    jelecz /= volume

    if !cs.fresh
        cs.jx, cs.jy, cs.jz = jelecx, jelecy, jelecz
        cs.fresh = true
    else
        plasmap2 = 3.0 / (p.r_ws * p.r_ws * p.r_ws)
        norma2 = delta_x * plasmap2 / (4.0 * π)
        cs.jx += norma2 * Eionx
        cs.jy += norma2 * Eiony
        cs.jz += norma2 * Eionz
    end

    power1 = jelecx * fex + jelecy * fey + jelecz * fez
    power2 = cs.jx * fex + cs.jy * fey + cs.jz * fez
    @printf(io, "%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            t, Eionx, Eiony, Eionz, cs.jx, cs.jy, cs.jz, jelecx, jelecy, jelecz, power1,
            power2)
    return io
end

"""
    DeltaUState()

`deltaU` accumulator: reference energy and running sum, measured every
`2 nb_step` steps after the pulse ramp (`each_cycle` in the C).
"""
mutable struct DeltaUState
    total_energy_ref::Float64
    total_integrate::Float64
end

DeltaUState() = DeltaUState(0.0, 0.0)

"""
    measure_deltaU!(acc, io, counter, delta_x, nb_step, n_relax, n_increase, etot)

Return `(delta_u, wrote)`, where `delta_u` is the value written for this step
(`0.0` otherwise) and `wrote` tells whether a line was written. The thermostat
removes `delta_u` from the momenta, the run loop counts the lines.
"""
function measure_deltaU!(acc::DeltaUState, io::IO, counter::Int, delta_x::Float64,
                         nb_step::Int, n_relax::Int, n_increase::Int, etot::Float64)
    ntt = counter - 1
    ntc = ntt - 2 * nb_step * (n_relax + n_increase)
    ntc < 0 && return 0.0, false
    quot, rem = divrem(ntc, 2 * nb_step)
    rem != 0 && return 0.0, false
    t = ntt * delta_x
    if quot == 0
        acc.total_energy_ref = etot
        acc.total_integrate = etot
        @printf(io, "%.17g\t%.17g\t%.17g\n", t, 0.0, acc.total_integrate)
        return 0.0, true
    end
    delta_u = etot - acc.total_energy_ref
    acc.total_integrate += delta_u
    @printf(io, "%.17g\t%.17g\t%.17g\n", t, delta_u, acc.total_integrate)
    return delta_u, true
end

"""
    dump_positions(io, st, t)

Write the electron positions as `t x y z` lines (full precision).
"""
function dump_positions(io::IO, st::PlasmaState, t::Float64)
    @inbounds for c in eachcol(st.pos)
        @printf(io, "%.17g\t%.17g\t%.17g\t%.17g\n", t, c[1], c[2], c[3])
    end
    return io
end

"""
    dump_momenta(io, st, t)

Write the electron momenta as `t ux uy uz` lines (full precision).
"""
function dump_momenta(io::IO, st::PlasmaState, t::Float64)
    @inbounds for c in eachcol(st.mom)
        @printf(io, "%.17g\t%.17g\t%.17g\t%.17g\n", t, c[1], c[2], c[3])
    end
    return io
end

"""
    dump_ions(io, st)

Write the fixed ion positions as `x y z` lines (full precision).
"""
function dump_ions(io::IO, st::PlasmaState)
    @inbounds for c in eachcol(st.ions)
        @printf(io, "%.17g\t%.17g\t%.17g\n", c[1], c[2], c[3])
    end
    return io
end
