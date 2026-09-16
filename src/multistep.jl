# Individual adaptive time-step integrator, ported from multistep.c.
#
# The scheme: one full RK4 step over the interval is the reference; two half
# RK4 steps give a finer estimate. Particles whose two estimates differ by more
# than the tolerance stay "unstable" and their sub-intervals are split again,
# while particles that became "stable" have their trajectory frozen by
# quadratic interpolation and are no longer integrated. The cost concentrates
# on close encounters, which is what the original method was designed for.
#
# In direct mode (`do_tree=0`) forces are evaluated per particle, in the C
# order: electron-electron, laser, then ions.

"""
    build_velocity_i!(vel, mom, i, relat)

Velocity `v = u/γ` of electron `i` (or `v = u` when `relat` is false).
"""
@inline function build_velocity_i!(vel::Matrix{Float64}, mom::Matrix{Float64}, i::Int,
                                   relat::Bool)
    ux = mom[1, i]
    uy = mom[2, i]
    uz = mom[3, i]
    if relat
        c2 = C_LIGHT * C_LIGHT
        u2 = ux * ux + uy * uy + uz * uz
        gammam1 = 1.0 / sqrt((u2 + c2) / c2)
        vel[1, i] = gammam1 * ux
        vel[2, i] = gammam1 * uy
        vel[3, i] = gammam1 * uz
    else
        vel[1, i] = ux
        vel[2, i] = uy
        vel[3, i] = uz
    end
    return vel
end

"""
    euler1_rec!(stable, t, dt, deb, vel, force, kpos, kmom, p, fm, tab, L, ions)

First Euler stage of RK4: velocities, forces and increments `k` for the unstable
particles.
"""
function euler1_rec!(stable::AbstractVector{Bool}, t::Float64, dt::Float64,
                     deb::ParticleSet, vel::Matrix{Float64}, force::Matrix{Float64},
                     kpos::ParticleSet, kmom::ParticleSet, p::PlasmaParams, fm::Forces,
                     tab::EwaldTable, L::Float64, ions::Matrix{Float64},
                     scheme::ForceScheme)
    prepare_stage!(scheme, deb, ions)
    n = size(deb.pos, 2)
    nt = Threads.nthreads()
    # Each unstable electron writes only its own columns of kpos/kmom (and of
    # `force`), so the loop parallelizes without a reduction.
    if par_worth(count(iszero, stable), n)
        chunk = cld(n, nt)
        Threads.@threads :static for c in 1:nt
            lo = (c - 1) * chunk + 1
            hi = min(c * chunk, n)
            lo <= hi && _euler1_chunk!(stable, t, dt, deb, vel, force, kpos, kmom, p, fm,
                                       tab, L, ions, scheme, lo, hi)
        end
    else
        _euler1_chunk!(stable, t, dt, deb, vel, force, kpos, kmom, p, fm, tab, L, ions,
                       scheme, 1, n)
    end
    return kpos, kmom
end

@inline function _euler1_chunk!(stable::AbstractVector{Bool}, t::Float64, dt::Float64,
                                deb::ParticleSet, vel::Matrix{Float64},
                                force::Matrix{Float64}, kpos::ParticleSet,
                                kmom::ParticleSet, p::PlasmaParams, fm::Forces,
                                tab::EwaldTable, L::Float64, ions::Matrix{Float64},
                                scheme::ForceScheme, lo::Int, hi::Int)
    @inbounds for i in lo:hi
        stable[i] && continue
        build_velocity_i!(vel, deb.mom, i, p.d_relat)
        build_force_i!(force, i, t, deb.pos, p, fm, tab, L, ions, scheme)
        for c in 1:3
            kpos.pos[c, i] = dt * vel[c, i]
            kmom.pos[c, i] = (dt / M_ELEC) * force[c, i]
        end
    end
    return nothing
end

"""
    euler2_rec!(stable, coef, deb, kpos, kmom, out)

`out = deb + coef * k`, for the unstable particles.
"""
function euler2_rec!(stable::AbstractVector{Bool}, coef::Float64, deb::ParticleSet,
                     kpos::ParticleSet, kmom::ParticleSet, out::ParticleSet)
    @inbounds for i in axes(deb.pos, 2)
        stable[i] && continue
        for c in 1:3
            out.pos[c, i] = deb.pos[c, i] + coef * kpos.pos[c, i]
            out.mom[c, i] = deb.mom[c, i] + coef * kmom.pos[c, i]
        end
    end
    return out
end

# Interpolation stencils used to freeze the stable particles inside an interval.
# Each entry gives the coefficients of `deb` (t), `mid` (t+dt/2) and `fin` (t+dt)
# for the values at `t+dt/4` (first triple) and `t+3dt/4` (second triple); the
# common 1/8 factor is applied last.
#  - :quad   quadratic Lagrange through the three RK4 nodes (the C scheme),
#  - :pwise  piecewise linear `deb-mid` then `mid-fin`,
#  - :simple straight line `deb-fin` (ignores `mid`).
# The two lower-order stencils exist to measure the coupling-order requirement
# of the recursion (see benchmarks/interp_order.jl).
const INTERP_COEFS = Dict{Symbol,NTuple{2,NTuple{3,Float64}}}(
    :quad => ((3.0, 6.0, -1.0), (-1.0, 6.0, 3.0)),
    :pwise => ((4.0, 4.0, 0.0), (0.0, 4.0, 4.0)),
    :simple => ((6.0, 0.0, 2.0), (2.0, 0.0, 6.0)),
)

@inline function _interp_lerp!(stable::AbstractVector{Bool}, deb::ParticleSet,
                               mid::ParticleSet, fin::ParticleSet, out::ParticleSet,
                               ca::Float64, cb::Float64, cc::Float64)
    ooh = 1.0 / 8.0
    @inbounds for i in axes(deb.pos, 2)
        stable[i] || continue
        for c in 1:3
            out.pos[c, i] = ooh * (ca * deb.pos[c, i] + cb * mid.pos[c, i] +
                                   cc * fin.pos[c, i])
            out.mom[c, i] = ooh * (ca * deb.mom[c, i] + cb * mid.mom[c, i] +
                                   cc * fin.mom[c, i])
        end
    end
    return out
end

"""
    interp_rec_fh!(stable, deb, mid, fin, out; interp = :quad)

Values at `t+dt/4` for the stable particles, by interpolation through
`(t, deb)`, `(t+dt/2, mid)` and `(t+dt, fin)`. `interp` selects the stencil
(`:quad`, `:pwise` or `:simple`, see `INTERP_COEFS`); the default is the
quadratic stencil of the C.
"""
function interp_rec_fh!(stable::AbstractVector{Bool}, deb::ParticleSet, mid::ParticleSet,
                        fin::ParticleSet, out::ParticleSet; interp::Symbol = :quad)
    ca, cb, cc = INTERP_COEFS[interp][1]
    return _interp_lerp!(stable, deb, mid, fin, out, ca, cb, cc)
end

"""
    interp_rec_sh!(stable, deb, mid, fin, out; interp = :quad)

Same as [`interp_rec_fh!`](@ref) at `t+3dt/4` (stencil reversed, as in the C).
"""
function interp_rec_sh!(stable::AbstractVector{Bool}, deb::ParticleSet, mid::ParticleSet,
                        fin::ParticleSet, out::ParticleSet; interp::Symbol = :quad)
    ca, cb, cc = INTERP_COEFS[interp][2]
    return _interp_lerp!(stable, deb, mid, fin, out, ca, cb, cc)
end

"""
    compare_rec!(stable, newstable, eps, ref, result)

Position-based stability test: `result` (two half steps) against `ref` (known
full step). Returns the number of still unstable particles; `newstable` is
updated in place.
"""
function compare_rec!(stable::AbstractVector{Bool}, newstable::AbstractVector{Bool},
                      eps::Float64, ref::ParticleSet, result::ParticleSet)
    population = 0
    eps2 = eps * eps
    @inbounds for i in axes(stable, 1)
        if stable[i]
            newstable[i] = true
            continue
        end
        dx = result.pos[1, i] - ref.pos[1, i]
        dy = result.pos[2, i] - ref.pos[2, i]
        dz = result.pos[3, i] - ref.pos[3, i]
        d2 = dx * dx + dy * dy + dz * dz
        if d2 > eps2
            population += 1
            newstable[i] = false
        else
            newstable[i] = true
        end
    end
    return population
end

"""
    RK4Work(n)

Work arrays of one RK4 call, reused between the two half steps.
"""
struct RK4Work
    p1::ParticleSet
    p2::ParticleSet
    p3::ParticleSet
    kp1::ParticleSet
    kp2::ParticleSet
    kp3::ParticleSet
    kp4::ParticleSet
    km1::ParticleSet
    km2::ParticleSet
    km3::ParticleSet
    km4::ParticleSet
    vel::Matrix{Float64}
    force::Matrix{Float64}
end

function RK4Work(n::Int)
    s() = particle_set(n)
    return RK4Work(s(), s(), s(), s(), s(), s(), s(), s(), s(), s(), s(), zeros(3, n),
                   zeros(3, n))
end

"""
    StepBuffers(n)

Scratch sets of one `step_rec!` call. One set per recursion depth is enough:
the recursion is depth-first, so at most one call per depth is live.
"""
struct StepBuffers
    mid_fh::ParticleSet
    mid_sh::ParticleSet
    int::ParticleSet
    out::ParticleSet
    old::ParticleSet
    in2::ParticleSet
    newstable::BitVector
end

StepBuffers(n::Int) = StepBuffers(particle_set(n), particle_set(n), particle_set(n),
                                  particle_set(n), particle_set(n), particle_set(n),
                                  falses(n))

"""
    MultistepWork(n, level_max = 15)

All buffers of one adaptive step, reusable across steps at constant `N`. The
entry state is `deb`; `mid_known`/`fin_known`/`ref` are the full-step reference
scratch. `levels` is indexed by recursion depth.
"""
struct MultistepWork
    deb::ParticleSet
    mid_known::ParticleSet
    fin_known::ParticleSet
    ref::ParticleSet
    out::ParticleSet
    stable::BitVector
    rk::RK4Work
    levels::Vector{StepBuffers}
    # Per-step recursion profile (same diagnostic as the C `level = k ; population
    # = n` printout): `visits[k]` counts the recursion nodes at level k and
    # `populations[k]` the total unstable-particle visits there. Reset by
    # `multistep_step!`.
    visits::Vector{Int}
    populations::Vector{Int}
    # Per-particle deepest level of the last step (level at which the particle
    # passed the stability test and was frozen), 0 if the step failed.
    particle_level::Vector{Int}
    # DFS preorder log of the recursion nodes of the last step: level and
    # population of unstable particles. The preorder sequence and the levels
    # reconstruct the dyadic interval tree (figures of method).
    node_levels::Vector{Int}
    node_pops::Vector{Int}
    # Optional per-particle trace (figures of method): for the indices with
    # `record[i] = true`, `traces[i]` accumulates the accepted RK4 nodes of the
    # particle as (t, x, y, z, level). Not reset by `multistep_step!`; the caller
    # clears it.
    record::BitVector
    traces::Vector{Vector{NTuple{5,Float64}}}
end

function MultistepWork(n::Int, level_max::Int = 15)
    levels = [StepBuffers(n) for _ in 1:level_max+2]
    return MultistepWork(particle_set(n), particle_set(n), particle_set(n),
                         particle_set(n), particle_set(n), falses(n), RK4Work(n), levels,
                         zeros(Int, level_max + 2), zeros(Int, level_max + 2),
                         zeros(Int, n), Int[], Int[], falses(n),
                         [NTuple{5,Float64}[] for _ in 1:n])
end

"""
    rk4_rec!(t1, dt, stable, deb, mid_known, fin_known, out, ws, p, fm, tab, L, ions)

One RK4 step of length `dt` from `deb`. For stable particles the stage values are
pinned to the known trajectory (`mid_known` at `t+dt/2`, `fin_known` at `t+dt`);
the output is `out`.
"""
function rk4_rec!(t1::Float64, dt::Float64, stable::AbstractVector{Bool},
                  deb::ParticleSet, mid_known::ParticleSet, fin_known::ParticleSet,
                  out::ParticleSet, ws::RK4Work, p::PlasmaParams, fm::Forces,
                  tab::EwaldTable, L::Float64, ions::Matrix{Float64},
                  scheme::ForceScheme)
    # k1: stage at t1 is pinned to mid_known for stable particles (as in the C)
    euler1_rec!(stable, t1, dt, deb, ws.vel, ws.force, ws.kp1, ws.km1, p, fm, tab, L, ions,
                scheme)
    euler2_rec!(stable, 0.5, deb, ws.kp1, ws.km1, ws.p1)
    copy_stable!(stable, ws.p1, mid_known)

    # k2: stage at t1+dt/2
    timei = t1 + 0.5 * dt
    euler1_rec!(stable, timei, dt, ws.p1, ws.vel, ws.force, ws.kp2, ws.km2, p, fm, tab, L,
                ions, scheme)
    euler2_rec!(stable, 0.5, deb, ws.kp2, ws.km2, ws.p2)
    copy_stable!(stable, ws.p2, mid_known)

    # k3: stage at t1+dt/2, pinned to fin_known for stable particles
    euler1_rec!(stable, timei, dt, ws.p2, ws.vel, ws.force, ws.kp3, ws.km3, p, fm, tab, L,
                ions, scheme)
    euler2_rec!(stable, 1.0, deb, ws.kp3, ws.km3, ws.p3)
    copy_stable!(stable, ws.p3, fin_known)

    # k4: stage at t1+dt
    euler1_rec!(stable, t1 + dt, dt, ws.p3, ws.vel, ws.force, ws.kp4, ws.km4, p, fm, tab,
                L, ions, scheme)

    # summation, in the C accumulation order
    oos = 1.0 / 6.0
    oot = 1.0 / 3.0
    @inbounds for i in axes(deb.pos, 2)
        stable[i] && continue
        for c in 1:3
            v = deb.pos[c, i] + oos * ws.kp1.pos[c, i]
            v += oot * ws.kp2.pos[c, i]
            v += oot * ws.kp3.pos[c, i]
            v += oos * ws.kp4.pos[c, i]
            out.pos[c, i] = v
            w = deb.mom[c, i] + oos * ws.km1.pos[c, i]
            w += oot * ws.km2.pos[c, i]
            w += oot * ws.km3.pos[c, i]
            w += oos * ws.km4.pos[c, i]
            out.mom[c, i] = w
        end
    end
    copy_stable!(stable, out, fin_known)
    return out
end

"""
    step_rec!(eps, t1, dt, stable, deb, mid_known, fin_known, out, ws, p, fm, tab, L,
              ions, level, level_max, failure, scheme; interp = :quad)

One recursive step: two half RK4 steps over `[t1, t1+dt]` compared with the known
full step `fin_known`; particles that fail are refined by recursion. `interp`
selects the stencil used to freeze the stable particles (see `INTERP_COEFS`).
"""
function step_rec!(eps::Float64, t1::Float64, dt::Float64,
                   stable::AbstractVector{Bool}, deb::ParticleSet,
                   mid_known::ParticleSet, fin_known::ParticleSet, out::ParticleSet,
                   ws::MultistepWork, p::PlasmaParams, fm::Forces, tab::EwaldTable,
                   L::Float64, ions::Matrix{Float64}, level::Base.RefValue{Int},
                   level_max::Int, failure::Base.RefValue{Bool}, scheme::ForceScheme;
                   interp::Symbol = :quad)
    level[] += 1
    level[] > level_max && (failure[] = true)
    push!(ws.node_levels, level[])
    push!(ws.node_pops, 0)
    node_i = length(ws.node_levels)
    rk = ws.rk
    buf = ws.levels[level[]+1]

    dts2 = dt * 0.5
    timei = t1 + dts2

    prepare_level!(scheme, stable)

    interp_rec_fh!(stable, deb, mid_known, fin_known, buf.mid_fh; interp = interp)
    interp_rec_sh!(stable, deb, mid_known, fin_known, buf.mid_sh; interp = interp)

    rk4_rec!(t1, dts2, stable, deb, buf.mid_fh, mid_known, buf.int, rk, p, fm, tab, L, ions,
             scheme)
    rk4_rec!(timei, dts2, stable, buf.int, buf.mid_sh, fin_known, buf.out, rk, p, fm, tab,
             L, ions, scheme)

    population = compare_rec!(stable, buf.newstable, eps, fin_known, buf.out)
    ws.visits[level[]] += 1
    ws.populations[level[]] += population
    ws.node_pops[node_i] = population
    @inbounds for i in eachindex(ws.stable)
        if !stable[i] && buf.newstable[i]
            ws.particle_level[i] = level[]
            if ws.record[i]
                push!(ws.traces[i], (t1 + dts2, buf.int.pos[1, i], buf.int.pos[2, i],
                                     buf.int.pos[3, i], Float64(level[])))
                push!(ws.traces[i], (t1 + dt, buf.out.pos[1, i], buf.out.pos[2, i],
                                     buf.out.pos[3, i], Float64(level[])))
            end
        end
    end
    # In the C, particles that just became stable are appended to the history
    # arrays here (yp/kount_i); only the figures of method use them.

    if !failure[] && population != 0
        interp_rec_fh!(buf.newstable, deb, buf.int, buf.out, buf.mid_fh; interp = interp)
        interp_rec_sh!(buf.newstable, deb, buf.int, buf.out, buf.mid_sh; interp = interp)

        epsnew = eps * 0.5
        copy_particles!(buf.old, buf.out)
        step_rec!(epsnew, t1, dts2, buf.newstable, deb, buf.mid_fh, buf.int, buf.in2, ws,
                  p, fm, tab, L, ions, level, level_max, failure, scheme; interp = interp)
        step_rec!(epsnew, timei, dts2, buf.newstable, buf.in2, buf.mid_sh, buf.old, out,
                  ws, p, fm, tab, L, ions, level, level_max, failure, scheme; interp = interp)
    else
        copy_particles!(out, buf.out)
    end

    level[] -= 1
    return out
end

"""
    multistep_step!(st, p, fm, tab, L, t1, dt, eps; level_max = 15, work = nothing,
                    scheme = DirectScheme(), interp = :quad)

One outer step `[t1, t1+dt]` with the individual adaptive scheme. Mirrors
`multistep_rec`, including the failure fallback onto two half steps. `interp`
selects the interpolation stencil of the frozen particles (see
`INTERP_COEFS`); the default is the quadratic stencil of the C.
"""
function multistep_step!(st::PlasmaState, p::PlasmaParams, fm::Forces, tab::EwaldTable,
                         L::Float64, t1::Float64, dt::Float64, eps::Float64;
                         level_max::Int = 15, work::Union{Nothing,MultistepWork} = nothing,
                         scheme::ForceScheme = DirectScheme(), interp::Symbol = :quad)
    ws = work === nothing ? MultistepWork(size(st.pos, 2), level_max) : work
    copyto!(ws.deb.pos, st.pos)
    copyto!(ws.deb.mom, st.mom)
    fill!(ws.stable, false)
    fill!(ws.visits, 0)
    fill!(ws.populations, 0)
    fill!(ws.particle_level, 0)
    empty!(ws.node_levels)
    empty!(ws.node_pops)
    ions = st.ions
    prepare_step!(scheme, ws.deb, ions)

    # full-step RK4 reference (all particles unstable)
    rk4_rec!(t1, dt, ws.stable, ws.deb, ws.mid_known, ws.fin_known, ws.ref, ws.rk, p, fm,
             tab, L, ions, scheme)

    failure = Ref(false)
    level = Ref(0)
    step_rec!(eps * dt, t1, dt, ws.stable, ws.deb, ws.mid_known, ws.ref, ws.out, ws, p,
              fm, tab, L, ions, level, level_max, failure, scheme; interp = interp)

    if !failure[]
        copyto!(st.pos, ws.out.pos)
        copyto!(st.mom, ws.out.mom)
    else
        multistep_step!(st, p, fm, tab, L, t1, dt * 0.5, eps; level_max = level_max,
                        work = ws, scheme = scheme, interp = interp)
        multistep_step!(st, p, fm, tab, L, t1 + dt * 0.5, dt * 0.5, eps;
                        level_max = level_max, work = ws, scheme = scheme, interp = interp)
    end
    return st
end
