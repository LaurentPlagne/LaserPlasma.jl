# Fixed-step plasma simulation loop (the C `verlet` scheme), with numerical
# outputs comparable to the oracle.

"""
    RunConfig(xfinal, nb_step, outdir)

Run request: final time, number of steps per half laser period, output directory.
"""
struct RunConfig
    xfinal::Float64
    nb_step::Int
    outdir::String
end

"""
    Checkpoint

On-disk checkpoint: the state and accumulators needed to resume a run at an outer
step boundary. The integrator is a one-step method, so this is exact.
"""
struct Checkpoint
    counter::Int
    t0::Float64
    pos::Matrix{Float64}
    mom::Matrix{Float64}
    ions::Matrix{Float64}
    jx::Float64
    jy::Float64
    jz::Float64
    fresh::Bool
    du_ref::Float64
    du_integrate::Float64
    nlines::NTuple{3,Int}
end

"""
    save_checkpoint(path, ck)

Write `ck` to `path` in the binary layout read back by [`load_checkpoint`](@ref)
(first `Int64` = counter, then the `Float64` time).
"""
function save_checkpoint(path::String, ck::Checkpoint)
    open(path, "w") do io
        write(io, Int64(ck.counter))
        write(io, ck.t0)
        write(io, Int64(size(ck.pos, 1)))
        write(io, Int64(size(ck.pos, 2)))
        write(io, ck.pos)
        write(io, ck.mom)
        write(io, Int64(size(ck.ions, 2)))
        write(io, ck.ions)
        write(io, ck.jx)
        write(io, ck.jy)
        write(io, ck.jz)
        write(io, ck.fresh)
        write(io, ck.du_ref)
        write(io, ck.du_integrate)
        for n in ck.nlines
            write(io, Int64(n))
        end
    end
    return path
end

"""
    load_checkpoint(path)

Read the checkpoint written by [`save_checkpoint`](@ref).
"""
function load_checkpoint(path::String)
    open(path, "r") do io
        counter = read(io, Int64)
        t0 = read(io, Float64)
        nr = read(io, Int64)
        nc = read(io, Int64)
        pos = Array{Float64}(undef, nr, nc)
        read!(io, pos)
        mom = Array{Float64}(undef, nr, nc)
        read!(io, mom)
        ni = Int(read(io, Int64))
        ions = Array{Float64}(undef, nr, ni)
        read!(io, ions)
        jx = read(io, Float64)
        jy = read(io, Float64)
        jz = read(io, Float64)
        fresh = read(io, Bool)
        du_ref = read(io, Float64)
        du_integrate = read(io, Float64)
        nlines = (read(io, Int64), read(io, Int64), read(io, Int64))
        return Checkpoint(counter, t0, pos, mom, ions, jx, jy, jz, fresh, du_ref,
                          du_integrate, nlines)
    end
end

"""
    truncate_lines(path, n)

Keep only the first `n` lines of a text output (used when resuming).
"""
function truncate_lines(path::String, n::Int)
    isfile(path) || return
    lines = readlines(path)
    length(lines) > n && write(path, join(lines[1:n], "\n") * "\n")
    return
end

"""
    run_multistep(p, cfg, eps; sample = 2, level_max = 15, thermostat = false,
                  checkpoint_every = 0, restart = false, use_tree = false,
                  theta = 0.5, usequad = 2, interp = :quad,
                  progress = (c, t) -> nothing)

Adaptive individual-step plasma run (the C `multistep_rec` scheme). All
diagnostics of this scheme refer to the state at the *start* of the step and are
labelled with the accumulated start time, as in the C: `energy.dat`, `power.dat`,
`deltaU` and the `pos_0001`/`mome_0001` snapshots.

With `use_tree = true`, forces come from the Barnes-Hut tree (C `do_tree = 1`):
the tree is rebuilt at each recursion level with the frozen-list semantics of
`docs/architecture/barnes-hut.md`, and the ions are bodies of the tree.

`interp` selects the stencil that freezes the stable particles inside a refined
interval: `:quad` (default, the quadratic stencil of the C), `:pwise` or
`:simple` (lower-order variants used to measure the coupling-order requirement
of the recursion, see `benchmarks/interp_order.jl`).
"""
function run_multistep(p::PlasmaParams, cfg::RunConfig, eps::Float64; sample::Int = 2,
                       level_max::Int = 15, thermostat::Bool = false,
                       checkpoint_every::Int = 0,
                       restart::Bool = false,
                       use_tree::Bool = false, theta::Float64 = 0.5, usequad::Int = 2,
                       interp::Symbol = :quad,
                       progress::Function = (counter, t0) -> nothing)
    L = box_length(p)
    tab = use_tree ? build_ewald_hessian_table(L) : build_ewald_table(L)
    scheme = use_tree ? TreeScheme(p; theta = theta, usequad = usequad) : DirectScheme()
    fm = Forces(p)
    st = PlasmaState(p)
    work = MultistepWork(nb_elec(p), level_max)
    cs = CurrentState()
    du = DeltaUState()
    ckpt_path = joinpath(cfg.outdir, "checkpoint.bin")
    mkpath(cfg.outdir)
    start_counter = 1
    t0 = 0.0
    n_ener = n_power = n_du = 0
    if restart && isfile(ckpt_path)
        ck = load_checkpoint(ckpt_path)
        st.pos .= ck.pos
        st.mom .= ck.mom
        st.ions .= ck.ions
        cs.jx, cs.jy, cs.jz, cs.fresh = ck.jx, ck.jy, ck.jz, ck.fresh
        du.total_energy_ref, du.total_integrate = ck.du_ref, ck.du_integrate
        start_counter = ck.counter + 1
        t0 = ck.t0
        n_ener, n_power, n_du = ck.nlines
        # drop any line written after the checkpoint (crash between checkpoints)
        truncate_lines(joinpath(cfg.outdir, "energy.dat"), n_ener)
        truncate_lines(joinpath(cfg.outdir, "power.dat"), n_power)
        truncate_lines(joinpath(cfg.outdir, "deltaU"), n_du)
    else
        init_plasma!(st, p; sample = sample)
    end
    delta_x = laser_half_period(p) / cfg.nb_step
    nsteps = trunc(Int, trunc(Int, cfg.xfinal) / delta_x)
    mode = start_counter > 1 ? "a" : "w"
    io_ener = open(joinpath(cfg.outdir, "energy.dat"), mode)
    io_power = open(joinpath(cfg.outdir, "power.dat"), mode)
    io_du = open(joinpath(cfg.outdir, "deltaU"), mode)
    io_pos = open(joinpath(cfg.outdir, "pos_0001"), mode)
    io_mom = open(joinpath(cfg.outdir, "mome_0001"), mode)
    ekin_ref = 1.5 * nb_elec(p) * p.kT
    for counter in start_counter:nsteps
        progress(counter, t0)
        ekin, epot, etot = energies(st, p, fm, tab, L)
        @printf(io_ener, "%.17g\t%.17g\t%.17g\t%.17g\n", t0, ekin, epot, etot)
        n_ener += 1
        measure_power!(cs, io_power, st, fm, tab, L, t0, delta_x, p)
        n_power += 1
        du_line, du_wrote = measure_deltaU!(du, io_du, counter, delta_x, cfg.nb_step,
                                            p.n_relax, p.n_increase, etot)
        du_wrote && (n_du += 1)
        # constant-temperature scheme (`multi_distribution` + `scale_veloc`):
        # the scaling factor is computed on the start-of-step state and applied
        # to the momenta after the integration, exactly as in the C
        sqalpha = 1.0
        if thermostat
            ntt = counter - 1
            ntc = ntt - 2 * cfg.nb_step * (p.n_relax + p.n_increase)
            uth = thermal_energy(st, p)
            if ntc < 0
                sqalpha = sqrt((ekin - (uth - ekin_ref)) / ekin)
            else
                quot, rem = divrem(ntc, 2 * cfg.nb_step)
                quot > 0 && rem == 0 && (sqalpha = sqrt((ekin - du_line) / ekin))
            end
        end
        if counter == nsteps
            dump_positions(io_pos, st, t0)
            dump_momenta(io_mom, st, t0)
        end
        multistep_step!(st, p, fm, tab, L, t0, delta_x, eps; level_max = level_max,
                        work = work, scheme = scheme, interp = interp)
        thermostat && (st.mom .*= sqalpha)
        wrap_in_box!(st, L)
        t0 += delta_x
        if checkpoint_every > 0 && counter % checkpoint_every == 0
            for io in (io_ener, io_power, io_du, io_pos, io_mom)
                flush(io)
            end
            save_checkpoint(ckpt_path,
                            Checkpoint(counter, t0, copy(st.pos), copy(st.mom),
                                       copy(st.ions), cs.jx, cs.jy, cs.jz, cs.fresh,
                                       du.total_energy_ref, du.total_integrate,
                                       (n_ener, n_power, n_du)))
        end
    end
    close(io_ener)
    close(io_power)
    close(io_du)
    close(io_pos)
    close(io_mom)
    open(joinpath(cfg.outdir, "ion_pos.dat"), "w") do io
        dump_ions(io, st)
    end
    return st
end

"""
    run_verlet(p, cfg; sample = 2)

Fixed-step plasma simulation loop (the C `verlet` scheme), with numerical outputs
comparable to the oracle. Records the state at the **end** of each step.
"""
function run_verlet(p::PlasmaParams, cfg::RunConfig; sample::Int = 2)
    L = box_length(p)
    tab = build_ewald_table(L)
    fm = Forces(p)
    st = PlasmaState(p)
    init_plasma!(st, p; sample = sample)
    f = zeros(3, nb_elec(p))
    delta_x = laser_half_period(p) / cfg.nb_step
    nsteps = trunc(Int, trunc(Int, cfg.xfinal) / delta_x)
    mkpath(cfg.outdir)
    io_ener = open(joinpath(cfg.outdir, "energy.dat"), "w")
    io_power = open(joinpath(cfg.outdir, "power.dat"), "w")
    io_du = open(joinpath(cfg.outdir, "deltaU"), "w")
    io_pos = open(joinpath(cfg.outdir, "pos_0001"), "w")
    io_mom = open(joinpath(cfg.outdir, "mome_0001"), "w")
    cs = CurrentState()
    du = DeltaUState()
    # the C accumulates time step by step (`x_start=x_end`), not by product
    t0 = 0.0
    for counter in 1:nsteps
        verlet_step!(st, f, p, fm, tab, L, t0, delta_x)
        t1 = t0 + delta_x
        t0 = t1
        ekin, epot, etot = energies(st, p, fm, tab, L)
        @printf(io_ener, "%.17g\t%.17g\t%.17g\t%.17g\n", t1, ekin, epot, etot)
        measure_power!(cs, io_power, st, fm, tab, L, t1, delta_x, p)
        measure_deltaU!(du, io_du, counter, delta_x, cfg.nb_step, p.n_relax, p.n_increase,
                        etot)
        # the C records the state before `back_in_cell` and wraps afterwards only
        if counter == nsteps
            dump_positions(io_pos, st, t1)
            dump_momenta(io_mom, st, t1)
        end
        wrap_in_box!(st, L)
    end
    close(io_ener)
    close(io_power)
    close(io_du)
    close(io_pos)
    close(io_mom)
    open(joinpath(cfg.outdir, "ion_pos.dat"), "w") do io
        dump_ions(io, st)
    end
    return st
end
