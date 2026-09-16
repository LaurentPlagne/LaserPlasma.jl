# method3d: figure of method — selected trajectories in 3D with the actual
# computation points (spheres) and the interpolated stretches (dashed), plus the
# dyadic subdivision tree of the most refined outer step.
#
# Two quiet particles and two violent ones are selected automatically from the
# recursion levels: the spheres are the accepted RK4 nodes (their density is the
# local time step), the dashes between them are the frozen/interpolated parts.
#
# Run: julia --project=. figures/method3d.jl

using LaserPlasma
using CairoMakie

const PC = LaserPlasma
include(joinpath(@__DIR__, "common.jl"))

function main()
    N, F0 = 50, 79.0
    nsteps = 8
    tol = 0.005
    p = PC.PlasmaParams(N, 1, 1.442, 6.935, 0.05, 0.05, F0, 3.0, 10, 2, 2, true)
    L = PC.box_length(p)
    tab = PC.build_ewald_table(L)
    fm = PC.Forces(p)
    st = PC.PlasmaState(p)
    PC.init_plasma!(st, p; sample = 2)
    w = PC.MultistepWork(N, 15)
    fill!(w.record, true)
    dt = (π / p.ω) / 8
    t0 = 8.37758

    maxlev = zeros(Int, N)
    traj = zeros(3, N, nsteps + 1)
    traj[:, :, 1] .= st.pos
    best_nodes = Int[]
    best_pops = Int[]
    best_time = t0
    for s in 1:nsteps
        PC.multistep_step!(st, p, fm, tab, L, t0 + (s - 1) * dt, dt, tol; work = w)
        traj[:, :, s + 1] .= st.pos
        maxlev .= max.(maxlev, w.particle_level)
        if sum(w.node_pops) >= sum(best_pops)
            best_nodes = copy(w.node_levels)
            best_pops = copy(w.node_pops)
            best_time = t0 + (s - 1) * dt
        end
    end

    order = sortperm(maxlev)
    quiet = order[1:2]
    violent = order[end-1:end]
    sel = vcat(quiet, violent)
    println("quiet: ", quiet, " (level ", maxlev[quiet], ")  strongly deflected: ", violent,
            " (level ", maxlev[violent], ")")

    tw = t0 + (nsteps - 2) * dt
    lvl_colors = [:gray35, :dodgerblue, :seagreen, :orange, :crimson, :purple]

    fig = Figure(size = (1240, 580))
    ax3 = Axis3(fig[1, 1], title = "trajectories (last 2 steps) — spheres: computed points; dashes: interpolation",
                xlabel = "x", ylabel = "y", zlabel = "z", aspect = (1, 1, 0.8))
    j0 = nsteps - 2 + 1   # state at the start of the displayed window
    for i in 1:N
        scatter!(ax3, traj[1, i, j0], traj[2, i, j0], traj[3, i, j0],
                 color = (:gray, 0.25), markersize = 2)
    end
    for i in sel
        pts = [pt for pt in w.traces[i] if pt[1] >= tw]
        isempty(pts) && continue
        lw = i in violent ? 2.0 : 3.0
        lines!(ax3, [p[2] for p in pts], [p[3] for p in pts], [p[4] for p in pts],
               color = (:black, 0.25), linewidth = 0.8)
        for k in 1:(length(pts) - 1)
            a, b = pts[k], pts[k + 1]
            lines!(ax3, [a[2], b[2]], [a[3], b[3]], [a[4], b[4]],
                   color = lvl_colors[clamp(round(Int, b[5]), 1, length(lvl_colors))],
                   linestyle = :dash, linewidth = lw, alpha = 0.95)
        end
        for pt in pts
            lev = clamp(round(Int, pt[5]), 1, length(lvl_colors))
            scatter!(ax3, pt[2], pt[3], pt[4], color = lvl_colors[lev],
                     marker = :circle, markersize = 4 + 2.2 * lev)
        end
    end

    ax2 = Axis(fig[1, 2], title = "subdivision of the most refined step\n(at t = " *
                                string(round(best_time, digits = 2)) * ")",
               xlabel = "t − t₀", ylabel = "level", yticks = 1:6)
    nodes = Tuple{Int,Float64,Float64,Int}[]
    function build!(i, level, tstart, dt_node)
        push!(nodes, (level, tstart, dt_node, best_pops[i]))
        i += 1
        if i <= length(best_nodes) && best_nodes[i] == level + 1
            i = build!(i, level + 1, tstart, dt_node / 2)
            i = build!(i, level + 1, tstart + dt_node / 2, dt_node / 2)
        end
        return i
    end
    build!(1, best_nodes[1], 0.0, dt)
    for (lev, ts, h, pop) in nodes
        poly!(ax2, Rect(ts, lev - 0.28, h, 0.56), color = pop,
              colormap = Reverse(:batlow), colorrange = (0, maximum(best_pops)))
        h > 0.15 * dt && text!(ax2, ts + h / 2, lev, text = string(pop),
                               align = (:center, :center), fontsize = 8)
    end
    Colorbar(fig[1, 3], limits = (0, maximum(best_pops)), colormap = Reverse(:batlow),
             label = "unstable particles")

    save_fig(fig, "method3d")
    return fig
end

main()
