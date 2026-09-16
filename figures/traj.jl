# traj79 / traj157: electron trajectories in the (x, y) plane over the ten
# laser cycles, with the cloud of all 2000 electrons of the run. The paths and
# the snapshots are the archived Grace sets.

using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

function make_traj(name::String, path::String, snapshot::String,
                   paths::Vector{String}; colors = [:red, :blue])
    g = read_xmgr(path)[1]
    getset(n) = g.sets[findfirst(s -> s.name == n, g.sets)].data
    cloud = getset(snapshot)
    fig = Figure(size = (620, 620))
    ax = Axis(fig[1, 1], aspect = DataAspect(), limits = (-14.6, 14.6, -14.6, 14.6),
              xlabel = "x", ylabel = "y")
    scatter!(ax, cloud[:, 1], cloud[:, 2], color = (:gray, 0.4), markersize = 3)
    for (i, n) in enumerate(paths)
        d = getset(n)
        lines!(ax, d[:, 1], d[:, 2], color = colors[i], linewidth = 1.2,
               label = "trajectory $i")
    end
    length(paths) > 1 && axislegend(ax; position = :lt)
    save_fig(fig, name)
    return fig
end

make_traj("traj79",
          figure_data("traj79.xmgr.gz"), "S2", ["S0"])
make_traj("traj157",
          figure_data("traj157.xmgr.gz"), "S1", ["S0", "S2"])
