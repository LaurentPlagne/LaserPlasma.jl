# fleur: zoom on one electron trajectory (the "flower") with the neighbouring
# particles of the r_ws = 1 DIRE run, from the archived Grace project.

using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

function main()
    g = read_xmgr(figure_data("fleur.xmgr.gz"))[1]
    curve = g.sets[findfirst(s -> s.name == "S0", g.sets)].data
    cloud = g.sets[findfirst(s -> s.name == "S1", g.sets)].data
    fig = Figure(size = (620, 560))
    ax = Axis(fig[1, 1], aspect = DataAspect(),
              limits = (-2.9, -2.3, -1.85, -1.6), xlabel = "x", ylabel = "y")
    scatter!(ax, cloud[:, 1], cloud[:, 2], color = (:gray, 0.5), markersize = 4)
    lines!(ax, curve[:, 1], curve[:, 2], color = :red, linewidth = 1)
    save_fig(fig, "fleur")
    return fig
end

main()
