# Redraw the two historical individual-step diagrams from the bundled Grace data.

using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

function main()
    g = read_xmgr(figure_data("indiv.agr.gz"))[1]
    fig = Figure(size = (700, 520))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(),
              title = "Trajectories and individual steps (2000 archive)")
    colors = (:royalblue, :orangered, :seagreen)
    for (i, s) in enumerate(g.sets[1:3])
        lines!(ax, s.data[:, 1], s.data[:, 2], color = colors[i], linewidth = 1.4,
               label = "particle $i")
        scatter!(ax, s.data[:, 1], s.data[:, 2], color = colors[i], markersize = 3)
    end
    axislegend(ax; position = :rt)
    save_fig(fig, "indiv")

    gs = read_xmgr(figure_data("indiv2.agr.gz"))
    fig2 = Figure(size = (700, 520))
    ax2 = Axis(fig2[1, 1], xlabel = "step number", ylabel = "time",
               title = "Computed times by particle (2000 archive)")
    for (i, s) in enumerate(gs[1].sets[2:4])
        lines!(ax2, s.data[:, 1], s.data[:, 2], color = colors[i], linewidth = 1.5,
               label = "particle $i")
        scatter!(ax2, s.data[:, 1], s.data[:, 2], color = colors[i], markersize = 3)
    end
    axislegend(ax2; position = :rb)
    save_fig(fig2, "indiv2")
    return fig, fig2
end

main()
