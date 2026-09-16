# fig5b: collision frequency ν_ei/ωp vs v0/vth at ω0 = 3 ωp (Γ = 0.1, Z = 1).
# The curves (Silin model, Silin model with corrected Coulomb log, Dawson and
# Oberman model) and the simulation points are the ones of the archived Grace
# project, redrawn with their recorded line styles and symbols.

using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

function main()
    g = read_xmgr(figure_data("fig5b.xmgr.gz"))[1]
    fig = Figure(size = (700, 560))
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10,
              limits = (0.1, 10.0, 1e-4, 0.2), xlabel = "v₀/v_th",
              ylabel = "ν_ei/ω_p",
              title = "Collision frequency\nΓ = 0.1 ; Z = 1 ; ω₀ = 3 ω_p")
    curves = [("S0", "Silin model", :dash),
              ("S1", "Silin model ; corrected Coulomb log", :solid),
              ("S2", "Dawson and Oberman model", :dot)]
    for (name, label, ls) in curves
        s = g.sets[findfirst(x -> x.name == name, g.sets)]
        lines!(ax, s.data[:, 1], s.data[:, 2], color = :black, linestyle = ls,
               linewidth = 2, label = label)
    end
    points = [("S3", "Simulation ; Si potential", :green, :utriangle),
              ("S4", "Simulation ; b_min = 0.144", :red, :rect),
              ("S6", "Simulation ; b_min = 0.05", :blue, :circle)]
    for (name, label, color, marker) in points
        s = g.sets[findfirst(x -> x.name == name, g.sets)]
        errorbars!(ax, s.data[:, 1], s.data[:, 2], s.data[:, 3], color = color)
        scatter!(ax, s.data[:, 1], s.data[:, 2], color = color, marker = marker,
                 markersize = 10, label = label)
    end
    axislegend(ax; position = :lb)
    save_fig(fig, "fig5b")
    return fig
end

main()
