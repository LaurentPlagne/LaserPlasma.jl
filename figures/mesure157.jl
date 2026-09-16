# mesure157: total energy E_tot(t) and absorbed energy ΣΔU_i(t) for the
# v0/vth = 0.2 run (N = 2000, ε = 0.144, F0 = 1.57, ω0 = 3 ωp, thermostat on).
# The curves are the archived Grace sets of the TREE_RL run (the power curve is
# `power.dat` column 11, the cumulative one is deltaU normalised by the archive
# plotting script).

using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

function main()
    graphs = read_xmgr(figure_data("mesure157.xmgr.gz"))
    etot = graphs[1].sets[1]
    power = graphs[2].sets[1]
    dU = graphs[2].sets[2]
    fig = Figure(size = (760, 700))
    ax1 = Axis(fig[1, 1], xlabel = "time", ylabel = "E_tot",
               title = "v₀/v_th = 0.2 ; N = 2000 ; ε = 0.144 ; ω₀ = 3 ω_p")
    lines!(ax1, etot.data[:, 1], etot.data[:, 2], color = :black, linewidth = 1.5)
    ax2 = Axis(fig[2, 1], xlabel = "time", ylabel = "ΣΔU_i")
    ylims!(ax2, -0.04, 0.04)
    lines!(ax2, power.data[:, 1], power.data[:, 2], color = :black, linewidth = 1,
           label = "P(t)")
    lines!(ax2, dU.data[:, 1], dU.data[:, 2], color = :red, linewidth = 1.5,
           label = "ΣΔU_i")
    axislegend(ax2; position = :rb)
    linkxaxes!(ax1, ax2)
    save_fig(fig, "mesure157")
    return fig
end

main()
