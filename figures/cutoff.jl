# cutoff: the three regularizations of the Coulomb potential compared —
# bare 1/r, Fourier-cutoff (2/π) Si(k_max r)/r, soft core 1/√(r²+ε²).

using LaserPlasma

const PC = LaserPlasma
include(joinpath(@__DIR__, "common.jl"))

function main()
    r = range(0.02, 2.0, length = 500)
    kmax = 10.0
    eps = 0.1
    coulomb = @. 1.0 / r
    cutoff = @. (2.0 / π) * PC.si(kmax * r) / r
    softcore = @. 1.0 / sqrt(r^2 + eps^2)

    fig = Figure(size = (620, 540))
    ax = Axis(fig[1, 1], xlabel = "r", ylabel = "potential",
              title = "Regularized Coulomb potentials")
    ylims!(ax, 0.0, 12.0)
    lines!(ax, r, coulomb, color = :black, linestyle = :dot, linewidth = 2,
           label = "1.0/r")
    lines!(ax, r, cutoff, color = :red, linewidth = 2,
           label = "(2.0/π) Si(k_max r)/r ; k_max = 10")
    lines!(ax, r, softcore, color = :royalblue, linestyle = :dash, linewidth = 2,
           label = "1/√(r²+ε²) ; ε = 0.1")
    axislegend(ax; position = :rt)

    save_fig(fig, "cutoff")
    return fig
end

main()
