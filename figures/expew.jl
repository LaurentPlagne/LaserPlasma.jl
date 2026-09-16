# expew3 / expew5: collision frequency ν_ei/ωp vs v0/vth for Γ = 0.1, Z = 1,
# ω0 = 3 and 5 ωp. Lines are the dielectric theory of `dielectric.jl` (Decker
# model, soft-core potentials ε = 0.05 and 0.144); markers with error bars are
# the archived March 2000 simulation points (N = 500 runs).

using LaserPlasma
using CairoMakie

const PC = LaserPlasma
include(joinpath(@__DIR__, "common.jl"))

const KT = 6.935
const VTH = sqrt(KT)

function theory_curves(omega; npts = 40)
    rapv = exp.(range(log(0.01), log(100.0), length = npts))
    soft(bmin) = [PC.decker_coll_freq(r * VTH, PC.DielectricTheory(KT, omega, 1.0, bmin);
                                      softcore = true) for r in rapv]
    bare = [PC.decker_coll_freq(r * VTH, PC.DielectricTheory(KT, omega, 1.0, 0.05))
            for r in rapv]
    return rapv, bare, soft(0.05), soft(0.144)
end

function make_figure(omega::Float64, name::String, simsets::Tuple{String,String})
    rapv, bare, soft05, soft144 = theory_curves(omega)
    raw = read_xmgr(figure_data(name * ".xmgr.gz"))[1]
    fig = Figure(size = (700, 560))
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10,
              limits = (0.1, 10.0, 1e-4, 0.1), xlabel = "v₀/v_th",
              ylabel = "ν_ei/ω_p",
              title = "Γ = 0.1 ; Z = 1 ; ω₀ = $omega ω_p")
    lines!(ax, rapv, bare, color = :black, linewidth = 1.5, label = "Decker model")
    lines!(ax, rapv, soft05, color = :orangered, linewidth = 1.5,
           label = "Soft-core potential ε = 0.05")
    lines!(ax, rapv, soft144, color = :royalblue, linewidth = 1.5,
           label = "Soft-core potential ε = 0.144")
    for (setname, bmin, color, marker) in zip(simsets, (0.05, 0.144),
                                              (:orangered, :royalblue), (:circle, :rect))
        s = raw.sets[findfirst(x -> x.name == setname, raw.sets)]
        errorbars!(ax, s.data[:, 1], s.data[:, 2], s.data[:, 3], color = color)
        scatter!(ax, s.data[:, 1], s.data[:, 2], color = color, marker = marker,
                 markersize = 9, label = "Simulation ; ε = $bmin")
    end
    axislegend(ax; position = :lb)
    save_fig(fig, name)
    return fig
end

make_figure(3.0, "expew3", ("S8", "S9"))
make_figure(5.0, "expew5", ("S6", "S7"))
