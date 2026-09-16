# Diel.w10: collision frequency from the dielectric theory, Γ=0.1, Z=1,
# ω0 = 10 ωp. Curves: Decker model (bare Coulomb with the k cutoff) and
# soft-core potentials ε = 0.05 and 0.144.

using LaserPlasma

const PC = LaserPlasma
include(joinpath(@__DIR__, "common.jl"))

function dielectric_curves(omega; npts = 60, rapv_lo = 0.01, rapv_hi = 100.0)
    kT = 6.935
    vth = sqrt(kT)
    rapv = exp.(range(log(rapv_lo), log(rapv_hi), length = npts))
    curve(softcore, bmin) = begin
        th = PC.DielectricTheory(kT, omega, 1.0, bmin)
        [PC.decker_coll_freq(r * vth, th; softcore = softcore) for r in rapv]
    end
    return rapv, curve(false, 0.05), curve(true, 0.05), curve(true, 0.144)
end

function main()
    rapv, decker, soft005, soft144 = dielectric_curves(10.0)
    fig = Figure(size = (760, 620))
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10,
              xlabel = "v₀/v_th", ylabel = "ν_ei/ω_p",
              title = "Dielectric Theory\nΓ = 0.1 ; Z = 1 ; ω₀ = 10 ω_p")
    lines!(ax, rapv, decker, color = :black, linewidth = 2, label = "Decker model")
    lines!(ax, rapv, soft005, color = :orangered, linewidth = 2,
           label = "Soft-core potential ε = 0.05")
    lines!(ax, rapv, soft144, color = :royalblue, linewidth = 2,
           label = "Soft-core potential ε = 0.144")
    axislegend(ax; position = :rb)
    save_fig(fig, "diel_w10")
    return fig
end

main()
