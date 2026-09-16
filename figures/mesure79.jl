# mesure79: absorbed power P(t) (column 11 of power.dat, `power1`) and
# cumulative absorbed energy ΣΔU_i for the v0/vth = 10 run
# (N=2000, ε=0.05, F0=79, ω=3). The archived data are bundled in figures/data;
# optionally overlay a Julia run by setting LASERPLASMA_RUN to its output folder.

include(joinpath(@__DIR__, "common.jl"))

function readcols(path)
    rows = [permutedims(parse.(Float64, split(l)))
            for l in read_lines(path) if !isempty(strip(l))]
    return reduce(vcat, rows)
end

function main()
    power = readcols(figure_data("mesure79_power.dat.gz"))
    du = readcols(figure_data("mesure79_deltaU.gz"))
    run_dir = get(ENV, "LASERPLASMA_RUN", "")
    has_run = !isempty(run_dir)
    run_power = has_run ? readcols(joinpath(run_dir, "power.dat")) : nothing
    run_du = has_run ? readcols(joinpath(run_dir, "deltaU")) : nothing

    fig = Figure(size = (760, 760))
    ax1 = Axis(fig[1, 1], xlabel = "time", ylabel = "P(t)",
               title = "v₀/v_th = 10 ; N = 2000 ; ε = 0.05 ; ω₀ = 3 ω_p")
    lines!(ax1, power[:, 1], power[:, 11], color = :black, linewidth = 1,
           label = "archive 2000 (tree)")
    if has_run
        lines!(ax1, run_power[:, 1], run_power[:, 11], color = :red, linewidth = 1,
               alpha = 0.6, label = "Julia run")
    end
    axislegend(ax1; position = :rb)

    ax2 = Axis(fig[2, 1], xlabel = "time", ylabel = "ΣΔUᵢ − ΣΔU₀",
               title = "absorbed energy per cycle")
    scatterlines!(ax2, du[:, 1], du[:, 3] .- du[1, 3], color = :black, marker = :circle,
              markersize = 8, label = "archive 2000 (tree)")
    if has_run
        scatterlines!(ax2, run_du[:, 1], run_du[:, 3] .- run_du[1, 3], color = :red,
              marker = :xcross, markersize = 8, label = "Julia run")
    end
    axislegend(ax2; position = :lt)
    linkxaxes!(ax1, ax2)

    save_fig(fig, "mesure79")
    return fig
end

main()
