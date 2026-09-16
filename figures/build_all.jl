# Rebuild every figure cited by the historical plasma slides, then refresh the
# images displayed by Documenter. Run: julia --project=. figures/build_all.jl
# The scientific plots use bundled data or the ported dielectric calculation;
# no N = 2000 simulation is started by this script.

const FIGURE_SCRIPTS = (
    "diagrams.jl", "cutoff.jl", "diel_w10.jl", "traj.jl",
    "mesure157.jl", "mesure79.jl", "fig5b.jl", "expew.jl",
    "fleur.jl", "indiv.jl", "method3d.jl",
)

const RESULT_FIGURES = (
    "cutoff", "diel_w10", "traj79", "traj157", "mesure157",
    "mesure79", "fig5b", "expew5", "expew3", "fleur",
    "indiv", "indiv2", "method3d",
)

const METHOD_FIGURES = (
    "intro", "field", "periodic_1", "periodic_2", "periodic_3",
    "ewald_1", "ewald_2", "tree_1", "tree_2",
    "traject_1", "traject_2", "traject_3", "traject_4",
    "traject_step1", "traject_step2", "traject_step3",
)

function run_figure_script(script)
    # Each historical script defines its own `main`, helper functions and
    # constants. A fresh module keeps their names from colliding.
    scope = Module(gensym(:HistoricalFigure))
    Core.eval(scope, :(include(path::AbstractString) = Base.include(@__MODULE__, path)))
    withenv("LASERPLASMA_RUN" => nothing) do
        Base.include(scope, joinpath(@__DIR__, script))
    end
end

function copy_figure(name, extension, asset_dir)
    filename = "$name.$extension"
    source = joinpath(@__DIR__, "out", filename)
    isfile(source) || error("Expected figure missing: $source")
    cp(source, joinpath(asset_dir, filename); force = true)
end

function check_slide_coverage()
    source = read(`gzip -dc $(joinpath(@__DIR__, "data", "plasma.tex.gz"))`, String)
    cited = Set(basename(match.captures[1]) for match in
                eachmatch(r"\{([^{}]+)\.eps\}", source))
    generated = Set((RESULT_FIGURES..., METHOD_FIGURES...))
    missing = setdiff(Set(replace(name, "Diel.w10" => "diel_w10") for name in cited),
                      generated)
    isempty(missing) || error("No Julia figure for slide sources: $(sort!(collect(missing)))")
    return length(cited)
end

function main()
    cited = check_slide_coverage()
    for script in FIGURE_SCRIPTS
        @info "Building historical figures" script
        run_figure_script(script)
    end
    asset_dir = joinpath(@__DIR__, "..", "docs", "src", "assets", "figures")
    mkpath(asset_dir)
    for name in RESULT_FIGURES
        copy_figure(name, "png", asset_dir)
        copy_figure(name, "pdf", asset_dir)
    end
    for name in METHOD_FIGURES
        copy_figure(name, "svg", asset_dir)
    end
    @info "Documenter gallery refreshed" slide_figures = cited plots = length(RESULT_FIGURES) diagrams = length(METHOD_FIGURES)
    return nothing
end

main()
