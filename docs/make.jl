using Documenter
using Documenter.Remotes: GitHub

# The guide is hand-written Markdown: it pulls in no docstrings, so the package
# itself is not a dependency of this environment and the docs build needs neither
# CairoMakie nor GLMakie. The figures it displays are committed under
# `src/assets/figures/` and refreshed by `figures/build_all.jl`.

makedocs(
    sitename = "LaserPlasma.jl",
    authors = "Laurent Plagne",
    # Declared explicitly rather than sniffed from `git remote`, so the build
    # works in a fresh clone, in CI, and before any remote exists.
    repo = GitHub("LaurentPlagne", "LaserPlasma.jl"),
    # Flat file names locally so `docs/build/index.html` opens directly; pretty
    # directory URLs when deployed.
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://laurentplagne.github.io/LaserPlasma.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Welcome" => "index.md",
        "Getting started" => "getting-started.md",
        "The model" => "model.md",
        "Algorithms" => "methods.md",
        "Figure gallery" => "gallery.md",
        "The method in pictures" => "method-illustrations.md",
        "Code guide" => "code-guide.md",
        "Glossary" => "glossary.md",
    ],
)

deploydocs(;
    repo = "github.com/LaurentPlagne/LaserPlasma.jl.git",
    devbranch = "main",
    push_preview = true,
)
