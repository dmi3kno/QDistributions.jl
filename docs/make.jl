# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, QDistributions

makedocs(
    modules = [QDistributions],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "dmi3kno",
    sitename = "QDistributions.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/dmi3kno/QDistributions.jl.git",
    push_preview = true
)
