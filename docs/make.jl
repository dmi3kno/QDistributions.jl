using Documenter
using QDistributions

makedocs(
    sitename = "QDistributions",
    format = Documenter.HTML(),
    modules = [QDistributions],
    checkdocs=:exports
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
