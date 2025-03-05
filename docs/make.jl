using Documenter
using Documenter.Remotes: GitHub
using DocumenterCitations
using MagneticFieldSHTCExpander

include("pages.jl")

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :authoryear,
)

makedocs(
    modules = [MagneticFieldSHTCExpander],
    sitename = "MagneticFieldSHTCExpander",
    repo = GitHub("abhro/MagneticFieldSHTCExpander.jl"),
    pages = pages,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/abhro/MagneticFieldSHTCExpander.jl.git",
    push_preview = true,
)
