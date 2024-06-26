using Documenter
using Documenter.Remotes: GitHub
using MagneticFieldSHTCExpander

include("pages.jl")

makedocs(
    modules = [MagneticFieldSHTCExpander],
    sitename = "MagneticFieldSHTCExpander",
    repo = GitHub("abhro/MagneticFieldSHTCExpander.jl"),
    pages = pages,
)

deploydocs(
    repo = "github.com/abhro/MagneticFieldSHTCExpander.jl.git",
    push_preview = true,
)
