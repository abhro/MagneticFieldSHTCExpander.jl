using Documenter
using MagneticFieldSHTCExpander
using Documenter.Remotes: GitHub

makedocs(
    modules = [MagneticFieldSHTCExpander],
    sitename = "MagneticFieldSHTCExpander",
    repo = GitHub("abhro/MagneticFieldSHTCExpander.jl"),
)

deploydocs(
    repo="github.com/abhro/MagneticFieldSHTCExpander.jl.git",
    push_preview=true,
)
