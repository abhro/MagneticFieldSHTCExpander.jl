using Documenter
using Documenter.Remotes: GitHub
using DocumenterCitations
using MagneticFieldSHTCExpander

include("pages.jl")

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :authoryear,
)

mathengine = Documenter.KaTeX(
    Dict(
        :macros => Dict(
            raw"\pdv"  => raw"\frac{∂#1}{∂#2}",
            raw"\dv"   => raw"\frac{d#1}{d#2}",
            raw"\dpdv" => raw"\dfrac{∂#1}{∂#2}",
            raw"\ddv"  => raw"\dfrac{d#1}{d#2}",
            raw"\TT"   => raw"\mathsf{T}",
            raw"\grad" => raw"\boldsymbol{∇}",
        )
    )
)

makedocs(
    modules = [MagneticFieldSHTCExpander],
    sitename = "MagneticFieldSHTCExpander",
    repo = GitHub("abhro/MagneticFieldSHTCExpander.jl"),
    pages = pages,
    plugins = [bib],
    format = Documenter.HTML(; mathengine)
)

deploydocs(
    repo = "github.com/abhro/MagneticFieldSHTCExpander.jl.git",
    push_preview = true,
)
