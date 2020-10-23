using PanelLag
using Documenter

makedocs(;
    modules=[PanelLag],
    authors="Jacob Adenbaum",
    repo="https://github.com/jacobadenbaum/PanelLag.jl/blob/{commit}{path}#L{line}",
    sitename="PanelLag.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jacobadenbaum.github.io/PanelLag.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jacobadenbaum/PanelLag.jl",
)
