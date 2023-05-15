push!(LOAD_PATH,"../src/")
using PlmDCA
using Documenter

makedocs(;
    modules=[PlmDCA],
    authors="Andrea Pagnani",
    clean=true,
    #repo="https://github.com/pagnani/ArDCA.jl/blob/{commit}{path}#L{line}",
    sitename="PlmDCA",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pagnani.github.io/PlmDCA",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pagnani/PlmDCA.git",
)
