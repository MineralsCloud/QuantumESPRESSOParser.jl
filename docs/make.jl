using QuantumESPRESSOParsers
using Documenter

makedocs(;
    modules=[QuantumESPRESSOParsers],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOParsers.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOParsers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOParsers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOParsers.jl",
)
