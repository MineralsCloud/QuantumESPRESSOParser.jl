using QuantumESPRESSOParser
using Documenter

makedocs(;
    modules=[QuantumESPRESSOParser],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOParser.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOParser.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOParser.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOParser.jl",
)
