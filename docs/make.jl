using Documenter, QuantumESPRESSOParsers

makedocs(;
    modules=[QuantumESPRESSOParsers],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MineralsCloud/QuantumESPRESSOParsers.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOParsers.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOParsers.jl",
)
