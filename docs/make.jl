using QuantumESPRESSOParser
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOParser, :DocTestSetup, :(using QuantumESPRESSOParser); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOParser],
    authors="Reno <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOParser.jl/blob/{commit}{path}#{line}",
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
