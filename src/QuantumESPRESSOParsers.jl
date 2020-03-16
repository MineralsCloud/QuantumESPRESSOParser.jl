module QuantumESPRESSOParsers

# Referenced from https://discourse.julialang.org/t/how-to-get-the-non-nothing-type-from-union-t-nothing/30523
nonnothingtype(::Type{T}) where {T} = Core.Compiler.typesubtract(T, Nothing)  # Should not be exported

include("Inputs/Inputs.jl")
include("Outputs/Outputs.jl")

end # module
