using QuantumESPRESSOParsers
using Test

@testset "QuantumESPRESSOParsers.jl" begin
    # Write your own tests here.
    include("Inputs/Namelists.jl")
    include("Outputs/PWscf.jl")
end
