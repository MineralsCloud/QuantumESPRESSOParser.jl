module PWscf

using Test

using QuantumESPRESSOBase.Inputs.PWscf
using QuantumESPRESSOParser.Inputs

@testset "Parse empty strings" begin
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
    @test tryparse(ControlNamelist, " ") === nothing
    @test tryparse(ControlNamelist, "&control/") == ControlNamelist()
    @test parse(ControlNamelist, "&control\n/\n") == ControlNamelist()
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
end

end
