"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData
using MLStyle: @match

using QuantumESPRESSOBase.Namelists: Namelist, to_dict
using QuantumESPRESSOBase.Namelists.PWscf

using QuantumESPRESSOParsers
using QuantumESPRESSOParsers.InputParsers.Namelists

export findnamelists, parsenamelists

function findnamelists(str::AbstractString)
    captured = map(x -> x.captures, eachmatch(QuantumESPRESSOParsers.NAMELIST_BLOCK_REGEX, str))
    dict = Dict{Symbol, String}()
    for (name, content) in captured
        @match uppercase(name) begin
            "CONTROL" => push!(dict, :ControlNamelist => content)
            "SYSTEM" => push!(dict, :SystemNamelist => content)
            "ELECTRONS" => push!(dict, :ElectronsNamelist => content)
            "CELL" => push!(dict, :CellNamelist => content)
            "IONS" => push!(dict, :IonsNamelist => content)
        end
    end
    return dict
end # function findnamelists

function parsenamelists(str::AbstractString)
    dict = findnamelists(str)
    return [parse(eval(key), value) for (key, value) in dict]
end # function parsenamelists

end
