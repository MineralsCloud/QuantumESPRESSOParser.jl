"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData
using MLStyle: @match

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const NAMELIST_BLOCK_REGEX = r"""
^ [ \t]* &(\S+) [ \t]* $\n  # match line w/ nmlst tag; save nmlst name
(
 [\S\s]*?                # match any line non-greedily
)                        # save the group of text between nmlst
^ [ \t]* / [ \t]* $\n    # match line w/ "/" as only non-whitespace char
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const NAMELIST_ITEM_REGEX = r"""
[ \t]* (\S+?) [ \t]*  # match and store key
=               # equals sign separates key and value
[ \t]* (\S+?) [ \t]*  # match and store value
[\n,]           # return or comma separates "key = value" pairs
"""mx

function findnamelists(str::AbstractString)
    captured = map(x -> x.captures, eachmatch(NAMELIST_BLOCK_REGEX, str))
    dict = Dict{Symbol,String}()
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

function lexnamelist(content::AbstractString)
    captured = map(x -> x.captures, eachmatch(NAMELIST_ITEM_REGEX, content))
    dict = Dict{String,Any}()
    for (key, value) in captured
        dict[key] = value
    end
    return dict
end # function lexnamelist

function Base.parse(T::Type{Namelist}, str::AbstractString)
    dict = findnamelists(str)
    return [parse(eval(key), value) for (key, value) in dict]
end # function Base.parse
function Base.parse(T::Type{<:Namelist}, content::AbstractString)
    regex = r"([\w\d]+)(?:\((\d+)\))?"
    dict = lexnamelist(content)
    result = Dict{Symbol,Any}()
    for (key, value) in dict
        captures = match(regex, key).captures
        k = Symbol(string(captures[1]))
        v = FortranData(string(value))
        # We need to parse a `FortranData` from `value` as type of the field of the namelist `T`.
        if !isnothing(captures[2])  # An entry with multiple values, e.g., `celldm(2) = 3.0`.
            # If `celldm` occurs before, push the new value, else create a vector of pairs.
            i = parse(Int, captures[2])
            v = parse(typeintersect(eltype(fieldtype(T, k)), Union{Int,Float64}), v)
            result[k] = (haskey(result, k) ? fillbyindex!(result[k], i, v) : fillbyindex!([], i, v))
        else  # Cases like `ntyp = 2`
            result[k] = parse(fieldtype(T, k), v)
        end
    end
    return T(T(), result)
end # function parsenamelist

function fillbyindex!(x::AbstractVector, index::Int, value::T) where {T}
    if isempty(x)
        x = Vector{Union{Missing,T}}(missing, index)
    else
        index > length(x) && append!(x, Vector{Union{Missing,T}}(missing, index - length(x)))
    end
    x[index] = value
    return x
end # function fillbyindex!

end
