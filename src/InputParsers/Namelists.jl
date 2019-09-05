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

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf

export lexnamelist

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const NAMELIST_ITEM_REGEX = r"""
[ \t]* (\S+?) [ \t]*  # match and store key
=               # equals sign separates key and value
[ \t]* (\S+?) [ \t]*  # match and store value
[\n,]           # return or comma separates "key = value" pairs
"""mx

function lexnamelist(content::AbstractString)
    captured = map(x -> x.captures, eachmatch(NAMELIST_ITEM_REGEX, content))
    dict = Dict{String, Any}()
    for (key, value) in captured
        dict[key] = value
    end
    return dict
end # function lexnamelist

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
            v = parse(typeintersect(eltype(fieldtype(T, k)), Union{Int, Float64}), v)
            result[k] = (haskey(result, k) ? fillbyindex!(result[k], i, v) : fillbyindex!([], i, v))
        else  # Cases like `ntyp = 2`
            result[k] = parse(fieldtype(T, k), v)
        end
    end
    final = merge(to_dict(T()), result)
    return T((final[f] for f in fieldnames(T))...)
end # function parsenamelist

function fillbyindex!(x::AbstractVector, index::Int, value::T) where {T}
    if isempty(x)
        x = Vector{Union{Missing, T}}(missing, index)
    else
        index > length(x) && append!(x, Vector{Union{Missing, T}}(missing, index - length(x)))
    end
    x[index] = value
    return x
end # function fillbyindex!

end
