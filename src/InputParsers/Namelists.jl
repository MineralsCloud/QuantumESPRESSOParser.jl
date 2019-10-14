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
using QuantumESPRESSOBase.Namelists.CP
using QuantumESPRESSOBase.Namelists.PHonon

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py.
const NAMELIST_BLOCK =
r"""
^ [ \t]* &(\S+) [ \t]* $\n  # match line w/ nmlst tag; save nmlst name
(
    [\S\s]*?                # match any line non-greedily
)                           # save the group of text between nmlst
^ [ \t]* / [ \t]* $\n       # match line w/ "/" as only non-whitespace char
"""mx
# This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py.
const NAMELIST_ITEM =
r"""
[ \t]* (?<key> \S+? )(?: (?<kind> [\(%]) (?<index> \w+) \)? )? [ \t]*  # match and store key
=                              # equals sign separates key and value
[ \t]* (?<value> \S+?) [ \t]*  # match and store value
[\n,]                          # return or comma separates "key = value" pairs
"""mx
const NAMELIST_HEADS = Dict{Any,String}(
    PWscf.ControlNamelist => "CONTROL",
    PWscf.SystemNamelist => "SYSTEM",
    PWscf.ElectronsNamelist => "ELECTRONS",
    PWscf.CellNamelist => "CELL",
    PWscf.IonsNamelist => "IONS",
    CP.ControlNamelist => "CONTROL",
    CP.SystemNamelist => "SYSTEM",
    CP.ElectronsNamelist => "ELECTRONS",
    CP.CellNamelist => "CELL",
    CP.IonsNamelist => "IONS",
    WannierNamelist => "WANNIER",
    PHNamelist => "INPUTPH",
    Q2RNamelist => "INPUT",
    MatdynNamelist => "INPUT",
    DynmatNamelist => "INPUT",
)

function Base.parse(T::Type{<:Namelist}, str::AbstractString)
    result = Dict{Symbol,Any}()
    for nml in eachmatch(NAMELIST_BLOCK, str)
        head, body = nml.captures
        !occursin(uppercase(head), NAMELIST_HEADS[T]) && continue
        for m in eachmatch(NAMELIST_ITEM, body)
            k = Symbol(m[:key])
            v = FortranData(string(m[:value]))
            # Parse a `FortranData` from `value` as type of the field of the namelist `T`
            if isnothing(m[:index])  # Cases like `ntyp = 2`
                result[k] = parse(fieldtype(T, k), v)
            else  # An entry with multiple values, e.g., `celldm(2) = 3.0`.
                if m[:kind] == '('
                    i = parse(Int, m[:index])
                    v = parse(Float64, v)  # TODO: This is tricky.
                    result[k] = if haskey(result, k)
                        # If `celldm` occurs before, push the new value, else create a vector of pairs.
                        fillbyindex!(result[k], i, v)
                    else
                        fillbyindex!([], i, v)
                    end
                else  # m[:kind] == '%'
                    i = string(m[:index])
                    # TODO: This is not finished!
                end
            end
        end
    end
    return isempty(result) ? nothing : T(T(), result)  # TODO: This does not dynamically change
end # function Base.parse

function fillbyindex!(x::AbstractVector, index::Int, value::T) where {T}
    if isempty(x)
        x = Vector{Union{Nothing,T}}(nothing, index)
    else
        index > length(x) && append!(
            x,
            Vector{Union{Nothing,T}}(nothing, index - length(x)),
        )
    end
    x[index] = value
    return x
end # function fillbyindex!

end
