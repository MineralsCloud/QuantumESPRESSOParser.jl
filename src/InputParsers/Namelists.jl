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

using QuantumESPRESSOParsers

# This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py.
const NAMELIST_ITEM = r"""
                      [ \t]* (?<key> \S+? )(?: (?<kind> [\(%]) (?<index> \w+) \)? )? [ \t]*  # match and store key
                      =                              # equals sign separates key and value
                      [ \t]* (?<value> \S+?) [ \t]*  # match and store value
                      [\n,]                          # return or comma separates "key = value" pairs
                      """mx
const NAMELIST_HEADS = Dict{Symbol,String}(
    :ControlNamelist => "CONTROL",
    :SystemNamelist => "SYSTEM",
    :ElectronsNamelist => "ELECTRONS",
    :CellNamelist => "CELL",
    :IonsNamelist => "IONS",
    :WannierNamelist => "WANNIER",
    :PHNamelist => "INPUTPH",
    :Q2RNamelist => "INPUT",
    :MatdynNamelist => "INPUT",
    :DynmatNamelist => "INPUT",
)

# This is an internal function and should not be exported.
function tryparse_internal(::Type{T}, str::AbstractString, raise::Bool) where {T<:Namelist}
    result = Dict{Symbol,Any}()
    head = NAMELIST_HEADS[nameof(T)]
    # This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py.
    NAMELIST_BLOCK = Regex("""
                           ^ [ \t]* &$head [ \t]* \$\n
                           (?<body>
                               [\\S\\s]*?
                           )
                           ^ [ \t]* / [ \t]* \$
                           """, "imx")
    m = match(NAMELIST_BLOCK, str)
    if isnothing(m)
        raise ? throw(Meta.ParseError("Namelist not found in string!")) : return
    end
    for item in eachmatch(NAMELIST_ITEM, m[:body])
        k = Symbol(item[:key])
        v = FortranData(string(item[:value]))
        # Parse a `FortranData` from `value` as type of the field of the namelist `T`
        if isnothing(item[:index])  # Cases like `ntyp = 2`
            result[k] = parse(fieldtype(T, k), v)
        else  # An entry with multiple values, e.g., `celldm(2) = 3.0`.
            if item[:kind] == "("  # Note: it cannot be `'('`. It will result in `false`!
                i = parse(Int, item[:index])
                i < 0 && throw(InvalidUserInput("Negative index found in $(item[:index])!"))
                S = QuantumESPRESSOParsers.nonnothingtype(eltype(fieldtype(T, k)))
                v = parse(S, v)
                arr = get(result, k, [])
                if i > length(arr)  # Works even if `x` is empty. If empty, `length(x)` will be `0`.
                    append!(arr, Vector{Union{Nothing,T}}(nothing, i - length(arr)))
                end
                arr[i] = v  # Now `index` cannot be greater than `length(x)` => a normal assignment
                result[k] = arr
            else  # item[:kind] == '%'
                i = string(item[:index])
                # TODO: This is not finished!
            end
        end
    end
    isempty(result) && @info("An empty Namelist found! Default values will be used!")
    # Works even if `result` is empty. If empty, it just returns the default `Namelist`.
    return T(; result...)
end # function tryparse_internal

function Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist}
    return tryparse_internal(T, str, false)
end # function Base.tryparse
function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    return tryparse_internal(T, str, true)
end # function Base.parse

end
