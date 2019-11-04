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

using QuantumESPRESSOParsers

# This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py.
const NAMELIST_ITEM = r"""
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
    head = NAMELIST_HEADS[T]
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
        @info("Namelist not found in string!")
        return
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
                i < 0 || throw(InvalidUserInput("Negative index found in $(item[:index])!"))
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
    if isempty(result)
        @info("Namelist found, but it is empty! Default values will be used!")
    end
    # Works even if `result` is empty. If empty, it just returns the default `Namelist`.
    return T(; result...)
end # function Base.parse

end
