"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using Compat: isnothing
using PyFortran90Namelists: FortranData
using QuantumESPRESSOBase.Inputs: Namelist, titleof

using QuantumESPRESSOParsers: nonnothingtype

struct InvalidUserInput
    msg::String
end

# From https://github.com/aiidateam/qe-tools/blob/570a648/qe_tools/parsers/qeinputparser.py#L315-L321
const NAMELIST_ITEM = r"""
                      [ \t]* (?<key>\w+?) (?: (?<kind>[\(%]) (?<index>\w+) \)? )? [ \t]*  # Match and store key
                      =                              # Equals sign separates key and value
                      [ \t]* (?<value>\S+?) [ \t]*  # Match and store value
                      [\n,]                          # Return or comma separates "key = value" pairs
                      """mx

# This is an internal function and should not be exported.
function tryparse_internal(::Type{T}, str::AbstractString, raise::Bool) where {T<:Namelist}
    result = Dict{Symbol,Any}()
    head = titleof(T)
    # From https://github.com/aiidateam/qe-tools/blob/570a648/qe_tools/parsers/qeinputparser.py#L305-L312
    NAMELIST_BLOCK = Regex(
        """
        ^ [ \\t]* &$head [ \\t]* \$  # Match `Namelist`'s name
        (?<body>
         [\\S\\s]*?  # Match any line non-greedily
        )            # Save the group of text between `Namelist`s
        ^ [ \\t]* \\/ [ \\t]* \$  # Match line with "/" as the only non-whitespace char
        """,
        "imx",
    )
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
                S = nonnothingtype(eltype(fieldtype(T, k)))
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

Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist} =
    tryparse_internal(T, str, false)
Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist} =
    tryparse_internal(T, str, true)

include("PWscf.jl")
# include("PHonon.jl")

end
