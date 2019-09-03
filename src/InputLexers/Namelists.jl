"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using QuantumESPRESSOParsers

export lexnamelist

function lexnamelist(content::AbstractString)
    captured = map(x -> x.captures, eachmatch(QuantumESPRESSOParsers.NAMELIST_ITEM_REGEX, content))
    dict = Dict{String, Any}()
    for (key, value) in captured
        dict[key] = value
    end
    return dict
end # function lexnamelist

end
