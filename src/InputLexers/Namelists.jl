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

end
