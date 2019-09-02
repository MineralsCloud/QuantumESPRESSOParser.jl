"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using MLStyle: @match

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

end
