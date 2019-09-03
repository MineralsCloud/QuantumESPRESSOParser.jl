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

export findnamelists, parsenamelists, findcards

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const ATOMIC_POSITIONS_BLOCK_REGEX = r"""
^ \s* ATOMIC_POSITIONS \s*                      # Atomic positions start with that string
[{(]? \s* (?P<units>\S+?)? \s* [)}]? \s* $\n    # The units are after the string in optional brackets
(?P<block>                                      # This is the block of positions
    (
        (
            \s*                                 # White space in front of the element spec is ok
            (
                [A-Za-z]+[A-Za-z0-9]{0,2}       # Element spec
                (
                    \s+                         # White space in front of the number
                    [-|+]?                      # Plus or minus in front of the number (optional)
                    (
                        (
                            \d*                 # optional decimal in the beginning .0001 is ok, for example
                            [\.]                # There has to be a dot followed by
                            \d+                 # at least one decimal
                        )
                        |                       # OR
                        (
                            \d+                 # at least one decimal, followed by
                            [\.]?               # an optional dot ( both 1 and 1. are fine)
                            \d*                 # And optional number of decimals (1.00001)
                        )                        # followed by optional decimals
                    )
                    ([E|e|d|D][+|-]?\d+)?       # optional exponents E+03, e-05
                ){3}                            # I expect three float values
                ((\s+[0-1]){3}\s*)?             # Followed by optional ifpos
                \s*                             # Followed by optional white space
                |
                \#.*                            # If a line is commented out, that is also ok
                |
                \!.*                            # Comments also with excl. mark in fortran
            )
            |                                   # OR
            \s*                                 # A line only containing white space
         )
        [\n]                                    # line break at the end
    )+                                          # A positions block should be one or more lines
)
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const CELL_PARAMETERS_BLOCK_REGEX = r"""
^ [ \t]*
CELL_PARAMETERS [ \t]*
[{(]? \s* (?P<units>[a-z]*) \s* [)}]? \s* [\n]
(?P<block>
(
    (
        \s*             # White space in front of the element spec is ok
        (
            # First number
            (
                [-|+]?   # Plus or minus in front of the number (optional)
                (\d*     # optional decimal in the beginning .0001 is ok, for example
                [\.]     # There has to be a dot followed by
                \d+)     # at least one decimal
                |        # OR
                (\d+     # at least one decimal, followed by
                [\.]?    # an optional dot
                \d*)     # followed by optional decimals
                ([E|e|d|D][+|-]?\d+)?  # optional exponents E+03, e-05, d0, D0
            
                (
                    \s+      # White space between numbers
                    [-|+]?   # Plus or minus in front of the number (optional)
                    (\d*     # optional decimal in the beginning .0001 is ok, for example
                    [\.]     # There has to be a dot followed by
                    \d+)     # at least one decimal
                    |        # OR
                    (\d+     # at least one decimal, followed by
                    [\.]?    # an optional dot
                    \d*)     # followed by optional decimals
                    ([E|e|d|D][+|-]?\d+)?  # optional exponents E+03, e-05, d0, D0
                ){2}         # I expect three float values
            )
            |
            \#
            |
            !            # If a line is commented out, that is also ok
        )
        .*               # I do not care what is after the comment or the vector
        |                # OR
        \s*              # A line only containing white space
     )
    [\n]                 # line break at the end
){3}                     # I need exactly 3 vectors
)
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const ATOMIC_SPECIES_BLOCK_REGEX = r"""
^ [ \t]* ATOMIC_SPECIES [ \t]* $\n
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* $\n?
 )+
)
"""imx
# This regex is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_SPECIAL_BLOCK_REGEX = r"""
^ [ \t]* K_POINTS [ \t]*
    [{(]? [ \t]* (?P<type>\S+?)? [ \t]* [)}]? [ \t]* $\n
^ [ \t]* \S+ [ \t]* $\n  # nks
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* $\n?
 )+
)
"""imx
# This regex is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_AUTOMATIC_BLOCK_REGEX = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* $\n
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+)
    [ \t]+ (\S+) [ \t]* $\n?
"""imx
# This regex is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_GAMMA_BLOCK_REGEX = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* $\n
"""imx

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

function findcards(str::AbstractString)
    captured = Dict{Symbol,String}()
    for (name, regex) in zip((:AtomicSpeciesCard, :AtomicPositionsCard), (ATOMIC_SPECIES_BLOCK_REGEX, ATOMIC_POSITIONS_BLOCK_REGEX))
        m = match(regex, str)
        @assert !isnothing(m) "Cannot find compulsory card $name. You must provide one!"
        push!(captured, name => m.captures)
    end
    for regex in (K_POINTS_AUTOMATIC_BLOCK_REGEX, K_POINTS_GAMMA_BLOCK_REGEX, K_POINTS_SPECIAL_BLOCK_REGEX)
        m = match(regex, str)
        isnothing(m.captures) ? continue : push!(captured, :KPointsCard => m.captures)
    end
    @assert !haskey(captured, :KPointsCard) "No `K_POINTS` card found! You must provide one!"
    for (name, regex) in zip((:CellParametersCard), (CELL_PARAMETERS_BLOCK_REGEX,))
        m = match(regex, str)
        # These are not compulsory cards, match failures will be allowed.
        push!(captured, name => m.captures)
    end
end # function findcards

end
