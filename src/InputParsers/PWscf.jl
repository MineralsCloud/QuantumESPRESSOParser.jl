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

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Namelists: Namelist, to_dict
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs.PWscf

using QuantumESPRESSOParsers
using QuantumESPRESSOParsers.InputParsers.Namelists

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
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
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
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_AUTOMATIC_BLOCK_REGEX = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* $\n
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+)
    [ \t]+ (\S+) [ \t]* $\n?
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_GAMMA_BLOCK_REGEX = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* $\n
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const ATOMIC_SPECIES_ITEM_REGEX = r"""
^ [ \t]* (?P<name>\S+) [ \t]+ (?P<mass>\S+) [ \t]+ (?P<pseudo>\S+)
    [ \t]* $\n?
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const ATOMIC_POSITIONS_ITEM_REGEX = r"""
^                                       # Linestart
[ \t]*                                  # Optional white space
(?P<name>[A-Za-z]+[A-Za-z0-9]{0,2})\s+   # get the symbol, max 3 chars, starting with a char
(?P<x>                                  # Get x
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                                  # Get y
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                                  # Get z
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]*
(?P<fx>[01]?)                           # Get fx
[ \t]*
(?P<fy>[01]?)                           # Get fx
[ \t]*
(?P<fz>[01]?)                           # Get fx
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_SPECIAL_ITEM_REGEX = r"""
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]* $\n?
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const CELL_PARAMETERS_ITEM_REGEX = r"""
^                        # Linestart
[ \t]*                   # Optional white space
(?P<x>                   # Get x
    [\-|\+]? ( \d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                   # Get y
    [\-|\+]? (\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                   # Get z
    [\-|\+]? (\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
"""mx

function Base.parse(::Type{<:AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK_REGEX, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    @assert !isnothing(m) "Cannot find card `ATOMIC_SPECIES`! Check your input!"
    content = m.captures[1]
    data = AtomicSpecies[]
    for matched in eachmatch(ATOMIC_SPECIES_ITEM_REGEX, content)
        captured = matched.captures
        atom, mass, pseudopotential = string(captured[1]), parse(Float64, FortranData(captured[2])), string(captured[3])
        push!(data, AtomicSpecies(atom, mass, pseudopotential))
    end
    return AtomicSpeciesCard(data)
end # function Base.parse
function Base.parse(T::Type{<:AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK_REGEX, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    @assert !isnothing(m) "Cannot find card `ATOMIC_POSITIONS`! Check your input!"
    option = string(m.captures[1])
    if isnothing(option)
        @warn "Not specifying units is DEPRECATED and will no longer be allowed in the future!"
        @info "No option is specified, 'alat' is assumed."
        option = "alat"
    end
    content = m.captures[2]
    data = AtomicPosition{String,Vector{Float64},Vector{Int}}[]
    for matched in eachmatch(ATOMIC_POSITIONS_ITEM_REGEX, content)
        # The `matched` cannot be a `nothing` since we have tested by the block regular expression
        captured = matched.captures
        # The `if_pos` field is optionally given by users. If they do not give, we provide the default values `1`.
        if_pos = map(x -> isempty(x) ? 1 : parse(Int, FortranData(x)), captured[11:13])
        # The `atom` and `pos` fields are mandatory. So we do not need special treatment.
        atom, pos = string(captured[1]),
            map(x -> parse(Float64, FortranData(x)), [captured[2], captured[5], captured[8]])
        push!(data, AtomicPosition(atom, pos, if_pos))
    end
    return AtomicPositionsCard(option, data)
end # function Base.parse
function Base.parse(::Type{<:KPointsCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK_REGEX, str)
    !isnothing(m) && return KPointsCard("gamma", [GammaPoint()])

    m = match(K_POINTS_AUTOMATIC_BLOCK_REGEX, str)
    if !isnothing(m)
        data = map(x -> parse(Int, FortranData(x)), m.captures)
        return KPointsCard("automatic", [MonkhorstPackGrid(data[1:3], data[4:6])])
    end

    m = match(K_POINTS_SPECIAL_BLOCK_REGEX, str)
    if !isnothing(m)
        option = m.captures[1]
        captured = m.captures[2]
        data = SpecialKPoint[]
        for matched in eachmatch(K_POINTS_SPECIAL_ITEM_REGEX, captured)
            # TODO: Match `nks`
            point = @match map(x -> parse(Float64, FortranData(x)), matched.captures) begin
                [coordinates..., weight] => SpecialKPoint(coordinates, weight)
            end
            push!(data, point)
        end
        return KPointsCard(option, data)
    end

    @info "Cannot find card `K_POINTS`!"
end # function Base.parse
function Base.parse(::Type{<:CellParametersCard}, str::AbstractString)
    m = match(CELL_PARAMETERS_BLOCK_REGEX, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    isnothing(m) && return nothing
    option = string(m.captures[1])
    if isnothing(option)
        @warn "Neither unit nor lattice parameter are specified. DEPRECATED, will no longer be allowed!"
        @info "'bohr' is assumed."
        option = "bohr"
    end
    content = m.captures[2]
    data = Matrix{Float64}(undef, 3, 3)
    for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM_REGEX, content))
        captured = matched.captures
        data[i, :] = map(x -> parse(Float64, FortranData(x)), [captured[1], captured[4], captured[7]])
    end
    return CellParametersCard(option, data)
end # function Base.parse

end
