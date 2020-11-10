"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using AbInitioSoftwareBase.Inputs: inputstring, groupname
using Compat: only
using PyFortran90Namelists: fparse
using QuantumESPRESSOBase.Inputs: Card
using QuantumESPRESSOBase.Inputs.PWscf:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    KPointsCard,
    GammaPointCard,
    KMeshCard,
    SpecialPointsCard,
    MonkhorstPackGrid,
    SpecialPoint,
    CellParametersCard,
    PWInput

export format_text, format_file

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
const ATOMIC_POSITIONS_BLOCK = r"""
^ \s* ATOMIC_POSITIONS \s*                      # Atomic positions start with that string
[{(]? \s* (?P<units>\S+?)? \s* [)}]? \s* $\R    # The units are after the string in optional brackets
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
        \R                                      # line break at the end
    )+                                          # A positions block should be one or more lines
)
"""imx
# This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py.
const CELL_PARAMETERS_BLOCK = r"""
^ [ \t]*
CELL_PARAMETERS [ \t]*
[{(]? \s* (?P<option>[a-z]*) \s* [)}]? \h* (?:[\#!].*)? \v  # Match option, do not match comment
(?P<data>
(?:
    \s*              # White space in front of the element spec is ok
    (?:
        [-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[ED][-+]?[0-9]+)?  # First number
        (?:
        \s+          # White space between numbers
        [-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[ED][-+]?[0-9]+)?
        ){2}         # I expect 3 numbers
        |
        \#           # If a line is commented out, that is also ok
        |
        !            # If a line is commented out, that is also ok
    )
    .*               # I do not care what is after the comment or the vector
    |                # Or
    ^\s*$            # A line only containing white space
){3}                 # I need exactly 3 vectors
)
"""imx
# This regular expression is referenced from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py#L683-L691
const ATOMIC_SPECIES_BLOCK = r"""
^ [ \t]* ATOMIC_SPECIES [ \t]* \R+
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ (?:[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+ \S+ [ \t]* \R?
 )+
)
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_SPECIAL_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]*
    [{(]? [ \t]* (?P<type>\S+?)? [ \t]* [)}]? [ \t]* \R+
^ [ \t]* \S+ [ \t]* \R+  # nks
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* \R+
 )+
)
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_AUTOMATIC_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* \R
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+)
    [ \t]+ (\S+) [ \t]* \R+
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_GAMMA_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* \R+
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
const ATOMIC_SPECIES_ITEM = r"""
^ [ \t]* (?P<name>\S+) [ \t]+ (?P<mass>[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+ (?P<pseudo>\S+)
    [ \t]* \R?
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
const ATOMIC_POSITIONS_ITEM = r"""
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
const K_POINTS_SPECIAL_ITEM = r"""
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]* \R?
"""mx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
const CELL_PARAMETERS_ITEM = r"""
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

function Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    if m !== nothing
        content = only(m.captures)
        return AtomicSpeciesCard(
            map(eachmatch(ATOMIC_SPECIES_ITEM, content)) do matched
                captured = matched.captures
                atom, mass, pseudopotential =
                    captured[1], fparse(Float64, captured[2]), captured[3]
                AtomicSpecies(atom, mass, pseudopotential)
            end,
        )
    end
end # function Base.tryparse
function Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    if m !== nothing
        if string(m.captures[1]) === nothing
            @warn "Not specifying units is DEPRECATED and will no longer be allowed in the future!"
            @info "No option is specified, 'alat' is assumed."
            option = "alat"
        else
            option = string(m.captures[1])
        end
        content = m.captures[2]
        return AtomicPositionsCard(
            map(eachmatch(ATOMIC_POSITIONS_ITEM, content)) do matched
                # The `matched` cannot be a `nothing` since we have tested by the block regular expression
                captured = matched.captures
                # The `if_pos` field is optionally given by users. If they do not give, we provide the default values `1`.
                if_pos = map(x -> isempty(x) ? 1 : fparse(Int, x), captured[11:13])
                # The `atom` and `pos` fields are mandatory. So we do not need special treatment.
                atom, pos = captured[1],
                map(x -> fparse(Float64, x), [captured[2], captured[5], captured[8]])
                AtomicPosition(atom, pos, if_pos)
            end,
            option,
        )
    end
end # function Base.tryparse
function Base.tryparse(::Type{GammaPointCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK, str)
    return m === nothing ? nothing : GammaPointCard()
end # function Base.tryparse
function Base.tryparse(::Type{KMeshCard}, str::AbstractString)
    m = match(K_POINTS_AUTOMATIC_BLOCK, str)
    if m !== nothing
        data = map(x -> fparse(Int, x), m.captures)
        return KMeshCard(MonkhorstPackGrid(data[1:3], data[4:6]))
    end
end # function Base.tryparse
function Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)
    m = match(K_POINTS_SPECIAL_BLOCK, str)
    if m !== nothing
        option = m.captures[1] === nothing ? "tpiba" : m.captures[1]
        return SpecialPointsCard(
            map(eachmatch(K_POINTS_SPECIAL_ITEM, m.captures[2])) do matched
                # TODO: Match `nks`
                SpecialPoint(map(x -> fparse(Float64, x), matched.captures)...)
            end,
            option,
        )
    end
end # function Base.tryparse
function Base.tryparse(::Type{KPointsCard}, str::AbstractString)
    for T in (GammaPointCard, KMeshCard, SpecialPointsCard)
        x = tryparse(T, str)
        if x !== nothing
            return x
        end
    end
end # function Base.tryparse
function Base.tryparse(::Type{CellParametersCard}, str::AbstractString)
    m = match(CELL_PARAMETERS_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    if m !== nothing
        option = string(m[:option])
        if isempty(option)
            @warn "Neither unit nor lattice parameter are specified. DEPRECATED, will no longer be allowed!"
            @info "'bohr' is assumed."
            option = "bohr"
        end
        content = m[:data]
        data = Matrix{Float64}(undef, 3, 3)
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM, content))
            captured = matched.captures
            data[i, :] =
                map(x -> fparse(Float64, x), [captured[1], captured[4], captured[7]])
        end
        return CellParametersCard(data, option)
    end
end # function Base.tryparse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Card}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find card `$(groupname(T))`!"))
    else
        return x
    end
end # function Base.parse
function Base.parse(::Type{PWInput}, str::AbstractString)
    args = []
    for T in (CellParametersCard,)  # ConstraintsCard, OccupationsCard, AtomicForcesCard
        x = tryparse(T, str)  # Optional cards
        if x !== nothing
            push!(args, x)
        end
    end
    for T in (AtomicSpeciesCard, AtomicPositionsCard, KPointsCard)  # Must-have cards, or else error
        push!(args, parse(T, str))
    end
    for T in
        (ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist)
        x = tryparse(T, str)
        push!(args, x === nothing ? T() : x)
    end
    return PWInput(args...)
end # function Base.parse

function format_file(filename::AbstractString; overwrite::Bool = true, kwargs...)
    text = read(filename, String)
    formatted_text = format_text(text; kwargs...)
    if overwrite
        open(filename, "w") do io
            write(io, formatted_text)
        end
    else
        println(formatted_text)
    end
end # function format_file

format_text(text::AbstractString) = inputstring(parse(PWInput, text))

end
