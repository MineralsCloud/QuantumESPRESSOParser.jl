export eachcellparameterscard, eachatomicpositionscard

# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/output_tau.f90#L47-L60.
const CELL_PARAMETERS_BLOCK_OUTPUT = r"""
CELL_PARAMETERS \h+
\( (?<option>\w+) =? \s* (?<alat>[-+]?[0-9]*\.[0-9]{8})? \) \h*  # Match `alat`: `F12.8`
(?<data>
    (?: \s*
        (?:
            [-+]?[0-9]*\.[0-9]+ \s*  # Match element
        ){3}  # I need exactly 3 elements per vector
    ){3}  # I need exactly 3 vectors
)
"""x
const CELL_PARAMETERS_ITEM_OUTPUT = r"""
\s*
([-+]?[0-9]*\.[0-9]+) \s*  # x
([-+]?[0-9]*\.[0-9]+) \s*  # y
([-+]?[0-9]*\.[0-9]+) \s*  # z
"""x
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/output_tau.f90#L64-L109.
const ATOMIC_POSITIONS_BLOCK_OUTPUT = r"""
ATOMIC_POSITIONS \h*                   # Atomic positions start with that string
\( (?<option>\w+) \)                   # Option of the card
(?<data>
    (?:
        \s*
        [A-Za-z]+[A-Za-z0-9]{0,2} \s+  # Atom spec
        (?:
            [-+]?[0-9]*\.[0-9]+ \s*  # Match element
        ){3}                           # I need exactly 3 floats per vector.
        (?:
            [-+]?[0-9]+ \s*
        ){0,3}                     # I need exactly 3 integers in `if_pos`, if there is any.
    )+
)
"""x
const ATOMIC_POSITIONS_ITEM_OUTPUT = r"""
\s*
([A-Za-z]+[A-Za-z0-9]{0,2}) \s+  # Atom spec
([-+]?[0-9]*\.[0-9]+) \s*  # x
([-+]?[0-9]*\.[0-9]+) \s*  # y
([-+]?[0-9]*\.[0-9]+) \s*  # z
([-+]?[0-9]+)? \s*            # if_pos(1)
([-+]?[0-9]+)? \s*            # if_pos(2)
([-+]?[0-9]+)? \s*            # if_pos(3)
"""x

function Base.parse(::Type{CellParametersCard}, str::AbstractString)
    obj = _tryparse(CellParametersCard, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.parse(::Type{AtomicPositionsCard}, str::AbstractString)
    obj = _tryparse(AtomicPositionsCard, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function _tryparse(::Type{CellParametersCard}, str::AbstractString)
    matched = match(CELL_PARAMETERS_BLOCK_OUTPUT, str)
    if isnothing(matched)
        return nothing
    else
        body, data = matched[:data], Matrix{Float64}(undef, 3, 3)  # Initialization
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM_OUTPUT, body))
            data[i, :] = map(x -> parse(Float64, x), matched.captures)
        end
        option = Symbol(matched[:option])
        if option == :alat
            alat = parse(Float64, matched[:alat])
            return CellParametersCard(alat * data, :bohr)
        else
            return CellParametersCard(data, option)
        end
    end
end
function _tryparse(::Type{AtomicPositionsCard}, str::AbstractString)
    matched = match(ATOMIC_POSITIONS_BLOCK_OUTPUT, str)
    if isnothing(matched)
        return nothing
    else
        option = Symbol(matched[1])
        body = matched[2]
        data = AtomicPosition[]
        for matched in eachmatch(ATOMIC_POSITIONS_ITEM_OUTPUT, body)
            captured = matched.captures
            if_pos = map(x -> x === nothing ? 1 : parse(Int, x), captured[5:7])
            atom, pos = captured[1], map(x -> parse(Float64, x), captured[2:4])
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        return AtomicPositionsCard(data, option)
    end
end

function Base.iterate(regexmatchiterator::EachParsed{CellParametersCard})
    regexmatchiterator = eachmatch(regexmatchiterator.regex, regexmatchiterator.string)
    iterated = iterate(regexmatchiterator)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return _tryparse(CellParametersCard, matched.match), state
    end
end
function Base.iterate(regexmatchiterator::EachParsed{CellParametersCard}, state)
    regexmatchiterator = eachmatch(regexmatchiterator.regex, regexmatchiterator.string)
    iterated = iterate(regexmatchiterator, state)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return _tryparse(CellParametersCard, matched.match), state
    end
end

eachcellparameterscard(str::AbstractString) =
    EachParsed{CellParametersCard}(CELL_PARAMETERS_BLOCK_OUTPUT, str)

eachatomicpositionscard(str::AbstractString) =
    EachParsed{AtomicPositionsCard}(ATOMIC_POSITIONS_BLOCK_OUTPUT, str)