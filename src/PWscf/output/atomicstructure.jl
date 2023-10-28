function tryparse_internal(::Type{CellParametersCard}, str::AbstractString)
    m = match(CELL_PARAMETERS_BLOCK_OUTPUT, str)
    return if m !== nothing
        body, data = m[:data], Matrix{Float64}(undef, 3, 3)  # Initialization
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM_OUTPUT, body))
            data[i, :] = map(x -> parse(Float64, x), matched.captures)
        end
        option = Symbol(m[:option])
        if option == :alat
            alat = parse(Float64, m[:alat])
            CellParametersCard(alat * data, :bohr)
        else
            CellParametersCard(data, option)
        end
    end
end # function tryparse_internal
function tryparse_internal(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK_OUTPUT, str)
    return if m !== nothing
        option = Symbol(m[1])
        body = m[2]
        data = AtomicPosition[]
        for matched in eachmatch(ATOMIC_POSITIONS_ITEM_OUTPUT, body)
            captured = matched.captures
            if_pos = map(x -> x === nothing ? 1 : parse(Int, x), captured[5:7])
            atom, pos = captured[1], map(x -> parse(Float64, x), captured[2:4])
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        AtomicPositionsCard(data, option)
    end
end # function tryparse_internal

const AtomicStructure = Union{CellParametersCard,AtomicPositionsCard}

function _parse(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    x = tryparse(T, str)
    return x === nothing ? throw(Meta.ParseError("cannot find `$(T)`!")) : x
end # function _parse

const REGEXOF = (
    CellParametersCard=CELL_PARAMETERS_BLOCK_OUTPUT,
    AtomicPositionsCard=ATOMIC_POSITIONS_BLOCK_OUTPUT,
)

tryparsefirst(::Type{T}, str::AbstractString) where {T<:AtomicStructure} =
    tryparse_internal(T, str)
parsefirst(::Type{T}, str::AbstractString) where {T<:AtomicStructure} = _parse(T, str)

function tryparseall(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    return map(eachmatch(REGEXOF[nameof(T)], str)) do x
        try
            tryparse_internal(T, x.match)
        catch
            nothing
        end
    end
end # function parseall
function parseall(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    return map(eachmatch(REGEXOF[nameof(T)], str)) do x
        try
            tryparse_internal(T, x.match)
        catch
            Meta.ParseError("Pass failed!")
        end
    end
end # function parseall

tryparselast(::Type{T}, str::AbstractString) where {T<:AtomicStructure} =
    tryparseall(T, str)[end]
parselast(::Type{T}, str::AbstractString) where {T<:AtomicStructure} = parseall(T, str)[end]

function _parsenext_internal(
    ::Type{T}, str::AbstractString, start::Integer, raise::Bool
) where {T}
    x = findnext(REGEXOF[nameof(T)], str, start)
    if x === nothing
        raise ? throw(Meta.ParseError("Nothing found for next!")) : return nothing
    end
    return tryparse_internal(T, str[x])
end # function parsenext
tryparsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, false)
parsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, true)

function tryparsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && return nothing
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && return nothing
    return tryparse_internal(T, m.match)
end # function parsefinal
function parsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && throw(Meta.ParseError("No final coordinates found!"))
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && throw(Meta.ParseError("No `CELL_PARAMETERS` found!"))
    return tryparse_internal(T, m.match)
end # function parsefinal
