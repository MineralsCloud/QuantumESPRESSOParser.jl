
using AbInitioSoftwareBase: Namelist, groupname
using PyFortran90Namelists: fparse
using QuantumESPRESSOBase.PHonon:
    SpecialPoint,
    QPointsCard,
    PhInput,
    Q2rInput,
    DynmatInput,
    MatdynInput,
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist

const Q_POINTS_SPECIAL_BLOCK_REGEX = r"""
^ [ \t]* qPointsSpecs [ \t]*$\R
^ [ \t]* \S+ [ \t]* $\R  # nqs
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* $\R?
 )+
)
"""imx
const Q_POINTS_SPECIAL_ITEM_REGEX = r"""
^ [ \t]*
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\h*
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\h*
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\h*
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\h*
"""imx

function Base.tryparse(::Type{QPointsCard}, str::AbstractString)
    m = match(Q_POINTS_SPECIAL_BLOCK_REGEX, str)
    if m === nothing
        @info "no q-points are provided."
        return nothing
    else
        captured = only(m.captures)
        data = map(eachmatch(Q_POINTS_SPECIAL_ITEM_REGEX, captured)) do matched
            # TODO: Match `nqs`
            SpecialPoint(map(x -> fparse(Float64, x), matched.captures))
        end
        return QPointsCard(data)
    end
end
function Base.tryparse(::Type{PhInput}, str::AbstractString)
    title_line = ""
    for line in split(str, r"\r\n|\r|\n")
        if !isempty(line)
            title_line = line
            break
        end
    end
    args = []
    for T in (PhNamelist, QPointsCard)
        push!(args, tryparse(T, str))
    end
    if all(x === nothing for x in args)
        return nothing
    else
        return PhInput(title_line, args[1], args[2])
    end
end
function Base.tryparse(::Type{Q2rInput}, str::AbstractString)
    x = parse(Q2rNamelist, str)
    if x !== nothing
        return Q2rInput(parse(Q2rNamelist, str))
    else
        return nothing
    end
end
function Base.parse(::Type{MatdynInput}, str::AbstractString)
    args = []
    for T in (MatdynNamelist, QPointsCard)
        push!(args, tryparse(T, str))
    end
    if all(x === nothing for x in args)
        return nothing
    else
        return MatdynInput(args[1], args[2])
    end
end
function Base.tryparse(::Type{DynmatInput}, str::AbstractString)
    if parse(DynmatNamelist, str)
        return DynmatInput(parse(DynmatNamelist, str))
    else
        return nothing
    end
end

function Base.parse(
    ::Type{T}, str::AbstractString
) where {T<:Union{QPointsCard,PhInput,Q2rInput,DynmatInput,MatdynInput}}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find card `$(groupname(T))`!"))
    else
        return x
    end
end
