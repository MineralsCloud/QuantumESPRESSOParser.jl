module PHonon

using Compat: only
using QuantumESPRESSOBase.Inputs: Namelist
using QuantumESPRESSOBase.Inputs.PHonon:
    SpecialKPoint,
    QPointsCard,
    PhInput,
    Q2rInput,
    DynmatInput,
    MatdynInput,
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist

using PyFortran90Namelists: FortranData

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

function Base.parse(::Type{QPointsCard}, str::AbstractString)
    m = match(Q_POINTS_SPECIAL_BLOCK_REGEX, str)
    if m === nothing
        @info "no q-points are provided."
        return
    else
        captured = only(m.captures)
        data = map(eachmatch(Q_POINTS_SPECIAL_ITEM_REGEX, captured)) do matched
            # TODO: Match `nqs`
            SpecialKPoint(map(x -> parse(Float64, FortranData(x)), matched.captures))
        end
        return QPointsCard(data)
    end
end # function Base.parse
function Base.parse(::Type{PhInput}, str::AbstractString)
    inputph = parse(PhNamelist, str)
    q_points = parse(QPointsCard, str)
    return PhInput(inputph, q_points)
end # function Base.parse
Base.parse(::Type{Q2rInput}, str::AbstractString) = Q2rInput(parse(Q2rNamelist, str))
function Base.parse(::Type{MatdynInput}, str::AbstractString)
    input = parse(MatdynNamelist, str)
    q_points = parse(QPointsCard, str)
    return MatdynInput(input, q_points)
end # function Base.parse
Base.parse(::Type{DynmatInput}, str::AbstractString) =
    DynmatInput(parse(DynmatNamelist, str))

end
