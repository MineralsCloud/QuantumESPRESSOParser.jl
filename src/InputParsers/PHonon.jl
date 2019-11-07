module PHonon

using Rematch: @match

using QuantumESPRESSOBase.Namelists: Namelist, to_dict
using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards.PHonon: SpecialQPoint, QPointsSpecsCard
using QuantumESPRESSOBase.Inputs.PHonon

using QuantumESPRESSOParsers.InputParsers.Namelists

const Q_POINTS_SPECIAL_BLOCK_REGEX = r"""
^ [ \t]* qPointsSpecs [ \t]*$\n
^ [ \t]* \S+ [ \t]* $\n  # nqs
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* $\n?
 )+
)
"""imx

const Q_POINTS_SPECIAL_ITEM_REGEX = r"""
^ [ \t]* 
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\s+
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\s+
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\s+
([\-|\+]?(?:\d*[\.]\d+ | \d+[\.]?\d*)
    (?:[E|e|d|D][+|-]?\d+)?)\s*
"""imx

function Base.parse(::Type{<:QPointsSpecsCard}, str::AbstractString)
    m = match(Q_POINTS_SPECIAL_BLOCK_REGEX, str)
    if !isnothing(m)
        captured = m.captures[1]
        data = SpecialQPoint[]
        for matched in eachmatch(Q_POINTS_SPECIAL_ITEM_REGEX, captured)
            # TODO: Match `nqs`
            point = @match map(x -> parse(Float64, FortranData(x)), matched.captures) begin
                [coordinates..., weight] => SpecialQPoint(coordinates, weight)
            end
            push!(data, point)
        end
        return QPointsSpecsCard(data)
    end
    @info "Cannot find card `Q_POINTS`!"
end # function Base.parse

function Base.parse(::Type{PHononInput}, str::AbstractString)
    inputph = parse(PHNamelist, str)
    q_points = parse(QPointsSpecsCard, str)
    return PHononInput(inputph, q_points)
end # function Base.parse

function Base.parse(::Type{Q2RInput}, str::AbstractString)
    input = parse(Q2RNamelist, str)
    return Q2RInput(input)
end # function Base.parse

function Base.parse(::Type{MatdynInput}, str::AbstractString)
    input = parse(MatdynNamelist, str)
    q_points = parse(QPointsSpecsCard, str)
    return MatdynInput(input, q_points)
end # function Base.parse

function Base.parse(::Type{DynmatInput}, str::AbstractString)
    input = parse(DynmatNamelist, str)
    return DynmatInput(input)
end # function Base.parse

end
