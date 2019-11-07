module PHonon

using Rematch: @match

using QuantumESPRESSOBase.Namelists: Namelist
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase.Cards.PHonon: QPointsSpecsCard
using QuantumESPRESSOBase.Inputs.PHonon

using QuantumESPRESSOParsers.Namelists
using QuantumESPRESSOParsers.Cards

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
