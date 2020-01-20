module PHonon

using QuantumESPRESSOBase.Namelists: Namelist
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase.Cards.PHonon: QPointsSpecsCard
using QuantumESPRESSOBase.Inputs.PHonon

using QuantumESPRESSOParsers.Namelists
using QuantumESPRESSOParsers.Cards

function Base.parse(::Type{T}, str::AbstractString) where {T<:PhInput}
    inputph = parse(PhNamelist, str)
    q_points = parse(QPointsSpecsCard, str)
    return T(inputph, q_points)
end # function Base.parse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Q2rInput}
    input = parse(Q2rNamelist, str)
    return T(input)
end # function Base.parse

function Base.parse(::Type{T}, str::AbstractString) where {T<:MatdynInput}
    input = parse(MatdynNamelist, str)
    q_points = parse(QPointsSpecsCard, str)
    return T(input, q_points)
end # function Base.parse

function Base.parse(::Type{T}, str::AbstractString) where {T<:DynmatInput}
    input = parse(DynmatNamelist, str)
    return T(input)
end # function Base.parse

end
