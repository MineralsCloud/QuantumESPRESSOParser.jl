"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using QuantumESPRESSOBase: asfieldname
using QuantumESPRESSOBase.Namelists.PWscf: ControlNamelist,
                                           SystemNamelist,
                                           ElectronsNamelist,
                                           IonsNamelist,
                                           CellNamelist
using QuantumESPRESSOBase.Cards.PWscf: AtomicSpeciesCard,
                                       AtomicPositionsCard,
                                       KPointsCard,
                                       CellParametersCard
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput

function Base.parse(::Type{<:PWscfInput}, str::AbstractString)
    dict = Dict{Symbol,Any}()
    for T in (CellParametersCard,)  # ConstraintsCard, OccupationsCard, AtomicForcesCard
        push!(dict, asfieldname(T) => tryparse(T, str))  # Optional cards, can be `nothing`
    end
    for T in (AtomicSpeciesCard, AtomicPositionsCard, KPointsCard)
        push!(dict, asfieldname(T) => parse(T, str))  # Must-have cards, or else error
    end
    for T in (ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist)
        nml = tryparse(T, str)
        push!(dict, asfieldname(T) => isnothing(nml) ? T() : nml)
    end
    return PWscfInput(; dict...)
end # function Base.parse

end
