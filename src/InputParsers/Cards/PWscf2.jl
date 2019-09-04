"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Fortran90Namelists.FortranToJulia
using QuantumESPRESSOBase.Cards.PWscf

const ATOMIC_SPECIES_ITEM_REGEX = r"""
^ [ \t]* (?P<name>\S+) [ \t]+ (?P<mass>\S+) [ \t]+ (?P<pseudo>\S+)
    [ \t]* $\n?
"""mx

function Base.parse(::Type{<:AtomicSpeciesCard}, str::AbstractString)
    data = AtomicSpecies[]
    for m in eachmatch(ATOMIC_SPECIES_ITEM_REGEX, str)
        captured = m.captures
        atom, mass, pseudopotential = string(captured[1]), parse(Float64, FortranData(captured[2])), string(captured[3])
        push!(data, AtomicSpecies(atom, mass, pseudopotential))
    end
    return AtomicSpeciesCard(data)
end # function Base.parse

end