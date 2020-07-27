"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using PyFortran90Namelists: FortranData, Parser
using QuantumESPRESSOBase.Inputs:
    Namelist, QuantumESPRESSOInputEntry, Input, titleof, inputstring

struct InvalidInput
    msg::String
end

function Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist}
    d::Dict{String,Any} = Parser().reads(str)
    return if haskey(d, lowercase(titleof(T)))
        dict = Dict(Symbol(k) => v for (k, v) in d[lowercase(titleof(T))])
        T(; dict...)
    end
end # function Base.tryparse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find namelist `$(titleof(T))`!"))
    else
        return x
    end
end # function Base.parse

include("PWscf.jl")
include("PHonon.jl")

end
