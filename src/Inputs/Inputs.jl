"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using AbInitioSoftwareBase.Inputs: Namelist, groupname
using PyFortran90Namelists: Parser

struct InvalidInput
    msg::String
end

function Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist}
    d::Dict{String,Any} = Parser().reads(str)
    return if haskey(d, lowercase(groupname(T)))
        dict = Dict(Symbol(k) => v for (k, v) in d[lowercase(groupname(T))])
        T(; dict...)
    end
end # function Base.tryparse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find namelist `$(groupname(T))`!"))
    else
        return x
    end
end # function Base.parse

include("PWscf.jl")
include("PHonon.jl")

end
