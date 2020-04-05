"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using Compat: isnothing
using FilePaths: AbstractPath
using PyFortran90Namelists: FortranData, Parser
using QuantumESPRESSOBase.Inputs: Namelist, titleof

struct InvalidInput
    msg::String
end

struct InputString <: AbstractString
    str::String
end

Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist} =
    tryparse(T, InputString(str))
function Base.tryparse(::Type{T}, str::InputString) where {T<:Namelist}
    str = str.str
    d::Dict{String,Any} = Parser().reads(str)
    return if haskey(d, lowercase(titleof(T)))
        dict = Dict(Symbol(k) => v for (k, v) in d[lowercase(titleof(T))])
        T(; dict...)
    end
end # function Base.tryparse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    x = tryparse(T, str)
    isnothing(x) ? throw(Meta.ParseError("cannot find namelist `$(titleof(T))`!")) : x
end # function Base.parse
Base.parse(::Type{T}, fp::AbstractPath) where {T<:Namelist} = parse(T, read(fp, String))

include("PWscf.jl")
# include("PHonon.jl")

end
