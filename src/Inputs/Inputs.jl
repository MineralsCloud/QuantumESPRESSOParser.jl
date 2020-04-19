"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using Compat: isnothing
using PyFortran90Namelists: FortranData, Parser
using QuantumESPRESSOBase.Inputs: Namelist, InputEntry, Input, titleof, qestring

export InputFile

struct InvalidInput
    msg::String
end

struct InputFile{A}
    source::A
end

Base.tryparse(::Type{T}, f::InputFile) where {T<:InputEntry} = tryparse(T, read(f))
function Base.tryparse(::Type{T}, str::AbstractString) where {T<:Namelist}
    d::Dict{String,Any} = Parser().reads(str)
    return if haskey(d, lowercase(titleof(T)))
        dict = Dict(Symbol(k) => v for (k, v) in d[lowercase(titleof(T))])
        T(; dict...)
    end
end # function Base.tryparse

Base.parse(::Type{T}, f::InputFile) where {T<:InputEntry} = parse(T, read(f))
function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    x = tryparse(T, str)
    isnothing(x) ? throw(Meta.ParseError("cannot find namelist `$(titleof(T))`!")) : x
end # function Base.parse

Base.read(f::InputFile{String}) = read(f.source, String)

Base.write(f::InputFile{String}, x::InputEntry) = write(f.source, qestring(x))
Base.write(f::InputFile{String}, x::Input) = write(f.source, qestring(x))

include("PWscf.jl")
# include("PHonon.jl")

end
