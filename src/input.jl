using AbInitioSoftwareBase: Namelist, groupname
using QuantumESPRESSOBase: QuantumESPRESSOInput
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
end

function Base.parse(::Type{T}, str::AbstractString) where {T<:Namelist}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find namelist `$(groupname(T))`!"))
    else
        return x
    end
end

# Idea from https://github.com/JuliaData/CSV.jl/blob/c3af297/src/CSV.jl#L64-L69
function Base.read(io::IO, ::Type{T}) where {T<:QuantumESPRESSOInput}
    str = read(io, String)
    return parse(T, str)
end
function Base.read(filename::AbstractString, ::Type{T}) where {T<:QuantumESPRESSOInput}
    str = read(filename, String)
    return parse(T, str)
end
