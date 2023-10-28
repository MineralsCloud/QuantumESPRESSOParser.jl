using Dates: Hour, Minute, Millisecond
using QuantumESPRESSOBase.PWscf
using VersionParsing: vparse

export Preamble,
    FFTGrid,
    IrreducibleBrillouinZone,
    IterationHead,
    UnconvergedEnergy,
    ConvergedEnergy,
    TimedItem
export parse_symmetries, parse_stress, parse_bands

struct ParseError <: Exception
    msg::String
end

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported

abstract type PWOutputItem end

include("regexes.jl")
include("once.jl")
include("each.jl")
include("atomicstructure.jl")
include("misc.jl")

function parse_symmetries(str::AbstractString)
    m = match(SYM_OPS, str)
    m === nothing && return nothing
    return num_sym_ops = isempty(m[:n]) ? 0 : parse(Int, m[:n])
end # function parse_symmetries

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses, kbar_stresses = Matrix{Float64}[], Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK, str)
        pressure, content = m.captures[1], m.captures[3]
        push!(pressures, parse(Float64, pressure))

        stress_atomic, stress_kbar = ntuple(_ -> Matrix{Float64}(undef, 3, 3), 2)
        for (i, line) in enumerate(split(content, r"\R+"))
            tmp = map(x -> parse(Float64, x), split(line, " "; keepempty=false))
            stress_atomic[i, :], stress_kbar[i, :] = tmp[1:3], tmp[4:6]
        end
        push!(atomic_stresses, stress_atomic)
        push!(kbar_stresses, stress_kbar)
    end
    return pressures, atomic_stresses, kbar_stresses
end # function parse_stress

# See https://github.com/QEF/q-e/blob/4132a64/PW/src/print_ks_energies.f90#L10.
function parse_bands(str::AbstractString)
    str == "Number of k-points >= 100: set verbosity='high' to print the bands." &&
        return nothing
    kpts, bands = nothing, nothing  # Initialization
    m = match(KS_ENERGIES_BLOCK, str)
    if m !== nothing
        kpts, bands = Vector{Float64}[], Vector{Float64}[]
        regex = if match(KS_ENERGIES_BANDS, str) === nothing
            KS_ENERGIES_BAND_ENERGIES
        else
            KS_ENERGIES_BANDS
        end
        for m in eachmatch(regex, str)
            push!(
                kpts, map(x -> parse(Float64, x[1]), eachmatch(Regex(GENERAL_REAL), m[:k]))
            )
            push!(
                bands,
                map(x -> parse(Float64, x[1]), eachmatch(Regex(GENERAL_REAL), m[:band])),
            )
        end
        len, nbnd = length(kpts), length(bands[1])
        kpts, bands = reshape(collect(Iterators.flatten(kpts)), len, 3),
        reshape(collect(Iterators.flatten(bands)), len, nbnd)
    end  # Keep them `nothing` if `m` is `nothing`
    return kpts, bands
end # function parse_bands

struct QuantumESPRESSOBaseVersion <: PWOutputItem
    version::VersionNumber
end

function Base.parse(::Type{QuantumESPRESSOBaseVersion}, str::AbstractString)
    obj = tryparse(QuantumESPRESSOBaseVersion, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{QuantumESPRESSOBaseVersion}, str::AbstractString)
    matched = match(PWSCF_VERSION, str)
    return isnothing(matched) ? nothing : vparse(matched[:version])
end

struct Parallelization <: PWOutputItem
    type::String
    cores::Maybe{Int64}
end

function Base.parse(::Type{Parallelization}, str::AbstractString)
    obj = tryparse(Parallelization, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{Parallelization}, str::AbstractString)
    matched = match(PARALLEL_INFO, str)
    if isnothing(matched)
        return nothing
    else
        return Parallelization(
            matched[:kind], isnothing(matched[:num]) ? 1 : parse(Int64, matched[:num])
        )
    end
end

struct FFTDimensions <: PWOutputItem
    ng::Int64
    nr1::Int64
    nr2::Int64
    nr3::Int64
end

function Base.parse(::Type{FFTDimensions}, str::AbstractString)
    obj = tryparse(FFTDimensions, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{FFTDimensions}, str::AbstractString)
    matched = match(FFT_DIMENSIONS, str)
    if isnothing(matched)
        return nothing
    else
        return FFTDimensions(map(Base.Fix1(parse, Int64), matched.captures)...)
    end
end

struct InputFile <: PWOutputItem
    name::String
end

function Base.parse(::Type{InputFile}, str::AbstractString)
    obj = tryparse(InputFile, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{InputFile}, str::AbstractString)
    matched = match(READING_INPUT_FROM, str)
    return isnothing(matched) ? nothing : InputFile(only(matched))
end
