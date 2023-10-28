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
export parse_symmetries,
    parse_stress,
    parse_bands,
    parse_version,
    parse_parallel_info,
    parse_fft_dimensions,
    parse_input_name,
    isoptimized,
    isjobdone

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

function parse_version(str::AbstractString)::Maybe{VersionNumber}
    m = match(PWSCF_VERSION, str)
    m !== nothing ? vparse(m[:version]) : return nothing
end # function parse_version

function parse_parallel_info(str::AbstractString)::Maybe{Tuple{String,Int}}
    m = match(PARALLEL_INFO, str)
    m === nothing && return nothing
    return m[:kind], m[:num] === nothing ? 1 : parse(Int, m[:num])
end # function parse_parallel_info

function parse_fft_dimensions(str::AbstractString)::Maybe{NamedTuple}
    m = match(FFT_DIMENSIONS, str)
    m === nothing && return nothing
    parsed = map(x -> parse(Int, x), m.captures)
    return (; zip((:ng, :nr1, :nr2, :nr3), parsed)...)
end # function parse_fft_dimensions

function parse_input_name(str::AbstractString)
    m = match(READING_INPUT_FROM, str)
    return m === nothing ? nothing : only(m)
end # function parse_input_name

isoptimized(str::AbstractString) =
    match(FINAL_COORDINATES_BLOCK, str) === nothing ? false : true

isjobdone(str::AbstractString) = match(JOB_DONE, str) !== nothing