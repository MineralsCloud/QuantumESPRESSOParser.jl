using Dates: Hour, Minute, Millisecond
using DataFrames: AbstractDataFrame, DataFrame, groupby
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
    parse_all_electron_energy,
    parse_energy_decomposition,
    parse_paw_contribution,
    parse_smearing_energy,
    parse_version,
    parse_parallel_info,
    parse_fft_dimensions,
    parse_input_name,
    isoptimized,
    isjobdone,
    tryparsefirst,
    parsefirst,
    tryparseall,
    parseall,
    tryparselast,
    parselast,
    tryparsenext,
    parsenext,
    tryparsefinal,
    parsefinal

struct ParseError <: Exception
    msg::String
end

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported

abstract type PWOutputItem end

include("regexes.jl")
include("once.jl")
include("each.jl")

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

function parse_all_electron_energy(str::AbstractString)
    df = DataFrame(; step=Int[], ae=Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        ae = if any(==(nothing), (m, m[:ae]))
            nothing
        else
            parse(Float64, match(Regex(FIXED_POINT_REAL), m[:ae])[1])
        end
        push!(df, [i ae])
    end
    return df
end # function parse_all_electron_energy

function parse_energy_decomposition(str::AbstractString)
    df = DataFrame(;
        step=Int[],
        onelectron=Maybe{Float64}[],
        hartree=Maybe{Float64}[],
        xc=Maybe{Float64}[],
        ewald=Maybe{Float64}[],
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if any(==(nothing), (m, m[:decomp]))
            ntuple(_ -> nothing, 4)
        else
            map(x -> parse(Float64, x[1]), eachmatch(Regex(FIXED_POINT_REAL), m[:decomp]))
        end
        push!(df, [i data...])
    end
    return df
end # function parse_energy_decomposition

function parse_paw_contribution(str::AbstractString)
    df = DataFrame(;
        step=Int[],
        one_electron=Maybe{Float64}[],
        hartree_ae=Maybe{Float64}[],
        hartree_ps=Maybe{Float64}[],
        xc_ae=Maybe{Float64}[],
        xc_ps=Maybe{Float64}[],
        eh=Maybe{Float64}[],
        exc=Maybe{Float64}[],
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if any(==(nothing), (m, m[:one]))
            ntuple(_ -> nothing, 6)
        else
            map(x -> parse(Float64, x[1]), eachmatch(Regex(FIXED_POINT_REAL), m[:one]))
        end
        push!(df, [i data...])
    end
    return df
end # function parse_paw_contribution

function parse_smearing_energy(str::AbstractString)
    df = DataFrame(; step=Int[], smearing=Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        smearing = if any(==(nothing), (m, m[:smearing]))
            nothing
        else
            parse(Float64, match(Regex(FIXED_POINT_REAL), m[:smearing])[1])
        end
        push!(df, [i smearing])
    end
    return df
end # function parse_smearing_energy

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

function tryparse_internal(::Type{CellParametersCard}, str::AbstractString)
    m = match(CELL_PARAMETERS_BLOCK_OUTPUT, str)
    return if m !== nothing
        body, data = m[:data], Matrix{Float64}(undef, 3, 3)  # Initialization
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM_OUTPUT, body))
            data[i, :] = map(x -> parse(Float64, x), matched.captures)
        end
        option = Symbol(m[:option])
        if option == :alat
            alat = parse(Float64, m[:alat])
            CellParametersCard(alat * data, :bohr)
        else
            CellParametersCard(data, option)
        end
    end
end # function tryparse_internal
function tryparse_internal(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK_OUTPUT, str)
    return if m !== nothing
        option = Symbol(m[1])
        body = m[2]
        data = AtomicPosition[]
        for matched in eachmatch(ATOMIC_POSITIONS_ITEM_OUTPUT, body)
            captured = matched.captures
            if_pos = map(x -> x === nothing ? 1 : parse(Int, x), captured[5:7])
            atom, pos = captured[1], map(x -> parse(Float64, x), captured[2:4])
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        AtomicPositionsCard(data, option)
    end
end # function tryparse_internal

const AtomicStructure = Union{CellParametersCard,AtomicPositionsCard}

function _parse(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    x = tryparse(T, str)
    return x === nothing ? throw(Meta.ParseError("cannot find `$(T)`!")) : x
end # function _parse

const REGEXOF = (
    CellParametersCard=CELL_PARAMETERS_BLOCK_OUTPUT,
    AtomicPositionsCard=ATOMIC_POSITIONS_BLOCK_OUTPUT,
)

tryparsefirst(::Type{T}, str::AbstractString) where {T<:AtomicStructure} =
    tryparse_internal(T, str)
parsefirst(::Type{T}, str::AbstractString) where {T<:AtomicStructure} = _parse(T, str)

function tryparseall(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    return map(eachmatch(REGEXOF[nameof(T)], str)) do x
        try
            tryparse_internal(T, x.match)
        catch
            nothing
        end
    end
end # function parseall
function parseall(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    return map(eachmatch(REGEXOF[nameof(T)], str)) do x
        try
            tryparse_internal(T, x.match)
        catch
            Meta.ParseError("Pass failed!")
        end
    end
end # function parseall

tryparselast(::Type{T}, str::AbstractString) where {T<:AtomicStructure} =
    tryparseall(T, str)[end]
parselast(::Type{T}, str::AbstractString) where {T<:AtomicStructure} = parseall(T, str)[end]

function _parsenext_internal(
    ::Type{T}, str::AbstractString, start::Integer, raise::Bool
) where {T}
    x = findnext(REGEXOF[nameof(T)], str, start)
    if x === nothing
        raise ? throw(Meta.ParseError("Nothing found for next!")) : return nothing
    end
    return tryparse_internal(T, str[x])
end # function parsenext
tryparsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, false)
parsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, true)

function tryparsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && return nothing
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && return nothing
    return tryparse_internal(T, m.match)
end # function parsefinal
function parsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && throw(Meta.ParseError("No final coordinates found!"))
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && throw(Meta.ParseError("No `CELL_PARAMETERS` found!"))
    return tryparse_internal(T, m.match)
end # function parsefinal
