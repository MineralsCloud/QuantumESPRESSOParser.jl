using Dates: Hour, Minute, Millisecond
using DataFrames: AbstractDataFrame, DataFrame, groupby
using QuantumESPRESSOBase.PWscf
using VersionParsing: vparse

export IrreducibleBrillouinZone, TimedItem

struct SubroutineError
    name::String
    cerr::String
    msg::String
end

struct ParseError <: Exception
    msg::String
end

export Diagonalization,
    Preamble,
    Davidson,
    ConjugateGradient,
    ProjectedPreconditionedConjugateGradient,
    parse_symmetries,
    parse_stress,
    parse_bands,
    parse_all_electron_energy,
    parse_energy_decomposition,
    parse_paw_contribution,
    parse_smearing_energy,
    parse_version,
    parse_parallel_info,
    parse_fft_dimensions,
    parse_electrons_energies,
    parse_input_name,
    isoptimized,
    isjobdone,
    eachstep,
    eachiteration,
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

include("regexes.jl")

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported

abstract type PWOutputParameter end

Base.@kwdef struct Preamble <: PWOutputParameter
    ibrav::Int
    alat::Float64
    omega::Float64
    nat::Int
    ntyp::Int
    nelec::Float64
    nelup::Maybe{Float64} = nothing
    neldw::Maybe{Float64} = nothing
    nbnd::Int
    ecutwfc::Float64
    ecutrho::Float64
    ecutfock::Maybe{Float64} = nothing
    conv_thr::Maybe{Float64} = nothing
    mixing_beta::Maybe{Float64} = nothing
    mixing_ndim::Maybe{Int} = nothing
    mixing_mode::Maybe{String} = nothing
    xc::String
    nstep::Maybe{Int} = nothing
end

struct FFTGrid <: PWOutputParameter
    type::String
    dense::Int64
    smooth::Int64
    PW::Int64
end

"""
    parse_fft_base_info(str::AbstractString)

Parse the FFT base information from `pw.x`'s output and return a `DataFrame`.

If there are more than one processors, the title is "Parallelization info" and three
rows, i.e., "Min", "Max", and "Sum" are printed. If not, the title is
"G-vector sticks info" and only the "Sum" row is printed. If no information is found,
return `nothing`. The `DataFrame` is grouped by "sticks" and "gvecs".
"""
function Base.parse(::Type{FFTGrid}, str::AbstractString)
    obj = tryparse(FFTGrid, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{FFTGrid}, str::AbstractString)
    matched = match(FFT_BASE_INFO, str)
    if isnothing(matched)
        @info("The FFT base info is not found!")
        return nothing
    end
    body = matched[:body]
    data = FFTGrid[]
    for line in split(body, r"\R+")  # Don’t want empty lines
        # "Min",4X,2I8,I7,12X,2I9,I8
        splitted = split(line, " "; keepempty=false)  # Don't want empty strings
        values = map(Base.Fix1(parse, Int64), splitted[2:7])
        push!(data, FFTGrid("sticks", values[1:3]...))
        push!(data, FFTGrid("G-vecs", values[4:6]...))
    end
    return data
end

function parse_symmetries(str::AbstractString)
    m = match(SYM_OPS, str)
    m === nothing && return nothing
    return num_sym_ops = isempty(m[:n]) ? 0 : parse(Int, m[:n])
end # function parse_symmetries

struct IrreducibleBrillouinZone <: PWOutputParameter
    cartesian::Maybe{Vector{SpecialPoint}}
    crystal::Maybe{Vector{SpecialPoint}}
end

function Base.parse(::Type{IrreducibleBrillouinZone}, str::AbstractString)
    obj = tryparse(IrreducibleBrillouinZone, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{IrreducibleBrillouinZone}, str::AbstractString)
    matched = match(K_POINTS_BLOCK, str)
    if isnothing(matched)
        @info("The k-points info is not found!")
        return nothing
    end
    nk = parse(Int64, matched[:nk])
    cartesian, crystal = map((:cartesian, :crystal)) do key
        if !isnothing(matched[key])
            points = map(eachmatch(K_POINTS_ITEM, matched[key])) do matched
                SpecialPoint(map(Base.Fix1(parse, Float64), matched.captures[begin:end])...)
            end
            @assert length(points) == nk
            points
        end
    end
    return IrreducibleBrillouinZone(cartesian, crystal)
end

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

struct EachStep
    iterator::Base.RegexMatchIterator
end

function Base.iterate(iter::EachStep)
    iterated = iterate(iter.iterator)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return matched.match, state
    end
end
function Base.iterate(iter::EachStep, state)
    iterated = iterate(iter.iterator, state)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return matched.match, state
    end
end

Base.eltype(::Type{EachStep}) = String

Base.IteratorSize(::Type{EachStep}) = Base.SizeUnknown()

eachstep(str::AbstractString) = EachStep(eachmatch(SELF_CONSISTENT_CALCULATION_BLOCK, str))

struct EachIteration
    iterator::Base.RegexMatchIterator
end

function Base.iterate(iter::EachIteration)
    iterated = iterate(iter.iterator)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return matched.match, state
    end
end
function Base.iterate(iter::EachIteration, state)
    iterated = iterate(iter.iterator, state)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return matched.match, state
    end
end

Base.eltype(::Type{EachIteration}) = String

Base.IteratorSize(::Type{EachIteration}) = Base.SizeUnknown()

eachiteration(str::AbstractString) = EachIteration(eachmatch(ITERATION_BLOCK, str))

function _iterationwise!(f::Function, df::AbstractDataFrame, str::AbstractString)
    # Loop relax steps
    for (i, scf) in enumerate(eachmatch(SELF_CONSISTENT_CALCULATION_BLOCK, str))
        # Loop scf iterations
        for (j, iter) in enumerate(eachmatch(ITERATION_BLOCK, scf[1]))
            push!(df, [i j f(iter[1])...])
        end
    end
    return df
end # function _iterationwise

struct IterationHead <: PWOutputParameter
    number::Int64
    ecut::Float64
    beta::Float64
end

function Base.parse(::Type{IterationHead}, str::AbstractString)
    obj = tryparse(IterationHead, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{IterationHead}, str::AbstractString)
    matched = match(ITERATION_HEAD, str)
    if isnothing(matched)
        return nothing
    else
        return IterationHead(
            parse(Int64, matched[1]), parse(Float64, matched[2]), parse(Float64, matched[3])
        )
    end
end

struct IterationTime <: PWOutputParameter
    time::Float64
end

function Base.parse(::Type{IterationTime}, str::AbstractString)
    obj = tryparse(IterationTime, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{IterationTime}, str::AbstractString)
    matched = match(TOTAL_CPU_TIME, str)
    if isnothing(matched)
        return nothing
    else
        return IterationTime(parse(Float64, matched[1]))
    end
end

abstract type DiagonalizationSolver end
struct Davidson <: DiagonalizationSolver end
struct ConjugateGradient <: DiagonalizationSolver end
struct ProjectedPreconditionedConjugateGradient <: DiagonalizationSolver end

struct Diagonalization <: PWOutputParameter
    solver::DiagonalizationSolver
    ethr::Float64
    avg_iter::Float64
end

function Base.parse(::Type{Diagonalization}, str::AbstractString)
    obj = tryparse(Diagonalization, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{Diagonalization}, str::AbstractString)
    matched = match(C_BANDS, str)
    if !isnothing(matched)
        return nothing
    else
        solver = if matched[:diag] == "Davidson diagonalization with overlap"
            Davidson()
        elseif matched[:diag] == "CG style diagonalization"
            ConjugateGradient()
        elseif matched[:diag] == "PPCG style diagonalization"
            ProjectedPreconditionedConjugateGradient()
        else
            throw(ParseError("unknown diagonalization style!"))
        end
        ethr, avg_iter = map(Base.Fix1(parse, parse), matched.captures[2:end])
        return Diagonalization(solver, ethr, avg_iter)
    end
end

struct UnconvergedEnergy <: PWOutputParameter
    total_energy::Float64
    harris_foulkes_estimate::Maybe{Float64}
    estimated_scf_accuracy::Float64
end

function Base.parse(::Type{UnconvergedEnergy}, str::AbstractString)
    obj = tryparse(UnconvergedEnergy, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{UnconvergedEnergy}, str::AbstractString)
    matched = match(UNCONVERGED_ELECTRONS_ENERGY, str)
    if isnothing(matched)
        return nothing
    else
        ɛ, hf, δ = map(_parser, matched.captures)
        return UnconvergedEnergy(ɛ, hf, δ)
    end
end

_parser(x) = isnothing(x) ? x : parse(Float64, x)

struct ConvergedEnergy <: PWOutputParameter
    total_energy::Float64
    harris_foulkes_estimate::Maybe{Float64}
    estimated_scf_accuracy::Float64
end

function Base.parse(::Type{ConvergedEnergy}, str::AbstractString)
    obj = tryparse(ConvergedEnergy, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{ConvergedEnergy}, str::AbstractString)
    matched = match(CONVERGED_ELECTRONS_ENERGY, str)
    if isnothing(matched)
        return nothing
    else
        ɛ, hf, δ = map(_parser, matched.captures[1:3])
        return ConvergedEnergy(ɛ, hf, δ)
    end
end

function _parse_nonconverged_energy(str::AbstractString)
    ɛ, hf, δ = nothing, nothing, nothing  # Initialization
    m = match(UNCONVERGED_ELECTRONS_ENERGY, str)
    if m !== nothing
        ɛ, hf, δ = map(x -> x === nothing ? x : parse(Float64, x), m.captures)
    end  # Keep them `nothing` if `m` is `nothing`
    return ɛ, hf, δ
end # function _parse_nonconverged_energy
function _parse_electrons_energies(str::AbstractString, ::Val{:nonconverged})
    df = DataFrame(;
        step=Int[],
        iteration=Int[],
        ɛ=Maybe{Float64}[],  # Total energy
        hf=Maybe{Float64}[],  # Harris-Foulkes estimate
        δ=Maybe{Float64}[],  # Estimated scf accuracy
    )
    return _iterationwise!(_parse_nonconverged_energy, df, str)
end # function _parse_electrons_energies
function _parse_electrons_energies(str::AbstractString, ::Val{:converged})
    df = DataFrame(;
        step=Int[],
        ɛ=Maybe{Float64}[],  # Total energy
        hf=Maybe{Float64}[],  # Harris-Foulkes estimate
        δ=Maybe{Float64}[],  # Estimated scf accuracy
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if m !== nothing
            map(x -> x === nothing ? x : parse(Float64, x), m.captures[1:3])
        else
            ntuple(_ -> nothing, 3)
        end  # Keep them `nothing` if `m` is `nothing`
        push!(df, [i data...])
    end
    return df
end # function _parse_electrons_energies
function _parse_electrons_energies(str::AbstractString, ::Val{:combined})
    converged = parse_electrons_energies(str, :converged)
    nonconverged = parse_electrons_energies(str, :nonconverged)
    # TODO: Very ugly hack
    m = 1  # Initial step number
    for (i, n) in enumerate(nonconverged.step)
        if n != m
            @assert(all(==(nothing), nonconverged[i - 1, 3:5]))
            # nonconverged[i - 1, 3:5] = converged[n, 2:4]  # Converged energies do not have `iteration` column
            nonconverged[i - 1, 3] = converged[n, 2]
            nonconverged[i - 1, 4] = converged[n, 3]
            nonconverged[i - 1, 5] = converged[n, 4]
        end
        m = n  # Save the last step number
    end
    if m == 1
        nonconverged[end, 3] = converged[end, 2]
        nonconverged[end, 4] = converged[end, 3]
        nonconverged[end, 5] = converged[end, 4]
    end
    return nonconverged
end # function _parse_electrons_energies
function parse_electrons_energies(str::AbstractString, option::Symbol)
    @assert(option ∈ (:combined, :converged, :nonconverged))
    return _parse_electrons_energies(str, Val(option))
end # function parse_electrons_energies

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

struct TimedItem <: PWOutputParameter
    name::String
    cpu::Millisecond
    wall::Millisecond
    calls::Maybe{Int64}
end

function Base.parse(::Type{TimedItem}, str::AbstractString)
    obj = tryparse(TimedItem, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
function Base.tryparse(::Type{TimedItem}, str::AbstractString)
    matched = match(TIMED_ITEM, str)
    if isnothing(matched)
        return nothing
    else
        name, cpu, wall = matched[1], parsetime(matched[2]), parsetime(matched[9])
        return TimedItem(
            name, cpu, wall, isnothing(matched[16]) ? nothing : parse(Int64, matched[16])
        )
    end
end

function parsetime(str::AbstractString)
    matched = match((HOURS_MINUTES), str)
    if !isnothing(matched)
        hours = parse(Int64, matched[1])
        minutes = parse(Int64, matched[2])
        return convert(Millisecond, Hour(hours) + Minute(minutes))
    end
    matched = match((MINUTES_SECONDS), str)
    if !isnothing(matched)
        minutes = parse(Int64, matched[1])
        seconds = parse(Float64, matched[2])
        return convert(Millisecond, Minute(minutes)) +
               Millisecond(round(Int64, 1000seconds))
    end
    matched = match((SECONDS), str)
    if !isnothing(matched)
        seconds = parse(Float64, matched[1])
        return Millisecond(round(Int64, 1000seconds))  # 1000 times a floating point number may not be an integer
    end
    throw(ParseError("unrecognized time format!"))
end

function parse_input_name(str::AbstractString)
    m = match(READING_INPUT_FROM, str)
    return m === nothing ? nothing : only(m)
end # function parse_input_name

isoptimized(str::AbstractString) =
    match(FINAL_COORDINATES_BLOCK, str) === nothing ? false : true

isjobdone(str::AbstractString) = match(JOB_DONE, str) !== nothing

# This is an internal function and should not be exported.
function Base.tryparse(::Type{Preamble}, str::AbstractString)
    dict = Dict{Symbol,Any}()
    m = match(SUMMARY_BLOCK, str)
    return if m !== nothing
        body = only(m.captures)
        f = (T, x) -> T == String ? string(x) : parse(T, x)
        for (field, regex, T) in zip(
            (
                :ibrav,
                :alat,
                :omega,
                :nat,
                :ntyp,
                :nelec,
                :nbnd,
                :ecutwfc,
                :ecutrho,
                :ecutfock,
                :conv_thr,
                :mixing_beta,
                :mixing_ndim,
                :xc,
                :nstep,
            ),
            (
                NUMBER_OF_ATOMS_PER_CELL,
                LATTICE_PARAMETER,
                UNIT_CELL_VOLUME,
                NUMBER_OF_ATOMS_PER_CELL,
                NUMBER_OF_ATOMIC_TYPES,
                NUMBER_OF_ELECTRONS,
                NUMBER_OF_KOHN_SHAM_STATES,
                KINETIC_ENERGY_CUTOFF,
                CHARGE_DENSITY_CUTOFF,
                CUTOFF_FOR_FOCK_OPERATOR,
                CONVERGENCE_THRESHOLD,
                MIXING_BETA,
                NUMBER_OF_ITERATIONS_USED,
                EXCHANGE_CORRELATION,
                NSTEP,
            ),
            (
                Int,
                Float64,
                Float64,
                Int,
                Int,
                Float64,
                Int,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Int,
                String,
                Int,
            ),
        )
            m = match(regex, body)
            if m !== nothing
                dict[field] = f(T, m[1])
            end
        end
        # 2 special cases
        let x = match(NUMBER_OF_ELECTRONS, body), y = match(NUMBER_OF_ITERATIONS_USED, body)
            if all(!=(nothing), x.captures[2:end])
                dict[:nelup], dict[:neldw] = parse.(Float64, x.captures[2:end])
            end
            if y !== nothing
                dict[:mixing_mode] = y[2]
            end
        end
        Preamble(; dict...)
    end
end # function Base.tryparse
function Base.tryparse(::Type{SubroutineError}, str::AbstractString)
    # According to my observation, a QE output can have at most one type of
    # `SubroutineError`. Warn me if there can be multiple types of errors.
    m = match(ERROR_BLOCK, str)
    return if m !== nothing
        # `tryparse` returns nothing if the string does not contain what we want,
        # while `parse` raises an error.
        body = strip(m[:body])
        # Referenced from https://stackoverflow.com/a/454919/3260253
        e, msg = map(strip, split(body, r"\R+"))
        m = match(ERROR_IN_ROUTINE, e)
        SubroutineError(m[1], m[2], msg)
    end
end # function Base.tryparse

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

function Base.parse(
    ::Type{T}, str::AbstractString
) where {T<:Union{Preamble,SubroutineError}}
    x = tryparse(T, str)
    return x === nothing ? throw(Meta.ParseError("cannot find `$(T)`!")) : x
end # function Base.parse

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
