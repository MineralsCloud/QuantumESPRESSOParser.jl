"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: only, @NamedTuple
# using Dates: DateTime, DateFormat
using DataFrames: AbstractDataFrame, DataFrame, groupby
using Parameters: @with_kw
using QuantumESPRESSOBase.Inputs.PWscf
using ReadableRegex:
    DEC_DIGIT_NUMBER,
    NON_SEPARATOR,
    WHITESPACE,
    maybe,
    either,
    char_in,
    one_or_more,
    zero_or_more,
    exactly,
    look_for,
    capture,
    @rs_str
using VersionParsing: vparse

using ..Outputs: SubroutineError

export Diagonalization,
    Preamble,
    FftDimensions,
    IrreducibleBrillouinZone,
    Davidson,
    ParallelInfo,
    ConjugateGradient,
    ProjectedPreconditionedConjugateGradient,
    parse_fft_base_info,
    parse_symmetries,
    parse_stress,
    parse_iteration_time,
    parse_bands,
    parse_all_electron_energy,
    parse_energy_decomposition,
    parse_paw_contribution,
    parse_smearing_energy,
    parse_version,
    parse_iteration_head,
    parse_electrons_energies,
    parse_clock,
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

include("regexes.jl")

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported

abstract type Diagonalization end
struct Davidson <: Diagonalization end
struct ConjugateGradient <: Diagonalization end
struct ProjectedPreconditionedConjugateGradient <: Diagonalization end

const FftDimensions = @NamedTuple begin
    ng::UInt
    nr1::UInt
    nr2::UInt
    nr3::UInt
end

const IrreducibleBrillouinZone = @NamedTuple begin
    cart::Maybe{SpecialPointsCard}
    cryst::Maybe{SpecialPointsCard}
end

const ParallelInfo = @NamedTuple begin
    flavor::String
    np::UInt
end

@with_kw struct Preamble
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
    xc::Maybe{String} = nothing  # Now is deprecated
    nstep::Maybe{Int} = nothing
end

"""
    parse_fft_base_info(str::AbstractString)

Parse the FFT base information from `pw.x`'s output and return a `DataFrame`.

If there are more than one processors, the title is "Parallelization info" and three
rows, i.e., "Min", "Max", and "Sum" are printed. If not, the title is
"G-vector sticks info" and only the "Sum" row is printed. If no information is found,
return `nothing`. The `DataFrame` is grouped by "sticks" and "gvecs".
"""
function parse_fft_base_info(str::AbstractString)::Maybe{AbstractDataFrame}
    df =
        DataFrame(kind = String[], stats = String[], dense = Int[], smooth = Int[], PW = [])
    m = match(FFT_BASE_INFO, str)
    if m === nothing
        @info("The FFT base info is not found!")
        return
    end
    body = m[:body]
    for line in split(body, r"\R+")  # Don’t want empty lines
        # "Min",4X,2I8,I7,12X,2I9,I8
        sp = split(line, " ", keepempty = false)  # Don't want empty strings
        numbers = map(x -> parse(Int, x), sp[2:7])
        push!(df, ["sticks" sp[1] numbers[1:3]...])
        push!(df, ["gvecs" sp[1] numbers[4:6]...])
    end
    return df
end # function parse_fft_base_info

function parse_symmetries(str::AbstractString)
    m = match(SYM_OPS, str)
    m === nothing && return
    num_sym_ops = isempty(m[:n]) ? 0 : parse(Int, m[:n])
end # function parse_symmetries

# Return `nothing`, `(cartesian_coordinates, nothing)`, `(nothing, crystal_coordinates)`, `(cartesian_coordinates, crystal_coordinates)`
function Base.tryparse(::Type{IrreducibleBrillouinZone}, str::AbstractString)
    block = match(K_POINTS_BLOCK, str)
    if block === nothing
        @info("The k-points info is not found!")
        return
    else
        result = map((:cart => "tpiba", :cryst => "crystal")) do (k, v)
            if block[k] !== nothing
                data = map(eachmatch(K_POINTS_ITEM, block[k])) do row
                    map(
                        Base.Fix1(parse, Float64),
                        (
                            split(row.captures[1], ' '; keepempty = false)...,
                            row.captures[2],
                        ),
                    )
                end
                @assert length(data) == parse(Int, block[:nk])
                SpecialPointsCard(data, v)
            else
                nothing
            end
        end
        return IrreducibleBrillouinZone(result)
    end
end # function parse_ibz

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses, kbar_stresses = Matrix{Float64}[], Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK, str)
        pressure, content = m.captures[1], m.captures[3]
        push!(pressures, parse(Float64, pressure))

        stress_atomic, stress_kbar = ntuple(_ -> Matrix{Float64}(undef, 3, 3), 2)
        for (i, line) in enumerate(split(content, r"\R+"))
            tmp = map(x -> parse(Float64, x), split(line, " ", keepempty = false))
            stress_atomic[i, :], stress_kbar[i, :] = tmp[1:3], tmp[4:6]
        end
        push!(atomic_stresses, stress_atomic)
        push!(kbar_stresses, stress_kbar)
    end
    return pressures, atomic_stresses, kbar_stresses
end # function parse_stress

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

function parse_iteration_head(str::AbstractString)
    df = DataFrame(
        step = Int[],  # Step number
        iteration = Int[],  # Iteration number
        ecut = Float64[],  # Cutoff energy
        β = Float64[],  # Mixing beta
    )
    return _iterationwise!(_parse_iteration_head, df, str)
end # function parse_iteration_head
# This is a helper function and should not be exported.
function _parse_iteration_head(str::AbstractString)
    head = match(ITERATION_HEAD, str)
    return map(x -> parse(Float64, x), head.captures[2:3])
end # function _parse_iteration_head

function parse_iteration_time(str::AbstractString)
    df = DataFrame(step = Int[], iteration = Int[], time = Float64[])
    return _iterationwise!(_parse_iteration_time, df, str)
end # function parse_iteration_time
# This is a helper function and should not be exported.
function _parse_iteration_time(str::AbstractString)
    return parse(Float64, match(TOTAL_CPU_TIME, str)[1])
end # function _parse_iteration_time

function parse_diagonalization(str::AbstractString)
    df = DataFrame(
        step = Int[],
        iteration = Int[],
        diag = Diagonalization[],  # Diagonalization style
        ethr = Float64[],  # Energy threshold
        avg = Float64[],  # Average # of iterations
    )
    return _iterationwise!(_parse_diagonalization, df, str)
end # function parse_diagonalization
# This is a helper function and should not be exported.
function _parse_diagonalization(str::AbstractString)
    solver, ethr, avg_iter = nothing, nothing, nothing  # Initialization
    m = match(C_BANDS, str)
    if m !== nothing
        solver = if m[:diag] == "Davidson diagonalization with overlap"
            Davidson()
        elseif m[:diag] == "CG style diagonalization"
            ConjugateGradient()
        elseif m[:diag] == "PPCG style diagonalization"
            ProjectedPreconditionedConjugateGradient()
        else
            error("unknown diagonalization style!")
        end
        ethr, avg_iter = map(x -> parse(Float64, x), m.captures[2:end])
    end  # Keep them `nothing` if `m` is `nothing`
    return solver, ethr, avg_iter
end # function _parse_diagonalization

function _parse_nonconverged_energy(str::AbstractString)
    ɛ, hf, δ = nothing, nothing, nothing  # Initialization
    m = match(UNCONVERGED_ELECTRONS_ENERGY, str)
    if m !== nothing
        ɛ, hf, δ = map(x -> parse(Float64, x), m.captures)
    end  # Keep them `nothing` if `m` is `nothing`
    return ɛ, hf, δ
end # function _parse_nonconverged_energy
function _parse_electrons_energies(str::AbstractString, ::Val{:nonconverged})
    df = DataFrame(
        step = Int[],
        iteration = Int[],
        ɛ = Maybe{Float64}[],  # Total energy
        hf = Maybe{Float64}[],  # Harris-Foulkes estimate
        δ = Maybe{Float64}[],  # Estimated scf accuracy
    )
    return _iterationwise!(_parse_nonconverged_energy, df, str)
end # function _parse_electrons_energies
function _parse_electrons_energies(str::AbstractString, ::Val{:converged})
    df = DataFrame(
        step = Int[],
        ɛ = Maybe{Float64}[],  # Total energy
        hf = Maybe{Float64}[],  # Harris-Foulkes estimate
        δ = Maybe{Float64}[],  # Estimated scf accuracy
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if m !== nothing
            map(x -> parse(Float64, x), m.captures[1:3])
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
            @assert(all(==(nothing), nonconverged[i-1, 3:5]))
            # nonconverged[i - 1, 3:5] = converged[n, 2:4]  # Converged energies do not have `iteration` column
            nonconverged[i-1, 3] = converged[n, 2]
            nonconverged[i-1, 4] = converged[n, 3]
            nonconverged[i-1, 5] = converged[n, 4]
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
    str == "Number of k-points >= 100: set verbosity='high' to print the bands." && return
    kpts, bands = nothing, nothing  # Initialization
    m = match(KS_ENERGIES_BLOCK, str)
    if m !== nothing
        kpts, bands = Vector{Float64}[], Vector{Float64}[]
        regex =
            match(KS_ENERGIES_BANDS, str) === nothing ? KS_ENERGIES_BAND_ENERGIES :
            KS_ENERGIES_BANDS
        for m in eachmatch(regex, str)
            push!(kpts, map(x -> parse(Float64, x[1]), eachmatch(Regex(EXP_REAL), m[:k])))
            push!(
                bands,
                map(x -> parse(Float64, x[1]), eachmatch(Regex(EXP_REAL), m[:band])),
            )
        end
        len, nbnd = length(kpts), length(bands[1])
        kpts, bands = reshape(Iterators.flatten(kpts) |> collect, len, 3),
        reshape(Iterators.flatten(bands) |> collect, len, nbnd)
    end  # Keep them `nothing` if `m` is `nothing`
    return kpts, bands
end # function parse_bands

function parse_all_electron_energy(str::AbstractString)
    df = DataFrame(step = Int[], ae = Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        ae = if any(==(nothing), (m, m[:ae]))
            nothing
        else
            parse(Float64, match(Regex(REAL), m[:ae])[1])
        end
        push!(df, [i ae])
    end
    return df
end # function parse_all_electron_energy

function parse_energy_decomposition(str::AbstractString)
    df = DataFrame(
        step = Int[],
        onelectron = Maybe{Float64}[],
        hartree = Maybe{Float64}[],
        xc = Maybe{Float64}[],
        ewald = Maybe{Float64}[],
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if any(==(nothing), (m, m[:decomp]))
            ntuple(_ -> nothing, 4)
        else
            map(x -> parse(Float64, x[1]), eachmatch(Regex(REAL), m[:decomp]))
        end
        push!(df, [i data...])
    end
    return df
end # function parse_energy_decomposition

function parse_paw_contribution(str::AbstractString)
    df = DataFrame(
        step = Int[],
        hartree_ae = Maybe{Float64}[],
        hartree_ps = Maybe{Float64}[],
        xc_ae = Maybe{Float64}[],
        xc_ps = Maybe{Float64}[],
        eh = Maybe{Float64}[],
        exc = Maybe{Float64}[],
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if any(==(nothing), (m, m[:one]))
            ntuple(_ -> nothing, 6)
        else
            map(x -> parse(Float64, x[1]), eachmatch(Regex(REAL), m[:one]))
        end
        push!(df, [i data...])
    end
    return df
end # function parse_paw_contribution

function parse_smearing_energy(str::AbstractString)
    df = DataFrame(step = Int[], smearing = Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        smearing = if any(==(nothing), (m, m[:smearing]))
            nothing
        else
            parse(Float64, match(Regex(REAL), m[:smearing])[1])
        end
        push!(df, [i smearing])
    end
    return df
end # function parse_smearing_energy

function parse_version(str::AbstractString)::Maybe{VersionNumber}
    m = match(PWSCF_VERSION, str)
    m !== nothing ? vparse(m[:version]) : return
end # function parse_version

function Base.tryparse(::Type{FftDimensions}, str::AbstractString)
    m = match(FFT_DIMENSIONS, str)
    if m === nothing
        return
    else
        ng = parse(Int, m.captures[1])
        nr = map(x -> parse(Int, x), split(m.captures[2], ','))
        return FftDimensions((ng, nr...))
    end
end # function parse_fft_dimensions

function parse_clock(str::AbstractString)::Maybe{AbstractDataFrame}
    m = match(TIME_BLOCK, str)
    m === nothing && return
    content = only(m.captures)

    info = DataFrame(
        subroutine = String[],
        item = String[],
        CPU = Float64[],
        wall = Float64[],
        calls = Int[],
    )
    for regex in [
        SUMMARY_TIME_BLOCK
        INIT_RUN_TIME_BLOCK
        ELECTRONS_TIME_BLOCK
        C_BANDS_TIME_BLOCK
        SUM_BAND_TIME_BLOCK
        EGTERG_TIME_BLOCK
        H_PSI_TIME_BLOCK
        GENERAL_ROUTINES_TIME_BLOCK
        PARALLEL_ROUTINES_TIME_BLOCK
    ]
        block = match(regex, content)
        if block !== nothing
            for m in eachmatch(TIME_ITEM, block[:body])
                push!(
                    info,
                    [block[:head] m[1] map(x -> parse(Float64, x), m.captures[2:4])...],
                )
            end
        end
    end
    # m = match(TERMINATED_DATE, content)
    # info["terminated date"] = parse(DateTime, m.captures[1], DateFormat("H:M:S"))
    return info
end # function parse_clock

function parse_input_name(str::AbstractString)
    m = match(READING_INPUT_FROM, str)
    return m === nothing ? nothing : only(m)
end # function parse_input_name

isoptimized(str::AbstractString) =
    match(FINAL_COORDINATES_BLOCK, str) === nothing ? false : true

isjobdone(str::AbstractString) = match(JOB_DONE, str) !== nothing

# This is an internal function and should not be exported.
function Base.tryparse(::Type{ParallelInfo}, str::AbstractString)
    m = match(PARALLEL_INFO, str)
    if m === nothing
        return
    else
        np = m.captures[2] === nothing ? 1 : parse(Int, m.captures[2])
        return ParallelInfo((m.captures[1], np))
    end
end
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
    m = match(CELL_PARAMETERS_BLOCK, str)
    return if m !== nothing
        body, data = m[:data], Matrix{Float64}(undef, 3, 3)  # Initialization
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM, body))
            data[i, :] = map(x -> parse(Float64, x), matched.captures)
        end
        if m[:option] == "alat"
            alat = parse(Float64, m[:alat])
            CellParametersCard(alat * data, "bohr")
        else
            CellParametersCard(data, m[:option])
        end
    end
end # function tryparse_internal
function tryparse_internal(::Type{AtomicPositionsCard}, str::AbstractString)
    atomic_positions = AtomicPositionsCard[]
    m = match(ATOMIC_POSITIONS_BLOCK, str)
    return if m !== nothing
        option = string(m[1])
        body = m[2]
        data = AtomicPosition[]
        for matched in eachmatch(ATOMIC_POSITIONS_ITEM, body)
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
    ::Type{T},
    str::AbstractString,
) where {T<:Union{Preamble,SubroutineError}}
    x = tryparse(T, str)
    x === nothing ? throw(Meta.ParseError("cannot find `$(T)`!")) : x
end # function Base.parse

function _parse(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    x = tryparse(T, str)
    x === nothing ? throw(Meta.ParseError("cannot find `$(T)`!")) : x
end # function _parse

const REGEXOF = (
    CellParametersCard = CELL_PARAMETERS_BLOCK,
    AtomicPositionsCard = ATOMIC_POSITIONS_BLOCK,
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
    return map(eachmatch(REGEXOF(T), str)) do x
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
    ::Type{T},
    str::AbstractString,
    start::Integer,
    raise::Bool,
) where {T}
    x = findnext(REGEXOF(T), str, start)
    if x === nothing
        raise ? throw(Meta.ParseError("Nothing found for next!")) : return
    end
    return tryparse_internal(T, str[x])
end # function parsenext
tryparsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, false)
parsenext(::Type{T}, str::AbstractString, start::Integer) where {T} =
    _parsenext_internal(T, str, start, true)

function tryparsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && return
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && return
    return tryparse_internal(T, m.match)
end # function parsefinal
function parsefinal(::Type{T}, str::AbstractString) where {T<:AtomicStructure}
    m = match(FINAL_COORDINATES_BLOCK, str)
    m === nothing && throw(Meta.ParseError("No final coordinates found!"))
    m = match(REGEXOF[nameof(T)], m.match)
    m === nothing && throw(Meta.ParseError("No `CELL_PARAMETERS` found!"))
    return tryparse_internal(T, m.match)
end # function parsefinal

end
