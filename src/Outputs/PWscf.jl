"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: isnothing
# using Dates: DateTime, DateFormat
using DataFrames: AbstractDataFrame, DataFrame, groupby
using Fortran90Namelists.FortranToJulia
using Rematch: @match
using Parameters: @with_kw
using QuantumESPRESSOBase.Cards.PWscf

using QuantumESPRESSOParsers

export DiagonalizationStyle,
       Preamble,
       DavidsonDiagonalization,
       CGDiagonalization,
       PPCGDiagonalization,
       parse_fft_base_info,
       parse_symmetries,
       parse_ibz,
       parse_stress,
       parse_iteration_time,
       parse_unconverged_energy,
       parse_bands,
       parse_converged_energy,
       parse_all_electron_energy,
       parse_energy_decomposition,
       parse_paw_contribution,
       parse_smearing_energy,
       parse_version,
       parse_parallel_info,
       parse_fft_dimensions,
       parse_atomic_positions,
       parse_iteration_head,
       parse_clock,
       whatinput,
       isrelaxed,
       isjobdone,
       tryparsefirst,
       parsefirst,
       tryparseall,
       parseall,
       tryparselast,
       parselast

include("regexes.jl")

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported

abstract type DiagonalizationStyle end
struct DavidsonDiagonalization <: DiagonalizationStyle end
struct CGDiagonalization <: DiagonalizationStyle end
struct PPCGDiagonalization <: DiagonalizationStyle end

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
    xc::String
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
    df = DataFrame(
        kind = String[],
        stats = String[],
        dense = Int[],
        smooth = Int[],
        PW = [],
    )
    m = match(FFT_BASE_INFO, str)
    if isnothing(m)
        @info("The FFT base info is not found!")
        return
    end
    body = m[:body]
    for line in split(body, r"[\r\n]+")  # Don’t want empty lines
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
    isnothing(m) && return
    num_sym_ops = isempty(m[:n]) ? 0 : parse(Int, m[:n])
end # function parse_symmetries

# Return `nothing`, `(cartesian_coordinates, nothing)`, `(nothing, crystal_coordinates)`, `(cartesian_coordinates, crystal_coordinates)`
function parse_ibz(str::AbstractString)::Maybe{Tuple}
    m = match(K_POINTS_BLOCK, str)
    if isnothing(m)
        @info("The k-points info is not found!")
        return
    end
    nk = parse(Int, m[:nk])
    result = []
    kinds = (:cart => "tpiba", :cryst => "crystal")
    for (k, v) in kinds
        if !isnothing(m[k])
            x = Matrix{Float64}(undef, nk, 4)
            for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, m[k]))
                x[i, :] = map(x -> parse(Float64, x), m.captures[1:end])
            end
            push!(result, KPointsCard(v, x))
        else
            push!(result, nothing)
        end
    end
    return Tuple(result)
end # function parse_ibz

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses, kbar_stresses = Matrix{Float64}[], Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK, str)
        pressure, content = m.captures[1], m.captures[3]
        push!(pressures, parse(Float64, pressure))

        stress_atomic, stress_kbar = ntuple(_ -> Matrix{Float64}(undef, 3, 3), 2)
        for (i, line) in enumerate(split(content, r"[\r\n]+"))
            tmp = map(x -> parse(Float64, x), split(line, " ", keepempty = false))
            stress_atomic[i, :], stress_kbar[i, :] = tmp[1:3], tmp[4:6]
        end
        push!(atomic_stresses, stress_atomic)
        push!(kbar_stresses, stress_kbar)
    end
    return pressures, atomic_stresses, kbar_stresses
end # function parse_stress

function parse_atomic_positions(str::AbstractString)::Vector{<:AtomicPositionsCard}
    atomic_positions = AtomicPositionsCard[]
    for m in eachmatch(ATOMIC_POSITIONS_BLOCK, str)
        unit = string(m.captures[1])
        content = m.captures[2]
        data = AtomicPosition[]

        for matched in eachmatch(ATOMIC_POSITIONS_ITEM, content)
            captured = matched.captures
            if_pos = map(x -> isempty(x) ? 1 : parse(Int, FortranData(x)), captured[11:13])
            atom, pos = string(captured[1]),
                map(
                    x -> parse(Float64, FortranData(x)),
                    [captured[3], captured[6], captured[9]],
                )
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        push!(atomic_positions, AtomicPositionsCard(unit, data))
    end
    return atomic_positions
end # parse_atomic_positions

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
        diag = DiagonalizationStyle[],  # Diagonalization style
        ethr = Float64[],  # Energy threshold
        avg = Float64[],  # Average # of iterations
    )
    return _iterationwise!(_parse_diagonalization, df, str)
end # function parse_diagonalization
# This is a helper function and should not be exported.
function _parse_diagonalization(str::AbstractString)
    solver, ethr, avg_iter = nothing, nothing, nothing  # Initialization
    m = match(C_BANDS, str)
    if !isnothing(m)
        solver = @match m[:diag] begin
            "Davidson diagonalization with overlap" => DavidsonDiagonalization()
            "CG style diagonalization" => CGDiagonalization()
            "PPCG style diagonalization" => PPCGDiagonalization()
        end
        ethr, avg_iter = map(x -> parse(Float64, x), m.captures[2:end])
    end  # Keep them `nothing` if `m` is `nothing`
    return solver, ethr, avg_iter
end # function _parse_diagonalization

function parse_unconverged_energy(str::AbstractString)
    df = DataFrame(
        step = Int[],
        iteration = Int[],
        ɛ = Maybe{Float64}[],  # Total energy
        hf = Maybe{Float64}[],  # Harris-Foulkes estimate
        δ = Maybe{Float64}[],  # Estimated scf accuracy
    )
    return _iterationwise!(_parse_unconverged_energy, df, str)
end # function parse_unconverged_energy
# This is a helper function and should not be exported.
function _parse_unconverged_energy(str::AbstractString)
    ɛ, hf, δ = nothing, nothing, nothing  # Initialization
    m = match(UNCONVERGED_ELECTRONS_ENERGY, str)
    if !isnothing(m)
        ɛ, hf, δ = map(x -> parse(Float64, x), m.captures)
    end  # Keep them `nothing` if `m` is `nothing`
    return ɛ, hf, δ
end # function _parse_unconverged_energy

# See https://github.com/QEF/q-e/blob/4132a64/PW/src/print_ks_energies.f90#L10.
function parse_bands(str::AbstractString)
    str == "Number of k-points >= 100: set verbosity='high' to print the bands." && return
    kpts, bands = nothing, nothing  # Initialization
    m = match(KS_ENERGIES_BLOCK, str)
    if !isnothing(m)
        kpts, bands = Vector{Float64}[], Vector{Float64}[]
        regex = isnothing(match(KS_ENERGIES_BANDS, str)) ? KS_ENERGIES_BAND_ENERGIES :
                KS_ENERGIES_BANDS
        for m in eachmatch(regex, str)
            push!(
                kpts,
                map(x -> parse(Float64, x[1]), eachmatch(Regex(GENERAL_REAL), m[:k])),
            )
            push!(
                bands,
                map(x -> parse(Float64, x[1]), eachmatch(Regex(GENERAL_REAL), m[:band])),
            )
        end
        len, nbnd = length(kpts), length(bands[1])
        kpts, bands = reshape(Iterators.flatten(kpts) |> collect, len, 3),
            reshape(Iterators.flatten(bands) |> collect, len, nbnd)
    end  # Keep them `nothing` if `m` is `nothing`
    return kpts, bands
end # function parse_bands

function parse_converged_energy(str::AbstractString)
    df = DataFrame(
        step = Int[],
        ɛ = Maybe{Float64}[],  # Total energy
        hf = Maybe{Float64}[],  # Harris-Foulkes estimate
        δ = Maybe{Float64}[],  # Estimated scf accuracy
    )
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        data = if !isnothing(m)
            map(x -> parse(Float64, x), m.captures[1:3])
        else
            ntuple(_ -> nothing, 3)
        end  # Keep them `nothing` if `m` is `nothing`
        push!(df, [i data...])
    end
    return df
end # function parse_converged_energy

function parse_all_electron_energy(str::AbstractString)
    df = DataFrame(step = Int[], ae = Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        ae = if any(isnothing, (m, m[:ae]))
            nothing
        else
            parse(Float64, match(Regex(FIXED_POINT_REAL), m[:ae])[1])
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
        data = if any(isnothing, (m, m[:decomp]))
            ntuple(_ -> nothing, 4)
        else
            map(
                x -> parse(Float64, x[1]),
                eachmatch(Regex(FIXED_POINT_REAL), m[:decomp]),
            )
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
        data = if any(isnothing, (m, m[:one]))
            ntuple(_ -> nothing, 6)
        else
            map(x -> parse(Float64, x[1]), eachmatch(Regex(FIXED_POINT_REAL), m[:one]))
        end
        push!(df, [i data...])
    end
    return df
end # function parse_paw_contribution

function parse_smearing_energy(str::AbstractString)
    df = DataFrame(step = Int[], smearing = Maybe{Float64}[])
    for (i, m) in enumerate(eachmatch(CONVERGED_ELECTRONS_ENERGY, str))
        smearing = if any(isnothing, (m, m[:smearing]))
            nothing
        else
            parse(Float64, match(Regex(FIXED_POINT_REAL), m[:smearing])[1])
        end
        push!(df, [i smearing])
    end
    return df
end # function parse_smearing_energy

function parse_version(str::AbstractString)::Maybe{String}
    m = match(PWSCF_VERSION, str)
    !isnothing(m) ? m[:version] : return
end # function parse_version

function parse_parallel_info(str::AbstractString)::Maybe{Tuple{String,Int}}
    m = match(PARALLEL_INFO, str)
    isnothing(m) && return
    return m[:kind], isnothing(m[:num]) ? 1 : parse(Int, m[:num])
end # function parse_parallel_info

function parse_fft_dimensions(str::AbstractString)::Maybe{Tuple{Int,NamedTuple}}
    m = match(FFT_DIMENSIONS, str)
    isnothing(m) && return
    parsed = map(x -> parse(Int, x), m.captures)
    return parsed[1], NamedTuple{(:nr1, :nr2, :nr3)}(parsed[2:end])
end # function parse_fft_dimensions

function parse_clock(str::AbstractString)::Maybe{AbstractDataFrame}
    m = match(TIME_BLOCK, str)
    isnothing(m) && return
    content = m.captures[1]

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
        isnothing(block) && continue
        for m in eachmatch(TIME_ITEM, block[:body])
            push!(info, [block[:head] m[1] map(x -> parse(Float64, x), m.captures[2:4])...])
        end
    end
    # m = match(TERMINATED_DATE, content)
    # info["terminated date"] = parse(DateTime, m.captures[1], DateFormat("H:M:S"))
    return info
end # function parse_clock

function whatinput(str::AbstractString)::Maybe{String}
    m = match(READING_INPUT_FROM, str)
    !isnothing(m) ? m[1] : return
end # function whatinput

function isrelaxed(str::AbstractString)::Bool
    isnothing(match(FINAL_COORDINATES_BLOCK, str)) ? false : true
end # function isrelaxed

isjobdone(str::AbstractString) = !isnothing(match(JOB_DONE, str))

# This is an internal function and should not be exported.
function tryparse_internal(::Type{T}, str::AbstractString, raise::Bool) where {T<:Preamble}
    arr = Pair{Symbol,Any}[]
    m = match(SUMMARY_BLOCK, str)
    if isnothing(m)
        raise ? throw(Meta.ParseError("No info found!")) : return
    end
    body = first(m.captures)
    for (field, regex) in (
        :ibrav => NUMBER_OF_ATOMS_PER_CELL,
        :alat => LATTICE_PARAMETER,
        :omega => UNIT_CELL_VOLUME,
        :nat => NUMBER_OF_ATOMS_PER_CELL,
        :ntyp => NUMBER_OF_ATOMIC_TYPES,
        :nelec => NUMBER_OF_ELECTRONS,
        :nbnd => NUMBER_OF_KOHN_SHAM_STATES,
        :ecutwfc => KINETIC_ENERGY_CUTOFF,
        :ecutrho => CHARGE_DENSITY_CUTOFF,
        :ecutfock => CUTOFF_FOR_FOCK_OPERATOR,
        :conv_thr => CONVERGENCE_THRESHOLD,
        :mixing_beta => MIXING_BETA,
        :mixing_ndim => NUMBER_OF_ITERATIONS_USED,
        :xc => EXCHANGE_CORRELATION,
        :nstep => NSTEP,
    )
        m = match(regex, body)
        if !isnothing(m)
            S = QuantumESPRESSOParsers.nonnothingtype(fieldtype(Preamble, field))
            push!(arr, field => (S <: AbstractString ? string : Base.Fix1(parse, S))(m[1]))
        end
    end
    # 2 special cases
    m = match(NUMBER_OF_ELECTRONS, body)
    if all(!isnothing, m.captures[2:end])
        push!(arr, zip([:nelup, :neldw], map(x -> parse(Float64, x), m.captures[2:end])))
    end
    m = match(NUMBER_OF_ITERATIONS_USED, body)
    if !isnothing(m)
        push!(arr, :mixing_mode => m[2])
    end
    return T(; arr...)
end # function tryparse_internal
function tryparse_internal(
    ::Type{T},
    str::AbstractString,
    raise::Bool,
) where {T<:SubroutineError}
    # According to my observation, a QE output can have at most one type of
    # `SubroutineError`. Warn me if there can be multiple types of errors.
    m = match(ERROR_BLOCK, str)
    if isnothing(m)
        # `tryparse` returns nothing if the string does not contain what we want,
        # while `parse` raises an error.
        raise ? throw(Meta.ParseError("No error found!")) : return
    end
    body = strip(m[:body])
    # Referenced from https://stackoverflow.com/a/454919/3260253
    e, msg = map(strip, split(body, r"[\r\n]+"))
    m = match(ERROR_IN_ROUTINE, e)
    return T(m[1], m[2], msg)
end # function tryparse_internal
function tryparse_internal(
    ::Type{T},
    str::AbstractString,
    raise::Bool,
) where {T<:CellParametersCard}
    m = match(CELL_PARAMETERS_BLOCK, str)
    if isnothing(m)
        raise ? throw(Meta.ParseError("Cannot find `CELL_PARAMETERS`!")) : return
    end
    alat = parse(Float64, m[:alat])
    body = m[:data]

    data = Matrix{Float64}(undef, 3, 3)
    for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM, body))
        data[i, :] = map(x -> parse(Float64, x), matched.captures)
    end
    return T("bohr", alat * data)
end # function tryparse_internal

const _INTERNAL_TYPES = Union{Preamble,SubroutineError,CellParametersCard}

function Base.tryparse(::Type{T}, str::AbstractString) where {T<:_INTERNAL_TYPES}
    return tryparse_internal(T, str, false)
end # function Base.tryparse
function Base.parse(::Type{T}, str::AbstractString) where {T<:_INTERNAL_TYPES}
    return tryparse_internal(T, str, true)
end # function Base.parse

# This is an internal function and should not be exported.
regexof(::Type{<:CellParametersCard})::Regex = CELL_PARAMETERS_BLOCK

tryparsefirst(::Type{T}, str::AbstractString) where {T} = tryparse(T, str)
parsefirst(::Type{T}, str::AbstractString) where {T} = parse(T, str)

function tryparseall(::Type{T}, str::AbstractString) where {T}
    return map(eachmatch(regexof(T), str)) do x
        try
            parse(T, x.match)
        catch
            nothing
        end
    end
end # function parseall
function parseall(::Type{T}, str::AbstractString) where {T}
    return map(eachmatch(regexof(T), str)) do x
        try
            parse(T, x.match)
        catch
            Meta.ParseError("Pass failed!")
        end
    end
end # function parseall

tryparselast(::Type{T}, str::AbstractString) where {T} = tryparseall(T, str)[end]
parselast(::Type{T}, str::AbstractString) where {T} = parseall(T, str)[end]

end
