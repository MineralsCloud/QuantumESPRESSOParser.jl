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
using MLStyle: @match
using Parameters: @with_kw
using QuantumESPRESSOBase.Cards.PWscf

using QuantumESPRESSOParsers: SubroutineError

export DiagonalizationStyle,
       Preamble,
       DavidsonDiagonalization,
       CGDiagonalization,
       PPCGDiagonalization,
       parse_fft_base_info,
       parse_ibz,
       parse_stress,
       parse_converged_energy,
       parse_version,
       parse_parallel_info,
       parse_fft_dimensions,
       parse_cell_parameters,
       parse_atomic_positions,
       parse_scf_calculation,
       parse_bands,
       parse_clock,
       whatinput,
       isrelaxed,
       isjobdone

include("regexes.jl")

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}  # Should not be exported
# Referenced from https://discourse.julialang.org/t/how-to-get-the-non-nothing-type-from-union-t-nothing/30523
nonnothingtype(::Type{T}) where {T} = Core.Compiler.typesubtract(T, Nothing)  # Should not be exported

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
        @info("The FFT base info is not found!") && return
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

# Return `nothing`, `(nothing, nothing)`, `(cartesian_coordinates, nothing)`, `(nothing, crystal_coordinates)`, `(cartesian_coordinates, crystal_coordinates)`
function parse_ibz(str::AbstractString)::Maybe{NamedTuple}
    m = match(K_POINTS_BLOCK, str)
    if isnothing(m)
        @info("The k-points info is not found!") && return
    end
    nk = parse(Int, m[:nk])
    result = []
    kinds = (:cart, :cryst)
    for k in kinds
        if !isnothing(m[k])
            x = Matrix{Float64}(undef, nk, 4)
            for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, m[k]))
                x[i, :] = map(x -> parse(Float64, x), m.captures[1:end])
            end
        else
            x = nothing
        end
        push!(result, x)
    end
    return NamedTuple{kinds}(result)
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

function parse_cell_parameters(str::AbstractString)::Vector{<:Matrix}
    cell_parameters = Matrix{Float64}[]
    for m in eachmatch(CELL_PARAMETERS_BLOCK, str)
        alat = parse(Float64, m.captures[1])
        content = m.captures[3]

        data = Matrix{Float64}(undef, 3, 3)
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM, content))
            captured = matched.captures
            data[i, :] = map(
                x -> parse(Float64, FortranData(x)),
                [captured[1], captured[4], captured[7]],
            )
        end
        push!(cell_parameters, alat * data)
    end
    return cell_parameters
end # function parse_cell_parameters

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

function parse_scf_calculation(str::AbstractString)
    df = DataFrame(
        n = Int[],  # Step number
        i = Int[],  # Iteration number
        ecut = Float64[],  # Cutoff energy
        diag = Maybe{DiagonalizationStyle}[],  # Diagonalization style
        ethr = Maybe{Float64}[],  # Energy threshold
        avg = Maybe{Float64}[],  # Average # of iterations
        β = Float64[],  # Mixing beta
        t = Float64[],  # Time
        ɛ = Maybe{Float64}[],  # Total energy
        hf = Maybe{Float64}[],  # Harris-Foulkes estimate
        δ = Maybe{Float64}[],  # Estimated scf accuracy
    )
    # (step counter, relax step)
    for (i, scf) in enumerate(eachmatch(SELF_CONSISTENT_CALCULATION_BLOCK, str))
        # (iteration counter, scf iteration)
        for (j, iter) in enumerate(eachmatch(ITERATION_BLOCK, scf[1]))
            body = iter[1]
            head = match(ITERATION_HEAD, body)
            n = parse(Int, head[1])  # Iteration number
            @assert(n == j, "Something went wrong when parsing iteration number!")
            ecut, β = map(x -> parse(Float64, x), head.captures[2:3])

            solver, ethr, avg = parse_diagonalization(body)

            t = parse(Float64, match(TOTAL_CPU_TIME, body)[1])

            ɛ, hf, δ = parse_unconverged_electrons_energy(body)

            push!(df, [i j ecut solver ethr avg β t ɛ hf δ])
        end
    end
    return groupby(df, :n)
end # function parse_scf_calculation

function parse_diagonalization(str::AbstractString)
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
end # function parse_diagonalization

function parse_unconverged_electrons_energy(str::AbstractString)
    ɛ, hf, δ = nothing, nothing, nothing  # Initialization
    m = match(UNCONVERGED_ELECTRONS_ENERGY, str)
    if !isnothing(m)
        ɛ, hf, δ = map(x -> parse(Float64, x), m.captures)
    end  # Keep them `nothing` if `m` is `nothing`
    return ɛ, hf, δ
end # function parse_UNCONVERGED_ELECTRONS_ENERGY

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
            push!(kpts, map(x -> parse(Float64, x[1]), eachmatch(Regex(GENERAL_REAL), m[:k])))
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
    ɛ, hf, δ = nothing, nothing, nothing  # Initialization
    m = match(CONVERGED_ELECTRONS_ENERGY, str)
    if !isnothing(m)
        ɛ, hf, δ = map(x -> parse(Float64, x), m.captures[1:3])
    else
        return
    end  # Keep them `nothing` if `m` is `nothing`
    regex = Regex(FIXED_POINT_REAL)
    ae = if !isnothing(m[:ae])  # 1 energy
        parse(Float64, match(regex, m[:ae])[1])
    else
        nothing
    end
    decomp = if !isnothing(m[:decomp])  # 4 energies
        map(x -> parse(Float64, x[1]), eachmatch(regex, m[:decomp]))
    else
        nothing
    end
    onecenter = if !isnothing(m[:one])  # 7 energies
        map(x -> parse(Float64, x[1]), eachmatch(regex, m[:one]))
    else
        nothing
    end
    smearing = if !isnothing(m[:smearing])  # 1 energy
        parse(Float64, match(regex, m[:smearing])[1])
    else
        nothing
    end
    return ɛ, hf, δ, ae, decomp, onecenter, smearing
end # function parse_converged_energy

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
    return parsed[1], NamedTuple{(:nr1,:nr2,:nr3)}(parsed[2:end])
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
            S = nonnothingtype(fieldtype(Preamble, field))
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
# This is an internal function and should not be exported.
function tryparse_internal(::Type{T}, str::AbstractString, raise::Bool) where {T<:SubroutineError}
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

const _INTERNAL_TYPES = Union{Preamble,SubroutineError}

function Base.tryparse(::Type{T}, str::AbstractString) where {T<:_INTERNAL_TYPES}
    return tryparse_internal(T, str, false)
end # function Base.tryparse
function Base.parse(::Type{T}, str::AbstractString) where {T<:_INTERNAL_TYPES}
    return tryparse_internal(T, str, true)
end # function Base.parse

end
