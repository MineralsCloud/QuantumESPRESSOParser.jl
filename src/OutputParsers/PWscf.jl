"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

# using Dates: DateTime, DateFormat
using DataFrames: DataFrame, GroupedDataFrame, groupby
using Fortran90Namelists.FortranToJulia
using QuantumESPRESSOBase.Cards.PWscf
using Compat: isnothing

using QuantumESPRESSOParsers: SubroutineError

export parse_summary,
       parse_fft_base_info,
       parse_ibz,
       parse_stress,
       parse_total_energy,
       parse_version,
       parse_parallel_info,
       parse_fft_dimensions,
       parse_cell_parameters,
       parse_atomic_positions,
       parse_scf_calculation,
       parse_clock,
       whatinput,
       isrelaxed,
       isjobdone,
       haserror

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}

include("regexes.jl")

const PATTERNS = [
    r"([0-9]+)\s*Sym\. Ops\., with inversion, found"i,
    r"starting charge(.*), renormalised to(.*)"i,
    r"the Fermi energy is\s*([-+]?[0-9]*\.?[0-9]+((:?[ed])[-+]?[0-9]+)?)\s*ev"i,
    r"The total energy is the sum of the following terms:"i,
    r"convergence has been achieved in\s*([0-9]+)\s*iterations"i,
    r"Forces acting on atoms \(cartesian axes, Ry\/au\):"i,
]

function parse_summary(str::AbstractString)
    dict = Dict{String,Any}()
    m = match(SUMMARY_BLOCK, str)
    if isnothing(m)
        @info("The head message is not found!") && return
    else
        content = first(m.captures)
    end

    for (f, r) in zip(
        [x -> parse(Int, x), x -> parse(Float64, x), string],
        [
         [
          BRAVAIS_LATTICE_INDEX
          NUMBER_OF_ATOMS_PER_CELL
          NUMBER_OF_ATOMIC_TYPES
          NUMBER_OF_KOHN_SHAM_STATES
          NUMBER_OF_ITERATIONS_USED
          NSTEP
         ],
         [
          LATTICE_PARAMETER
          UNIT_CELL_VOLUME
          NUMBER_OF_ELECTRONS  # TODO: This one is special.
          KINETIC_ENERGY_CUTOFF
          CHARGE_DENSITY_CUTOFF
          CUTOFF_FOR_FOCK_OPERATOR
          CONVERGENCE_THRESHOLD
          MIXING_BETA
         ],
         [EXCHANGE_CORRELATION],
        ],
    )
        for regex in r
            m = match(regex, content)
            if !isnothing(m)
                push!(dict, m.captures[1] => f(m.captures[2]))
            end
        end
    end
    return dict
end # function parse_summary

"""
    parse_fft_base_info(str::AbstractString)

Parse the FFT base information from `pw.x`'s output and return a `GroupedDataFrame`.

If there are more than one processors, the title is "Parallelization info" and three
rows, i.e., "Min", "Max", and "Sum" are printed. If not, the title is
"G-vector sticks info" and only the "Sum" row is printed. If no information is found,
return `nothing`. The `DataFrame` is grouped by "sticks" and "gvecs".
"""
function parse_fft_base_info(str::AbstractString)::Maybe{GroupedDataFrame}
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
    else
        body = m[:body]
    end
    for line in split(body, '\n')
        # "Min",4X,2I8,I7,12X,2I9,I8
        sp = split(strip(line), r"\s+")
        numbers = map(x -> parse(Int, x), sp[2:7])
        push!(df, ["sticks" sp[1] numbers[1:3]...])
        push!(df, ["gvecs" sp[1] numbers[4:6]...])
    end
    return groupby(df, :kind)
end # function parse_fft_base_info

# Return `nothing`, `(nothing, nothing)`, `(cartesian_coordinates, nothing)`, `(nothing, crystal_coordinates)`, `(cartesian_coordinates, crystal_coordinates)`
function parse_ibz(str::AbstractString)::Maybe{Tuple}
    m = match(K_POINTS_BLOCK, str)
    if isnothing(m)
        @info("The k-points info is not found!") && return
    end
    nk = parse(Int, m[:nk])

    if !isnothing(m[:cart])
        cartesian_coordinates = zeros(nk, 4)
        for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, m[:cart]))
            cartesian_coordinates[i, :] = map(x -> parse(Float64, x), m.captures[2:5])
        end
        @assert(size(cartesian_coordinates)[1] == nk)
    else
        # If there is no Cartesian coordinates, there must be no crystal coordinates.
        @info("Cartesian coordinates is `nothing`!")
        cartesian_coordinates = nothing
    end

    if !isnothing(m[:cryst])
        crystal_coordinates = zeros(nk, 4)
        for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, m[:cryst]))
            crystal_coordinates[i, :] = map(x -> parse(Float64, x), m.captures[2:5])
        end
        @assert(size(crystal_coordinates)[1] == nk)
    else
        # Only Cartesian coordinates present.
        crystal_coordinates = nothing
    end
    return cartesian_coordinates, crystal_coordinates
end # function parse_ibz

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses, kbar_stresses = Matrix{Float64}[], Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK, str)
        pressure, content = m.captures[1], m.captures[3]
        push!(pressures, parse(Float64, pressure))

        stress_atomic, stress_kbar = ntuple(_ -> Matrix{Float64}(undef, 3, 3), 2)
        for (i, line) in enumerate(split(content, '\n'))
            tmp = map(x -> parse(Float64, x), split(strip(line), " ", keepempty = false))
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
        step = Int[],
        i = Int[],
        ecut = Float64[],
        β = Float64[],
        t = Float64[],
        ɛ = Maybe{Float64}[],
        hf = Maybe{Float64}[],
        δ = Maybe{Float64}[],
    )
    for (j, m) in enumerate(eachmatch(SELF_CONSISTENT_CALCULATION_BLOCK, str))
        for (k, n) in enumerate(eachmatch(ITERATION_BLOCK, m.captures |> first))
            head = match(ITERATION_NUMBER, n.captures[1])
            isnothing(head) && continue
            i = parse(Int, head.captures[1])
            @assert(i == k)
            ecut, β = map(x -> parse(Float64, x), head.captures[2:3])
            t = parse(Float64, match(TOTAL_CPU_TIME, n.captures[1])[1])

            if !isnothing(n.captures[2])
                body = n.captures[2]
                ɛ, hf, δ = map(x -> parse(Float64, x), match(UNCONVERGED_ELECTRONS_ENERGY, body).captures)
            else
                ɛ, hf, δ = nothing, nothing, nothing
            end
            push!(df, [j k ecut β t ɛ hf δ])
        end
    end
    return groupby(df, :step)
end # function parse_scf_calculation

# See https://github.com/QEF/q-e/blob/4132a64/PW/src/print_ks_energies.f90#L10.
function parse_ks_energy(str::AbstractString)

end # function parse_ks_energy

function parse_total_energy(str::AbstractString)::Vector{Float64}
    result = Float64[]
    for m in eachmatch(
        r"!\s+total energy\s+=\s*([-+]?[0-9]*\.?[0-9]+((:?[ed])[-+]?[0-9]+)?)\s*Ry"i,
        str,
    )
        push!(result, parse(Float64, FortranData(m.captures[1])))
    end
    return result
end # function parse_total_energy

function parse_version(str::AbstractString)::Maybe{String}
    m = match(PWSCF_VERSION, str)
    !isnothing(m) ? m[:version] : return
end # function parse_version

function parse_parallel_info(str::AbstractString)::Maybe{Tuple{String,Int}}
    m = match(PARALLEL_INFO, str)
    isnothing(m) && return
    return m[:kind], isnothing(m[:num]) ? 1 : parse(Int, m[:num])
end # function parse_parallel_info

function parse_fft_dimensions(str::AbstractString)
    m = match(FFT_DIMENSIONS, str)
    !isnothing(m) ? map(x -> parse(Int, x), m.captures) : return
end # function parse_fft_dimensions

function parse_clock(str::AbstractString)::Maybe{GroupedDataFrame}
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
    return groupby(info, :subroutine)
end # function parse_clock

function whatinput(str::AbstractString)::Maybe{String}
    m = match(READING_INPUT_FROM, str)
    if isnothing(m)
        @info("The input file name is not found!") && return
    else
        return m[1]
    end
end # function whatinput

function isrelaxed(str::AbstractString)::Bool
    isnothing(match(FINAL_COORDINATES_BLOCK, str)) ? false : true
end # function isrelaxed

isjobdone(str::AbstractString) = !isnothing(match(JOB_DONE, str))

function haserror(str::AbstractString)
    isnothing(match(ERROR_IDENTIFIER, str)) ? false : true
end # function haserror

function Base.parse(::Type{T}, str::AbstractString) where {T<:SubroutineError}
    errors = T[]
    if haserror(str)
        for e in eachmatch(ERROR_BLOCK, str)
            body = strip(e[:body])
            s, msg = split(body, '\n')
            m = match(ERROR_IN_ROUTINE, s)
            push!(errors, T(m[1], m[2], strip(msg)))
        end
        return errors
    end
    return
end # function Base.parse

end
