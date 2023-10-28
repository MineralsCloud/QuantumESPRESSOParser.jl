using Dates: Hour, Minute, Millisecond
using VersionParsing: vparse

export QuantumESPRESSOBaseVersion,
    Parallelization,
    FFTDimensions,
    InputFile,
    SubroutineError,
    Preamble,
    FFTGrid,
    IrreducibleBrillouinZone

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

struct SubroutineError <: PWOutputItem
    name::String
    cerr::String
    msg::String
end

function Base.parse(::Type{SubroutineError}, str::AbstractString)
    obj = tryparse(SubroutineError, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
end
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
end

Base.@kwdef struct Preamble <: PWOutputItem
    ibrav::Int8
    alat::Float64
    omega::Float64
    nat::Int64
    ntyp::Int64
    nelec::Float64
    nelup::Maybe{Float64} = nothing
    neldw::Maybe{Float64} = nothing
    nbnd::Int64
    ecutwfc::Float64
    ecutrho::Float64
    ecutfock::Maybe{Float64} = nothing
    conv_thr::Maybe{Float64} = nothing
    mixing_beta::Maybe{Float64} = nothing
    mixing_ndim::Maybe{Int64} = nothing
    mixing_mode::Maybe{String} = nothing
    xc::String
    nstep::Maybe{Int64} = nothing
end

function Base.parse(::Type{Preamble}, str::AbstractString)
    obj = tryparse(Preamble, str)
    isnothing(obj) ? throw(ParseError("no matched string found!")) : return obj
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
end

struct FFTGrid <: PWOutputItem
    type::String
    dense::Int64
    smooth::Int64
    PW::Int64
end

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

struct IrreducibleBrillouinZone <: PWOutputItem
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
