export eachstep,
    eachiteration,
    eachiterationhead,
    eachiterationtime,
    eachdiagonalization,
    eachunconvergedenergy,
    eachconvergedenergy,
    eachtimeditem

abstract type Each end

struct EachStep <: Each
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

struct EachIteration <: Each
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

struct EachParsed{T} <: Each
    regex::Regex
    string::String
end

function Base.iterate(regexmatchiterator::EachParsed{T}) where {T}
    regexmatchiterator = eachmatch(regexmatchiterator.regex, regexmatchiterator.string)
    iterated = iterate(regexmatchiterator)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return tryparse(T, matched.match), state
    end
end
function Base.iterate(regexmatchiterator::EachParsed{T}, state) where {T}
    regexmatchiterator = eachmatch(regexmatchiterator.regex, regexmatchiterator.string)
    iterated = iterate(regexmatchiterator, state)
    if isnothing(iterated)
        return nothing
    else
        matched, state = iterated
        return tryparse(T, matched.match), state
    end
end

Base.eltype(::Type{EachParsed{T}}) where {T} = T

Base.IteratorSize(::Type{<:EachParsed}) = Base.SizeUnknown()

struct IterationHead <: PWOutputItem
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

eachiterationhead(str::AbstractString) = EachParsed{IterationHead}(ITERATION_HEAD, str)

struct IterationTime <: PWOutputItem
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

eachiterationtime(str::AbstractString) = EachParsed{IterationTime}(TOTAL_CPU_TIME, str)

abstract type DiagonalizationSolver end
struct Davidson <: DiagonalizationSolver end
struct ConjugateGradient <: DiagonalizationSolver end
struct ProjectedPreconditionedConjugateGradient <: DiagonalizationSolver end

struct Diagonalization <: PWOutputItem
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
    if isnothing(matched)
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

eachdiagonalization(str::AbstractString) = EachParsed{Diagonalization}(C_BANDS, str)

Base.@kwdef struct UnconvergedEnergy <: PWOutputItem
    total::Float64
    harris_foulkes_estimate::Maybe{Float64} = nothing
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
        total, harris_foulkes_estimate, estimated_scf_accuracy = map(
            _parser, matched.captures
        )
        return UnconvergedEnergy(total, harris_foulkes_estimate, estimated_scf_accuracy)
    end
end

eachunconvergedenergy(str::AbstractString) =
    EachParsed{UnconvergedEnergy}(UNCONVERGED_ELECTRONS_ENERGY, str)

_parser(x) = isnothing(x) ? x : parse(Float64, x)

Base.@kwdef struct ConvergedEnergy <: PWOutputItem
    total::Float64
    harris_foulkes_estimate::Maybe{Float64} = nothing
    estimated_scf_accuracy::Float64
    all_electron::Maybe{Float64} = nothing
    one_electron::Maybe{Float64} = nothing
    hartree::Maybe{Float64} = nothing
    xc::Maybe{Float64} = nothing
    ewald::Maybe{Float64} = nothing
    one_center_paw::Maybe{Float64} = nothing
    paw_hartree_ae::Maybe{Float64} = nothing
    paw_hartree_ps::Maybe{Float64} = nothing
    paw_xc_ae::Maybe{Float64} = nothing
    paw_xc_ps::Maybe{Float64} = nothing
    total_e_h_paw::Maybe{Float64} = nothing
    total_e_xc_paw::Maybe{Float64} = nothing
    smearing::Maybe{Float64} = nothing
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
        total, harris_foulkes_estimate, estimated_scf_accuracy = map(
            _parser,
            (matched[1], matched[2], matched[4]),  # Not `matched[3]`!
        )
        all_electron = if !isnothing(matched[:ae])
            parse(Float64, only(match(Regex(FIXED_POINT_REAL), matched[:ae])))
        end
        one_electron, hartree, xc, ewald = if !isnothing(matched[:decomp])
            map(
                Base.Fix1(parse, Float64) ∘ only,
                eachmatch(Regex(FIXED_POINT_REAL), matched[:decomp]),
            )
        end
        # one_center_paw,
        # paw_hartree_ae, paw_hartree_ps, paw_xc_ae, paw_xc_ps, total_e_h_paw,
        # total_e_xc_paw = if !isnothing(matched[:decomp])
        #     map(
        #         Base.Fix1(parse, Float64) ∘ only,
        #         eachmatch(Regex(FIXED_POINT_REAL), matched[:decomp]),
        #     )
        # end
        smearing = if !isnothing(matched[:smearing])
            parse(Float64, only(match(Regex(FIXED_POINT_REAL), matched[:smearing])))
        end
        return ConvergedEnergy(;
            total,
            harris_foulkes_estimate,
            estimated_scf_accuracy,
            all_electron,
            one_electron,
            hartree,
            xc,
            ewald,
            # one_center_paw,
            # paw_hartree_ae,
            # paw_hartree_ps,
            # paw_xc_ae,
            # paw_xc_ps,
            # total_e_h_paw,
            # total_e_xc_paw,
            smearing,
        )
    end
end

eachconvergedenergy(str::AbstractString) =
    EachParsed{ConvergedEnergy}(CONVERGED_ELECTRONS_ENERGY, str)

struct TimedItem <: PWOutputItem
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

eachtimeditem(str::AbstractString) = EachParsed{TimedItem}(TIMED_ITEM, str)

# See https://docs.julialang.org/en/v1/manual/types/#man-custom-pretty-printing
function Base.show(io::IO, iter::Each)
    if get(io, :compact, false)
        print(IOContext(io, :limit => true, :compact => true), summary(iter), "(...)")
    else
        print(io, summary(iter), "(...)")
    end
end
