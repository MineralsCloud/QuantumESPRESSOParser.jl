export eachstep, eachiteration

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

struct UnconvergedEnergy <: PWOutputItem
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

struct ConvergedEnergy <: PWOutputItem
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
