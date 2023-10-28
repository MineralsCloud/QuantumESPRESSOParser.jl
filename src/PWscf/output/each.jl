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
