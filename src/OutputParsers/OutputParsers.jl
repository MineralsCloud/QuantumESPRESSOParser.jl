"""
# module OutputParsers



# Examples

```jldoctest
julia>
```
"""
module OutputParsers

export LineNumber, LineRange

struct LineNumber{T<:Integer}
    n::T
end

struct LineRange{T<:AbstractRange}
    r::T
end

include("PWscf.jl")

end
