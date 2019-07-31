"""
# module OutputParsers



# Examples

```jldoctest
julia>
```
"""
module OutputParsers

export Linemark

struct Linemark
    p::Regex
    n::Int
end

include("PW.jl")

end
