"""
# module Utils



# Examples

```jldoctest
julia>
```
"""
module Utils

export split_to_ranges

function split_to_ranges(array::AbstractVector{T}) where {T <: Integer}
    backup = [array[2:end] .- one(T); typemax(T)]
    [range(s; stop = e) for (s, e) in zip(array, backup)]
end # function split_to_ranges

end
