"""
# module Outputs



# Examples

```jldoctest
julia>
```
"""
module Outputs

export SubroutineError

struct SubroutineError
    name::String
    cerr::String
    msg::String
end

include("PWscf.jl")

end
