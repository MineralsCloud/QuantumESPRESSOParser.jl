"""
# module Outputs



# Examples

```jldoctest
julia>
```
"""
module Outputs

export SubroutineError, OutputFile

struct SubroutineError
    name::String
    cerr::String
    msg::String
end

struct OutputFile{A}
    source::A
end

Base.read(f::OutputFile{String}) = read(f.source, String)

include("PWscf.jl")

end
