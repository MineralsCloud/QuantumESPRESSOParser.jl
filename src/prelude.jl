struct SubroutineError
    name::String
    cerr::String
    msg::String
end

struct ParseFailure <: Exception
    msg::String
end
