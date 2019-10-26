struct SubroutineError <: Exception
    name::String
    cerr::String
    msg::String
end
