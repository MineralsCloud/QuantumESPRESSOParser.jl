export isoptimized, isjobdone

const FINAL_COORDINATES_BLOCK = r"""
Begin final coordinates
(\X+?)
End final coordinates
"""

isoptimized(str::AbstractString) =
    match(FINAL_COORDINATES_BLOCK, str) === nothing ? false : true

isjobdone(str::AbstractString) = match(JOB_DONE, str) !== nothing
