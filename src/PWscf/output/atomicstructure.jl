export eachcellparameterscard, eachatomicpositionscard

eachcellparameterscard(str::AbstractString) =
    EachParsed{CellParametersCard}(CELL_PARAMETERS_BLOCK, str)

eachatomicpositionscard(str::AbstractString) =
    EachParsed{AtomicPositionsCard}(ATOMIC_POSITIONS_BLOCK, str)
