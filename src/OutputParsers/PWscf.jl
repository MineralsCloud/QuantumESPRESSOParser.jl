"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Fortran90Namelists.FortranToJulia
using QuantumESPRESSOBase.Cards.PWscf
using Compat: isnothing

export parse_head,
       parse_parallelization_info,
       parse_k_points,
       parse_stress,
       parse_total_energy,
       parse_qe_version,
       parse_processors_num,
       parse_fft_dimensions,
       parse_cell_parameters,
       parse_atomic_positions,
       parse_scf_calculation,
       isjobdone

# See https://gist.github.com/singularitti/e9e04c501ddfe40ba58917a754707b2e
const INTEGER = raw"([+-]?\d+)"
const FIXED_POINT_REAL = raw"([-+]?\d*\.\d+|\d+\.?\d*)"
const REAL_WITH_EXPONENT = raw"([-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?[0-9]+)?)"
# The following format is from https://github.com/QEF/q-e/blob/7357cdb/PW/src/summary.f90#L100-L119.
const HEAD_BLOCK = r"(bravais-lattice index\X+?)\s*celldm"i  # Match between "bravais-lattice index" and any of the "celldm" pattern, `+?` means un-greedy matching (required)
# 'bravais-lattice index     = ',I12
const BRAVAIS_LATTICE_INDEX = Regex(raw"(bravais-lattice index)\s*=\s*" * INTEGER, "i")
# 'lattice parameter (alat)  = ',F12.4,'  a.u.'
const LATTICE_PARAMETER = Regex(
    raw"(lattice parameter \(alat\))\s*=\s*" * FIXED_POINT_REAL,
    "i",
)
# 'unit-cell volume          = ',F12.4,' (a.u.)^3'
const UNIT_CELL_VOLUME = Regex(raw"(unit-cell volume)\s+=\s*" * FIXED_POINT_REAL, "i")
# 'number of atoms/cell      = ',I12
const NUMBER_OF_ATOMS_PER_CELL = Regex(raw"(number of atoms\/cell)\s+=\s*" * INTEGER, "i")
# 'number of atomic types    = ',I12
const NUMBER_OF_ATOMIC_TYPES = Regex(raw"(number of atomic types)\s*=\s*" * INTEGER, "i")
# 'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')'
const NUMBER_OF_ELECTRONS = Regex(raw"(number of electrons)\s*=\s*" * FIXED_POINT_REAL *
                                  raw"(?:\(up:\s*" * FIXED_POINT_REAL * raw", down:\s*" *
                                  FIXED_POINT_REAL * raw"\))?")
# 'number of Kohn-Sham states= ',I12
const NUMBER_OF_KOHN_SHAM_STATES = Regex(
    raw"(number of Kohn-Sham states)\s*=\s*" * INTEGER,
    "i",
)
# 'kinetic-energy cutoff     = ',F12.4,'  Ry'
const KINETIC_ENERGY_CUTOFF = Regex(
    raw"(kinetic-energy cutoff)\s*=\s*" * FIXED_POINT_REAL * raw"\s+Ry",
    "i",
)
# 'charge density cutoff     = ',F12.4,'  Ry'
const CHARGE_DENSITY_CUTOFF = Regex(
    raw"(charge density cutoff)\s*=\s*" * FIXED_POINT_REAL * raw"\s+Ry",
    "i",
)
# 'cutoff for Fock operator  = ',F12.4,'  Ry'
const CUTOFF_FOR_FOCK_OPERATOR = Regex(
    raw"(cutoff for Fock operator)\s*=\s*" * FIXED_POINT_REAL * raw"\s+Ry",
    "i",
)
# 'convergence threshold     = ',1PE12.1
const CONVERGENCE_THRESHOLD = Regex(
    raw"(convergence threshold)\s*=\s*" * REAL_WITH_EXPONENT,
    "i",
)
# 'mixing beta               = ',0PF12.4
const MIXING_BETA = Regex(raw"(mixing beta)\s*=\s*" * FIXED_POINT_REAL, "i")
# 'number of iterations used = ',I12,2X,A,' mixing'
const NUMBER_OF_ITERATIONS_USED = Regex(
    raw"(number of iterations used)\s*=\s*" * INTEGER,
    "i",
)
const EXCHANGE_CORRELATION = r"(Exchange-correlation)\s*=\s*(.*)"i
# "nstep                     = ",I12
const NSTEP = Regex(raw"(nstep)\s*=\s*" * INTEGER, "i")
const PARALLELIZATION_INFO_BLOCK = r"""Parallelization info
\s*--------------------
\s*sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
(\X+?)
\s*bravais-lattice index"""im
const K_POINTS_BLOCK = r"""
number of k points=\s*(\d+)\X+?
\s*cart\. coord\. in units 2pi\/alat\s*
(\X+?)
\s*cryst\. coord\.\s*
(\X+?)
\s*Dense  grid"""im
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L353-L354.
# '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)'
const K_POINTS_ITEM = r"k\(\s*(\d+)\s*\) = \(\s*([-+]?\d*\.\d+|\d+\.?\d*)\s*([-+]?\d*\.\d+|\d+\.?\d*)\s*([-+]?\d*\.\d+|\d+\.?\d*)\s*\), wk =\s*([-+]?\d*\.\d+|\d+\.?\d*)"i
const CELL_PARAMETERS_BLOCK = r"""
^ [ \t]*
CELL_PARAMETERS [ \t]*
\(?\w+\s*=\s*[\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?\)? \s* [\n]
(
(
\s*
(
[\-|\+]? ( \d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?\s*
){3}[\n]
){3}
)
"""imx
const CELL_PARAMETERS_ITEM = r"""
^                        # Linestart
[ \t]*                   # Optional white space
(?P<x>                   # Get x
    [\-|\+]? ( \d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                   # Get y
    [\-|\+]? (\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                   # Get z
    [\-|\+]? (\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
"""mx
const ATOMIC_POSITIONS_BLOCK = r"""
^ \s* ATOMIC_POSITIONS \s*                      # Atomic positions start with that string
[{(]? \s* (?P<units>\S+?)? \s* [)}]? \s* $\n    # The units are after the string in optional brackets
(?P<block>                                      # This is the block of positions
    (
        (
            \s*                                 # White space in front of the element spec is ok
            (
                [A-Za-z]+[A-Za-z0-9]{0,2}       # Element spec
                (
                    \s+                         # White space in front of the number
                    [-|+]?                      # Plus or minus in front of the number (optional)
                    (
                        (
                            \d*                 # optional decimal in the beginning .0001 is ok, for example
                            [\.]                # There has to be a dot followed by
                            \d+                 # at least one decimal
                        )
                        |                       # OR
                        (
                            \d+                 # at least one decimal, followed by
                            [\.]?               # an optional dot ( both 1 and 1. are fine)
                            \d*                 # And optional number of decimals (1.00001)
                        )                        # followed by optional decimals
                    )
                    ([E|e|d|D][+|-]?\d+)?       # optional exponents E+03, e-05
                ){3}                            # I expect three float values
                ((\s+[0-1]){3}\s*)?             # Followed by optional ifpos
                \s*                             # Followed by optional white space
                |
                \#.*                            # If a line is commented out, that is also ok
                |
                \!.*                            # Comments also with excl. mark in fortran
            )
            |                                   # OR
            \s*                                 # A line only containing white space
         )
        [\n]                                    # line break at the end
    )+                                          # A positions block should be one or more lines
)
"""imx
const ATOMIC_POSITIONS_ITEM = r"""
^                                       # Linestart
[ \t]*                                  # Optional white space
(?P<name>[A-Za-z]+[A-Za-z0-9]{0,2})\s+   # get the symbol, max 3 chars, starting with a char
(?P<x>                                  # Get x
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                                  # Get y
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                                  # Get z
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]*
(?P<fx>[01]?)                           # Get fx
[ \t]*
(?P<fy>[01]?)                           # Get fx
[ \t]*
(?P<fz>[01]?)                           # Get fx
"""mx
const STRESS_BLOCK = r"""
^[ \t]*
total\s+stress\s*\(Ry\/bohr\*\*3\)\s+
\(kbar\)\s+P=\s*([\-|\+]? (?: \d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?)
[\n]
(
(?:
\s*
(?:
[\-|\+]? (?: \d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?[ \t]*
){6}
){3}
)
"""imx
const SELF_CONSISTENT_CALCULATION_BLOCK = r"Self-consistent Calculation(\X+?)End of self-consistent calculation"i
const ITERATION_BLOCK = r"(iteration #\X+?secs)\s*(total energy\X+?estimated scf accuracy    <.*)?"i
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L920-L921.
# '     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2
const ITERATION_NUMBER_ITEM = Regex(raw"iteration #\s*" * INTEGER * raw"\s+ecut=\s*" * FIXED_POINT_REAL * raw" Ry\s+beta=\s*" * FIXED_POINT_REAL, "i")
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L917-L918.
# '     total cpu time spent up to now is ',F10.1,' secs'
const TOTAL_CPU_TIME = Regex(raw"total cpu time spent up to now is\s*" * FIXED_POINT_REAL * raw"\s* secs", "i")
const JOB_DONE = r"JOB DONE\."i

const PATTERNS = [
    r"Program PWSCF v\.(\d\.\d+\.?\d?)"i,
    r"Parallelization info"i,
    r"(\d+)\s*Sym\. Ops\., with inversion, found"i,
    r"number of k points=\s*(\d+)\s*(.*)width \(Ry\)=\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)"i,
    r"starting charge(.*), renormalised to(.*)"i,
    r"total cpu time spent up to now is\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*secs"i,
    r"End of self-consistent calculation"i,
    r"the Fermi energy is\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*ev"i,
    r"The total energy is the sum of the following terms:"i,
    r"convergence has been achieved in\s*(\d+)\s*iterations"i,
    r"Forces acting on atoms \(cartesian axes, Ry\/au\):"i,
    r"Writing output data file\s*(.*)"i,
    r"This run was terminated on:\s*(.*)\s+(\w+)"i,
]

function parse_head(str::AbstractString)
    dict = Dict{String,Any}()
    m = match(HEAD_BLOCK, str)
    if isnothing(m)
        @info("The head message is not found!")
        return
    else
        content = first(m.captures)
    end

    function _parse_by(f::Function, r::AbstractVector)
        for regex in r
            m = match(regex, content)
            if !isnothing(m)
                push!(dict, m.captures[1] => f(m.captures[2]))
            end
        end
    end # function _parse_by

    _parse_by(
        x -> parse(Int, x),
        [
         BRAVAIS_LATTICE_INDEX
         NUMBER_OF_ATOMS_PER_CELL
         NUMBER_OF_ATOMIC_TYPES
         NUMBER_OF_KOHN_SHAM_STATES
         NUMBER_OF_ITERATIONS_USED
         NSTEP
        ],
    )
    _parse_by(
        x -> parse(Float64, x),
        [
         LATTICE_PARAMETER
         UNIT_CELL_VOLUME
         NUMBER_OF_ELECTRONS  # TODO: This one is special.
         KINETIC_ENERGY_CUTOFF
         CHARGE_DENSITY_CUTOFF
         CUTOFF_FOR_FOCK_OPERATOR
         CONVERGENCE_THRESHOLD
         MIXING_BETA
        ],
    )
    _parse_by(string, [EXCHANGE_CORRELATION])
    return dict
end # function parse_head

function parse_parallelization_info(str::AbstractString)
    sticks, gvecs = Dict{String,NamedTuple}(), Dict{String,NamedTuple}()
    m = match(PARALLELIZATION_INFO_BLOCK, str)
    if isnothing(m)
        @info("The parallelization info is not found!")
        return
    else
        content = first(m.captures)
    end

    for line in split(content, '\n')
        # The following format is from https://github.com/QEF/q-e/blob/7357cdb/Modules/fft_base.f90#L73-L90.
        # "Min",4X,2I8,I7,12X,2I9,I8
        sp = split(strip(line), r"\s+")
        numbers = map(x -> parse(Int, x), sp[2:7])
        push!(sticks, sp[1] => (dense = numbers[1], smooth = numbers[2], PW = numbers[3]))
        push!(gvecs, sp[1] => (dense = numbers[4], smooth = numbers[5], PW = numbers[6]))
    end
    return sticks, gvecs
end # function parse_parallelization_info

function parse_k_points(str::AbstractString)
    m = match(K_POINTS_BLOCK, str)
    if isnothing(m)
        @info("The k-points info is not found!")
        return
    else
        nks, cartesian, crystal = m.captures
    end
    nks = parse(Int, nks)

    cartesian_coordinates, crystal_coordinates = zeros(nks, 4), zeros(nks, 4)
    for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, cartesian))
        cartesian_coordinates[i, :] = map(x -> parse(Float64, x), m.captures[2:5])
    end
    for (i, m) in enumerate(eachmatch(K_POINTS_ITEM, crystal))
        crystal_coordinates[i, :] = map(x -> parse(Float64, x), m.captures[2:5])
    end
    @assert(size(cartesian_coordinates)[1] == size(crystal_coordinates)[1] == nks)
    return cartesian_coordinates, crystal_coordinates
end # function parse_k_points

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses = Matrix{Float64}[]
    kbar_stresses = Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK, str)
        pressure, content = m.captures[1], m.captures[3]
        push!(pressures, parse(Float64, pressure))

        stress_atomic = Matrix{Float64}(undef, 3, 3)
        stress_kbar = Matrix{Float64}(undef, 3, 3)
        for (i, line) in enumerate(split(content, '\n'))
            tmp = map(x -> parse(Float64, x), split(strip(line), " ", keepempty = false))
            stress_atomic[i, :], stress_kbar[i, :] = tmp[1:3], tmp[4:6]
        end
        push!(atomic_stresses, stress_atomic)
        push!(kbar_stresses, stress_kbar)
    end
    return pressures, atomic_stresses, kbar_stresses
end # function parse_stress

function parse_cell_parameters(str::AbstractString)
    cell_parameters = Matrix{Float64}[]
    for m in eachmatch(CELL_PARAMETERS_BLOCK, str)
        alat = parse(Float64, m.captures[1])
        content = m.captures[3]

        data = Matrix{Float64}(undef, 3, 3)
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM, content))
            captured = matched.captures
            data[i, :] = map(
                x -> parse(Float64, FortranData(x)),
                [captured[1], captured[4], captured[7]],
            )
        end
        push!(cell_parameters, alat * data)
    end
    return cell_parameters
end # function parse_cell_parameters

function parse_atomic_positions(str::AbstractString)
    atomic_positions = AtomicPositionsCard[]
    for m in eachmatch(ATOMIC_POSITIONS_BLOCK, str)
        unit = string(m.captures[1])
        content = m.captures[2]
        data = AtomicPosition[]

        for matched in eachmatch(ATOMIC_POSITIONS_ITEM, content)
            captured = matched.captures
            if_pos = map(
                x -> isempty(x) ? 1 : parse(Int, FortranData(x)),
                captured[11:13],
            )
            atom, pos = string(captured[1]),
                map(
                    x -> parse(Float64, FortranData(x)),
                    [captured[3], captured[6], captured[9]],
                )
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        push!(atomic_positions, AtomicPositionsCard(unit, data))
    end
    return atomic_positions
end # parse_atomic_positions

function parse_scf_calculation(str::AbstractString)
    scf_calculations = []
    for m in eachmatch(SELF_CONSISTENT_CALCULATION_BLOCK, str)
        iterations = Dict{String,Any}[]
        for n in eachmatch(ITERATION_BLOCK, m.captures |> first)
            d = Dict{String,Any}()
            head = match(ITERATION_NUMBER_ITEM, n.captures[1])
            isnothing(head) && continue
            d["iteration"] = parse(Int, head.captures[1])
            d["ecut"], d["beta"] = map(x -> parse(Float64, x), head.captures[2:3])
            time = match(TOTAL_CPU_TIME, n.captures[1])
            d["time"] = parse(Float64, time.captures[1])

            if !isnothing(n.captures[2])
                body = n.captures[2]
                e = parse(Float64, match(r"total energy\s+=\s*([-+]?\d*\.\d+|\d+\.?\d*)"i, body).captures[1])
                hf = parse(Float64, match(r"Harris-Foulkes estimate\s+=\s*([-+]?\d*\.\d+|\d+\.?\d*)"i, body).captures[1])
                ac = parse(Float64, match(Regex(raw"estimated scf accuracy\s+<\s*" * REAL_WITH_EXPONENT, "i"), body).captures[1])
                d["total energy"] = e
                d["Harris-Foulkes estimate"] = hf
                d["estimated scf accuracy"] = ac
            end
            push!(iterations, d)
        end
        push!(scf_calculations, iterations)
    end
    return scf_calculations
end # function parse_scf_calculation

function parse_total_energy(str::AbstractString)
    result = Float64[]
    for m in eachmatch(
        r"!\s+total energy\s+=\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*Ry"i,
        str,
    )
        push!(result, parse(Float64, FortranData(m.captures[1])))
    end
    return result
end # function parse_total_energy

function parse_qe_version(str::AbstractString)
    m = match(r"Program PWSCF v\.(\d\.\d+\.?\d?)"i, str)
    !isnothing(m) ? string(m.captures[1]) : return
end # function parse_qe_version

function parse_processors_num(str::AbstractString)
    m = match(
        r"(?:Parallel version \((.*)\), running on\s+(\d+)\s+processor|Serial version)"i,
        str,
    )
    isnothing(m) && return
    isnothing(m.captures) && return "Serial version", 1
    return string(m.captures[1]), parse(Int, m.captures[2])
end # function parse_processors_num

function parse_fft_dimensions(str::AbstractString)
    m = match(
        r"Dense  grid:\s*(\d+)\s*G-vectors     FFT dimensions: \((.*),(.*),(.*)\)"i,
        str,
    )
    !isnothing(m) ? map(x -> parse(Int, x), m.captures) : return
end # function parse_fft_dimensions

isjobdone(str::AbstractString) = !isnothing(match(JOB_DONE, str))

end
