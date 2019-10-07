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
       parse_stress,
       parse_total_energy,
       parse_qe_version,
       parse_processors_num,
       parse_fft_dimensions,
       parse_cell_parameters,
       parse_atomic_positions,
       isjobdone

const HEAD_BLOCK_REGEX = r"(bravais-lattice index\X+?)\s*celldm"i  # Match between "bravais-lattice index" and any of the "celldm" pattern, `+?` means un-greedy matching (required)
const BRAVAIS_LATTICE_INDEX_REGEX = r"bravais-lattice index\s+=\s*(-?\d+)"i
const LATTICE_PARAMETER_REGEX = r"lattice\s+parameter\s+\(alat\)\s+=\s*([\-|\+]? (?: \d*[\.]\d+ | \d+[\.]?\d*)    ([E|e|d|D][+|-]?\d+)?)\s*\w+"ix
const UNIT_CELL_VOLUME_REGEX = r"unit-cell\s+volume\s+=\s*([\-|\+]? (?: \d*[\.]\d+ | \d+[\.]?\d*)    ([E|e|d|D][+|-]?\d+)?)\s*\("ix
const NUMBER_OF_ATMOS_CELL_REGEX = r"number\s+of\s+atoms\/cell\s+=\s*(-?\d+)"i
const NUMBER_OF_ATOMIC_TYPES_REGEX = r"number\s+of\s+atomic\s+types\s*=\s*(-?\d+)"i
const NUMBER_OF_ELECTRONS_REGEX = r"number\s+of\s+electrons\s*=\s*(-?\d+[\.]\d+)"i
const NUMBER_OF_KOHN_SHAM_STATES_REGEX = r"number\s+of\s+Kohn-Sham\s+states\s*=\s*(-?\d+)"i
const KINETIC_ENERGY_CUTOFF_REGEX = r"kinetic-energy\s+cutoff\s*=\s*(-?\d+[\.]\d+)\s*Ry"i
const CHARGE_DENSITY_CUTOFF_REGEX = r"charge\s+density\s+cutoff\s*=\s*(-?\d+[\.]\d+)\s*Ry"i
const CONVERGENCE_THRESHOLD_REGEX = r"convergence\s+threshold\s+=\s*([\-|\+]? (?: \d*[\.]\d+ | \d+[\.]?\d*)    [E|e|d|D][+|-]?\d+)"ix
const MIXING_BETA_REGEX = r"mixing\s+beta\s+=\s*([+-]?(?:\d*\.\d+|\d+\.?\d*)([ED][+-]?\d+)?)"i
const NUMBER_OF_ITERATIONS_USED_REGEX = r"number of iterations used\s*=\s*(\d+)\s*plain\s*(\d*)mixing"i
const EXCHANGE_CORRELATION_REGEX = r"Exchange-correlation\s*=\s*(.*)"i
const NSTEP_REGEX = r"nstep\s+=\s*(\d+)"i
const CELL_PARAMETERS_BLOCK_REGEX = r"""
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
const CELL_PARAMETERS_ITEM_REGEX = r"""
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
const ATOMIC_POSITIONS_BLOCK_REGEX = r"""
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
const ATOMIC_POSITIONS_ITEM_REGEX = r"""
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
const STRESS_BLOCK_REGEX = r"""
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
const JOB_DONE_REGEX = r"JOB DONE\."i

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
    m = match(HEAD_BLOCK_REGEX, str)
    isnothing(m) && return
    content = m.captures[1]
    dict = Dict{String,Any}()
    m = match(BRAVAIS_LATTICE_INDEX_REGEX, content)
    isnothing(m) || (dict["bravais-lattice index"] = parse(Int, m.captures[1]))
    m = match(LATTICE_PARAMETER_REGEX, content)
    isnothing(m) || (dict["lattice parameter"] = parse(Float64, FortranData(m.captures[1])))
    m = match(UNIT_CELL_VOLUME_REGEX, content)
    isnothing(m) || (dict["unit-cell volume"] = parse(Float64, FortranData(m.captures[1])))
    m = match(NUMBER_OF_ATMOS_CELL_REGEX, content)
    isnothing(m) || (dict["number of atoms/cell"] = parse(Int, m.captures[1]))
    m = match(NUMBER_OF_ATOMIC_TYPES_REGEX, content)
    isnothing(m) || (dict["number of atomic types"] = parse(Int, m.captures[1]))
    m = match(NUMBER_OF_ELECTRONS_REGEX, content)
    isnothing(m) || (dict["number of electrons"] = parse(Float64, m.captures[1]))
    m = match(NUMBER_OF_KOHN_SHAM_STATES_REGEX, content)
    isnothing(m) || (dict["number of Kohn-Sham states"] = parse(Int, m.captures[1]))
    m = match(KINETIC_ENERGY_CUTOFF_REGEX, content)
    isnothing(m) || (dict["kinetic energy cutoff"] = parse(Float64, FortranData(m.captures[1])))
    m = match(CHARGE_DENSITY_CUTOFF_REGEX, content)
    isnothing(m) || (dict["charge density cutoff"] = parse(Float64, FortranData(m.captures[1])))
    m = match(CONVERGENCE_THRESHOLD_REGEX, content)
    isnothing(m) || (dict["convergence threshold"] = parse(
        Float64,
        FortranData(m.captures[1])
    ))
    m = match(MIXING_BETA_REGEX, content)
    isnothing(m) || (dict["mixing beta"] = parse(Float64, FortranData(m.captures[1])))
    m = match(NUMBER_OF_ITERATIONS_USED_REGEX, content)
    isnothing(m) || (dict["number of iterations used"] = parse(Int, m.captures[1]))
    m = match(EXCHANGE_CORRELATION_REGEX, content)
    isnothing(m) || (dict["Exchange-correlation"] = m.captures[1] |> string)
    m = match(NSTEP_REGEX, content)
    isnothing(m) || (dict["nstep"] = parse(Int, m.captures[1]))
    return dict
end # function read_head

function parse_stress(str::AbstractString)
    pressures = Float64[]
    atomic_stresses = Matrix{Float64}[]
    kbar_stresses = Matrix{Float64}[]
    for m in eachmatch(STRESS_BLOCK_REGEX, str)
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
    for m in eachmatch(CELL_PARAMETERS_BLOCK_REGEX, str)
        alat = parse(Float64, m.captures[1])
        content = m.captures[3]

        data = Matrix{Float64}(undef, 3, 3)
        for (i, matched) in enumerate(eachmatch(CELL_PARAMETERS_ITEM_REGEX, content))
            captured = matched.captures
            data[i, :] = map(
                x -> parse(Float64, FortranData(x)),
                [captured[1], captured[4], captured[7]]
            )
        end
        push!(cell_parameters, alat * data)
    end
    return cell_parameters
end # function read_cell_parameters

function parse_atomic_positions(str::AbstractString)
    atomic_positions = AtomicPositionsCard[]
    for m in eachmatch(ATOMIC_POSITIONS_BLOCK_REGEX, str)
        unit = string(m.captures[1])
        content = m.captures[2]
        data = AtomicPosition[]

        for matched in eachmatch(ATOMIC_POSITIONS_ITEM_REGEX, content)
            captured = matched.captures
            if_pos = map(x -> isempty(x) ? 1 : parse(Int, FortranData(x)), captured[11:13])
            atom, pos = string(captured[1]),
                map(
                    x -> parse(Float64, FortranData(x)),
                    [captured[3], captured[6], captured[9]]
                )
            push!(data, AtomicPosition(atom, pos, if_pos))
        end
        push!(atomic_positions, AtomicPositionsCard(unit, data))
    end
    return atomic_positions
end

function parse_total_energy(str::AbstractString)
    result = Float64[]
    for m in eachmatch(
        r"!\s+total energy\s+=\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*Ry"i,
        str
    )
        push!(result, parse(Float64, FortranData(m.captures[1])))
    end
    return result
end # function read_total_energy

function parse_qe_version(line::AbstractString)
    m = match(r"Program PWSCF v\.(\d\.\d+\.?\d?)"i, line)
    isnothing(m) && error("Match error!")
    return "$(parse(Float64, FortranData(m.captures[1])))"
end # function read_qe_version

function parse_processors_num(line::AbstractString)
    m = match(
        r"(?:Parallel version \((.*)\), running on\s+(\d+)\s+processor|Serial version)"i,
        line
    )
    isnothing(m) && error("Match error!")
    isnothing(m.captures) && return "Serial version"
    return m.captures[1], parse(Int, m.captures[2])
end # function read_processors_num

function parse_fft_dimensions(line::AbstractString)
    m = match(
        r"Dense  grid:\s*(\d+)\s*G-vectors     FFT dimensions: \((.*),(.*),(.*)\)"i,
        line
    )
    isnothing(m) && error("Match error!")
    return map(x -> parse(Int, FortranData(x)), m.captures)
end # function read_fft_dimensions

isjobdone(str::AbstractString) = !isnothing(match(JOB_DONE_REGEX, str))

end
