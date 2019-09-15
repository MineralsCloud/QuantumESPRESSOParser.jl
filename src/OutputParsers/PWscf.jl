"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Fortran90Namelists.FortranToJulia

using Compat: isnothing

export read_head,
       read_stress,
       read_total_energy,
       read_qe_version,
       read_processors_num,
       read_fft_dimensions,
       read_cell_parameters,
       isjobdone

const HEAD_BLOCK_REGEX = r"""
(bravais-lattice[\s\w\d\.\(\)\-\/_=^\[\]]+?)  # Match block start with "bravais-lattice", `+?` means un-greedy matching
(?=^\s*celldm)                                # Do not match any of the "celldm" pattern, must be un-greedy
"""imx
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

function read_head(str::AbstractString)
    m = match(HEAD_BLOCK_REGEX, str)
    isnothing(m) && return
    content = m.captures[1]
    dict = Dict{String,Any}()
    m = match(BRAVAIS_LATTICE_INDEX_REGEX, content)
    isnothing(m) || (dict["bravais-lattice index"] = parse(Int, m.captures[1]))
    m = match(LATTICE_PARAMETER_REGEX, str)
    isnothing(m) || (dict["lattice parameter"] = parse(Float64, FortranData(m.captures[1])))
    m = match(UNIT_CELL_VOLUME_REGEX, str)
    isnothing(m) || (dict["unit-cell volume"] = parse(Float64, FortranData(m.captures[1])))
    m = match(NUMBER_OF_ATMOS_CELL_REGEX, str)
    isnothing(m) || (dict["number of atoms/cell"] = parse(Int, m.captures[1]))
    m = match(NUMBER_OF_ATOMIC_TYPES_REGEX, str)
    isnothing(m) || (dict["number of atomic types"] = parse(Int, m.captures[1]))
    m = match(NUMBER_OF_ELECTRONS_REGEX, str)
    isnothing(m) || (dict["number of electrons"] = parse(Float64, m.captures[1]))
    m = match(NUMBER_OF_KOHN_SHAM_STATES_REGEX, str)
    isnothing(m) || (dict["number of Kohn-Sham states"] = parse(Int, m.captures[1]))
    m = match(KINETIC_ENERGY_CUTOFF_REGEX, str)
    isnothing(m) || (dict["kinetic energy cutoff"] = parse(Float64, FortranData(m.captures[1])))
    m = match(CHARGE_DENSITY_CUTOFF_REGEX, str)
    isnothing(m) || (dict["charge density cutoff"] = parse(Float64, FortranData(m.captures[1])))
    m = match(CONVERGENCE_THRESHOLD_REGEX, str)
    isnothing(m) || (dict["convergence threshold"] = parse(
        Float64,
        FortranData(m.captures[1])
    ))
    m = match(MIXING_BETA_REGEX, str)
    isnothing(m) || (dict["mixing beta"] = parse(Float64, FortranData(m.captures[1])))
    m = match(NUMBER_OF_ITERATIONS_USED_REGEX, str)
    isnothing(m) || (dict["number of iterations used"] = parse(Int, m.captures[1]))
    m = match(EXCHANGE_CORRELATION_REGEX, str)
    isnothing(m) || (dict["Exchange-correlation"] = m.captures[1] |> string)
    m = match(NSTEP_REGEX, str)
    isnothing(m) || (dict["nstep"] = parse(Int, m.captures[1]))
    return dict
end # function read_head

function read_stress(str::AbstractString)
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

function read_cell_parameters(str::AbstractString)
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

function read_total_energy(str::AbstractString)
    result = Float64[]
    for m in eachmatch(
        r"!\s+total energy\s+=\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*Ry"i,
        str
    )
        push!(result, parse(Float64, FortranData(m.captures[1])))
    end
    return result
end # function read_total_energy

function read_qe_version(line::AbstractString)
    m = match(r"Program PWSCF v\.(\d\.\d+\.?\d?)"i, line)
    isnothing(m) && error("Match error!")
    return "$(parse(Float64, FortranData(m.captures[1])))"
end # function read_qe_version

function read_processors_num(line::AbstractString)
    m = match(
        r"(?:Parallel version \((.*)\), running on\s+(\d+)\s+processor|Serial version)"i,
        line
    )
    isnothing(m) && error("Match error!")
    isnothing(m.captures) && return "Serial version"
    return m.captures[1], parse(Int, m.captures[2])
end # function read_processors_num

function read_fft_dimensions(line::AbstractString)
    m = match(
        r"Dense  grid:\s*(\d+)\s*G-vectors     FFT dimensions: \((.*),(.*),(.*)\)"i,
        line
    )
    isnothing(m) && error("Match error!")
    return map(x -> parse(Int, FortranData(x)), m.captures)
end # function read_fft_dimensions

isjobdone(str::AbstractString) = !isnothing(match(JOB_DONE_REGEX, str))

end
