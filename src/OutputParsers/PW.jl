"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using QuantumESPRESSOParsers.Utils
using QuantumESPRESSOParsers.OutputParsers

export mark_lines, mark_ranges, PATTERNS

function mark_lines(patterns, lines)
    marks = LineNumber[]
    patterns = copy(PATTERNS)  # Cannot use `PATTERNS` directly, it will be modified!
    for (n, line) in enumerate(lines)
        for (i, p) in enumerate(patterns)
            if occursin(p, line)
                push!(marks, LineNumber(n))
                deleteat!(patterns, i)
            end
        end
    end
    marks
end # function mark_lines

function mark_ranges(patterns, lines)
    [LineRange(r) for r in split_to_ranges(collect(getfield(x, :n) for x in mark_lines(patterns, lines)))]
end # function mark_ranges

const PATTERNS = [
    r"Program PWSCF v\.(\d+)\.?(\d+)"i,
    r"Parallel version \((.*)\), running on\s+(\d+)\s+processor"i,
    r"Parallelization info"i,
    r"bravais-lattice index"i,
    r"(\d+)\s*Sym\. Ops\., with inversion, found"i,
    r"number of k points=\s*(\d+)\s*(.*)width \(Ry\)=\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)"i,
    r"Dense  grid:\s*(\d+)\s*G-vectors     FFT dimensions: \((.*),(.*),(.*)\)"i,
    r"starting charge(.*), renormalised to(.*)"i,
    r"total cpu time spent up to now is\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*secs"i,
    r"End of self-consistent calculation"i,
    r"the Fermi energy is\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*ev"i,
    r"!    total energy              =\s*([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?)\s*Ry"i,
    r"The total energy is the sum of the following terms:"i,
    r"convergence has been achieved in\s*(\d+)\s*iterations"i,
    r"Forces acting on atoms \(cartesian axes, Ry\/au\):"i,
    r"Computing stress \(Cartesian axis\) and pressure"i,
    r"Writing output data file\s*(.*)"i,
    r"This run was terminated on:\s*(.*)\s+(\w+)"i,
    r"JOB DONE\."i
]

end
