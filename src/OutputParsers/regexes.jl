# See https://gist.github.com/singularitti/e9e04c501ddfe40ba58917a754707b2e
const INTEGER = raw"([-+]?[0-9]+)"
const FIXED_POINT_REAL = raw"([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)"
const GENERAL_REAL = raw"([-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[eE][-+]?[0-9]+)?)"
const EQUAL_SIGN = raw"\s*=\s*"

# This format is from https://github.com/QEF/q-e/blob/4132a64/Modules/environment.f90#L215-L224.
const PARALLEL_INFO = r"(?<kind>(?:Parallel version [^,]*|Serial version))(?:, running on\s*(?<num>[0-9]+) processors)?"i
const PWSCF_VERSION = r"Program PWSCF v\.(?<version>[0-9]\.[0-9]+\.?[0-9]?)"i
const FFT_DIMENSIONS = r"Dense  grid:\s*([0-9]+)\s*G-vectors     FFT dimensions: \((.*),(.*),(.*)\)"i
# The following format is from https://github.com/QEF/q-e/blob/7357cdb/PW/src/summary.f90#L100-L119.
const HEAD_BLOCK = r"(bravais-lattice index\X+?)\s*celldm"i  # Match between "bravais-lattice index" & the 1st of the "celldm"s, `+?` means un-greedy matching (required)
# 'bravais-lattice index     = ',I12
const BRAVAIS_LATTICE_INDEX = Regex("(bravais-lattice index)" * EQUAL_SIGN * INTEGER, "i")
# 'lattice parameter (alat)  = ',F12.4,'  a.u.'
const LATTICE_PARAMETER = Regex(
    raw"(lattice parameter \(alat\))" * EQUAL_SIGN * FIXED_POINT_REAL,
    "i",
)
# 'unit-cell volume          = ',F12.4,' (a.u.)^3'
const UNIT_CELL_VOLUME = Regex("(unit-cell volume)" * EQUAL_SIGN * FIXED_POINT_REAL, "i")
# 'number of atoms/cell      = ',I12
const NUMBER_OF_ATOMS_PER_CELL = Regex(
    raw"(number of atoms\/cell)" * EQUAL_SIGN * INTEGER,
    "i",
)
# 'number of atomic types    = ',I12
const NUMBER_OF_ATOMIC_TYPES = Regex("(number of atomic types)" * EQUAL_SIGN * INTEGER, "i")
# 'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')'
const NUMBER_OF_ELECTRONS = Regex("(number of electrons)" * EQUAL_SIGN * FIXED_POINT_REAL *
                                  raw"(?:\(up:\s*" * FIXED_POINT_REAL * raw", down:\s*" *
                                  FIXED_POINT_REAL * raw"\))?")
# 'number of Kohn-Sham states= ',I12
const NUMBER_OF_KOHN_SHAM_STATES = Regex(
    "(number of Kohn-Sham states)" * EQUAL_SIGN * INTEGER,
    "i",
)
# 'kinetic-energy cutoff     = ',F12.4,'  Ry'
const KINETIC_ENERGY_CUTOFF = Regex(
    "(kinetic-energy cutoff)" * EQUAL_SIGN * FIXED_POINT_REAL * "\\s+Ry",
    "i",
)
# 'charge density cutoff     = ',F12.4,'  Ry'
const CHARGE_DENSITY_CUTOFF = Regex(
    "(charge density cutoff)" * EQUAL_SIGN * FIXED_POINT_REAL * "\\s+Ry",
    "i",
)
# 'cutoff for Fock operator  = ',F12.4,'  Ry'
const CUTOFF_FOR_FOCK_OPERATOR = Regex(
    "(cutoff for Fock operator)" * EQUAL_SIGN * FIXED_POINT_REAL * "\\s+Ry",
    "i",
)
# 'convergence threshold     = ',1PE12.1
const CONVERGENCE_THRESHOLD = Regex(
    "(convergence threshold)" * EQUAL_SIGN * GENERAL_REAL,
    "i",
)
# 'mixing beta               = ',0PF12.4
const MIXING_BETA = Regex("(mixing beta)" * EQUAL_SIGN * FIXED_POINT_REAL, "i")
# 'number of iterations used = ',I12,2X,A,' mixing'
const NUMBER_OF_ITERATIONS_USED = Regex(
    "(number of iterations used)" * EQUAL_SIGN * INTEGER,
    "i",
)
const EXCHANGE_CORRELATION = r"(Exchange-correlation)\s*=\s*(.*)"i
# "nstep                     = ",I12
const NSTEP = Regex("(nstep)" * EQUAL_SIGN * INTEGER, "i")
const PARALLELIZATION_INFO_BLOCK = Regex("""Parallelization info
\\s*--------------------
\\s*sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
(\\s*Min.*
\\s*Max.*
\\s*Sum.*)
""", "im")
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L341-L381.
const K_POINTS_BLOCK = r"""
number of k points=\s*(?<nk>[0-9]+)\h*(?<metainfo>.*)
\s*(?:cart\. coord\. in units 2pi\/alat\s*(?<cart>\X+?)^\s*$|Number of k-points >= 100: set verbosity='high' to print them\.)
\s*(?:cryst\. coord\.\s*(?<cryst>\X+?)
\s*Dense  grid)?"""im
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L353-L354.
# '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)'
const K_POINTS_ITEM = r"k\(\s*([0-9]+)\s*\) = \(\s*([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)\s*([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)\s*([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)\s*\), wk =\s*([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)"i
const CELL_PARAMETERS_BLOCK = r"""
^ [ \t]*
CELL_PARAMETERS [ \t]*
\(?\w+\s*=\s*[\-|\+]?([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?\)? \s* [\n]
(
(
\s*
(
[\-|\+]? ( [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?\s*
){3}[\n]
){3}
)
"""imx
const CELL_PARAMETERS_ITEM = r"""
^                        # Linestart
[ \t]*                   # Optional white space
(?P<x>                   # Get x
    [\-|\+]? ( [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
)
[ \t]+
(?P<y>                   # Get y
    [\-|\+]? ([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
)
[ \t]+
(?P<z>                   # Get z
    [\-|\+]? ([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
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
                            [0-9]*                 # optional decimal in the beginning .0001 is ok, for example
                            [\.]                # There has to be a dot followed by
                            [0-9]+                 # at least one decimal
                        )
                        |                       # OR
                        (
                            [0-9]+                 # at least one decimal, followed by
                            [\.]?               # an optional dot ( both 1 and 1. are fine)
                            [0-9]*                 # And optional number of decimals (1.00001)
                        )                        # followed by optional decimals
                    )
                    ([E|e|d|D][+|-]?[0-9]+)?       # optional exponents E+03, e-05
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
    [\-|\+]?([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
)
[ \t]+
(?P<y>                                  # Get y
    [\-|\+]?([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
)
[ \t]+
(?P<z>                                  # Get z
    [\-|\+]?([0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?
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
\(kbar\)\s+P=\s*([\-|\+]? (?: [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?)
[\n]
(
(?:
\s*
(?:
[\-|\+]? (?: [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?[ \t]*
){6}
){3}
)
"""imx
const SELF_CONSISTENT_CALCULATION_BLOCK = r"Self-consistent Calculation(\X+?)End of self-consistent calculation"i
const ITERATION_BLOCK = r"(iteration #\X+?secs)\s*(total energy\X+?estimated scf accuracy    <.*)?"i
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L920-L921.
# '     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2
const ITERATION_NUMBER_ITEM = Regex(
    "iteration #\\s*" * INTEGER * "\\s+ecut=\\s*" * FIXED_POINT_REAL * " Ry\\s+beta=\\s*" *
    FIXED_POINT_REAL,
    "i",
)
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L917-L918.
# '     total cpu time spent up to now is ',F10.1,' secs'
const TOTAL_CPU_TIME = Regex(
    "total cpu time spent up to now is\\s*" * FIXED_POINT_REAL * "\\s* secs",
    "i",
)
const TIME_BLOCK = r"(init_run\X+?This run was terminated on:.*)"i
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L29-L33.
const SUMMARY_TIME_BLOCK = r"""
(?<head>)
(?<body>
init_run\s+:.*
\s*electrons\s+:.*
\s*(?:update_pot\s+.*)?  # This does not always exist.
\s*(?:forces\s+:.*)?     # This does not always exist.
\s*(?:stress\s+:.*)?     # This does not always exist.
)
"""imx
const TIME_ITEM = Regex(
    raw"\s*([\w0-9:]+)\s+:\s*" * FIXED_POINT_REAL * "s\\sCPU\\s*" * FIXED_POINT_REAL *
    raw"s\sWALL\s\(\s*([+-]?[0-9]+)\scalls\)",
    "i",
)
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L35-L36.
const INIT_RUN_TIME_BLOCK = r"Called by (?<head>init_run):(?<body>\X+?)^\s*$"im
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L53-L54.
const ELECTRONS_TIME_BLOCK = r"Called by (?<head>electrons):(?<body>\X+?)^\s*$"im
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L78-L79.
const C_BANDS_TIME_BLOCK = r"Called by (?<head>c_bands):(?<body>\X+?)^\s*$"im
const SUM_BAND_TIME_BLOCK = r"Called by (?<head>sum_band):(?<body>\X+?)^\s*$"im
const EGTERG_TIME_BLOCK = r"Called by (?<head>\*egterg):(?<body>\X+?)^\s*$"im
const H_PSI_TIME_BLOCK = r"Called by (?<head>h_psi):(?<body>\X+?)^\s*$"im
const GENERAL_ROUTINES_TIME_BLOCK = r"(?<head>General routines)(?<body>\X+?)^\s*$"im
const PARALLEL_ROUTINES_TIME_BLOCK = r"(?<head>Parallel routines)(?<body>\X+?)^\s*$"im
const TERMINATED_DATE = r"This run was terminated on:(.+)"i  # TODO: Date
const JOB_DONE = r"JOB DONE\."i
