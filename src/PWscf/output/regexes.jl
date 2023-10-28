using ReadableRegex: ANY, lazy_zero_or_more, capture, either, maybe, @rs_str

# See https://gist.github.com/singularitti/e9e04c501ddfe40ba58917a754707b2e
const INTEGER = raw"([-+]?[0-9]+)"
const FIXED_POINT_REAL = raw"([-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)"
const GENERAL_REAL = raw"([-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[eE][-+]?[0-9]+)?)"
const EQUAL_SIGN = raw"\s*=\s*"

# This format is from https://github.com/QEF/q-e/blob/4132a64/Modules/environment.f90#L215-L224.
const PARALLEL_INFO = r"(?<kind>(?:Parallel version [^,]*|Serial version))(?:, running on\s*(?<num>[0-9]+) processors)?"
const READING_INPUT_FROM = r"(?:Reading input from \s*(.*|standard input))"
const PWSCF_VERSION = r"Program PWSCF v\.(?<version>[0-9]\.[0-9]+\.?[0-9]?)"
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L374-L375.
const FFT_DIMENSIONS = Regex(
    "Dense  grid:\\s*$INTEGER\\s+G-vectors\\s+FFT dimensions: \\(\\s*$INTEGER,\\s*$INTEGER,\\s*$INTEGER\\)",
)
# The following format is from https://github.com/QEF/q-e/blob/7357cdb/PW/src/summary.f90#L100-L119.
const SUMMARY_BLOCK = r"(bravais-lattice index\X+?)\s*celldm"  # Match between "bravais-lattice index" & the 1st of the "celldm"s, `+?` means un-greedy matching (required)
# 'bravais-lattice index     = ',I12
const BRAVAIS_LATTICE_INDEX = Regex("bravais-lattice index$EQUAL_SIGN$INTEGER")
# 'lattice parameter (alat)  = ',F12.4,'  a.u.'
const LATTICE_PARAMETER = Regex("lattice parameter \\(alat\\)$EQUAL_SIGN$FIXED_POINT_REAL")
# 'unit-cell volume          = ',F12.4,' (a.u.)^3'
const UNIT_CELL_VOLUME = Regex("unit-cell volume$EQUAL_SIGN$FIXED_POINT_REAL")
# 'number of atoms/cell      = ',I12
const NUMBER_OF_ATOMS_PER_CELL = Regex("number of atoms\\/cell$EQUAL_SIGN$INTEGER")
# 'number of atomic types    = ',I12
const NUMBER_OF_ATOMIC_TYPES = Regex("number of atomic types$EQUAL_SIGN$INTEGER")
# 'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')'
const NUMBER_OF_ELECTRONS = Regex(
    "number of electrons$EQUAL_SIGN$FIXED_POINT_REAL" *
    "(?:\\(up:\\s*$FIXED_POINT_REAL, down:\\s*$FIXED_POINT_REAL\\))?",
)
# 'number of Kohn-Sham states= ',I12
const NUMBER_OF_KOHN_SHAM_STATES = Regex("number of Kohn-Sham states$EQUAL_SIGN$INTEGER")
# 'kinetic-energy cutoff     = ',F12.4,'  Ry'
const KINETIC_ENERGY_CUTOFF = Regex(
    "kinetic-energy cutoff$EQUAL_SIGN$FIXED_POINT_REAL\\s+Ry"
)
# 'charge density cutoff     = ',F12.4,'  Ry'
const CHARGE_DENSITY_CUTOFF = Regex(
    "charge density cutoff$EQUAL_SIGN$FIXED_POINT_REAL\\s+Ry"
)
# 'cutoff for Fock operator  = ',F12.4,'  Ry'
const CUTOFF_FOR_FOCK_OPERATOR = Regex(
    "cutoff for Fock operator$EQUAL_SIGN$FIXED_POINT_REAL\\s+Ry"
)
# 'convergence threshold     = ',1PE12.1
const CONVERGENCE_THRESHOLD = Regex("convergence threshold$EQUAL_SIGN$GENERAL_REAL")
# 'mixing beta               = ',0PF12.4
const MIXING_BETA = Regex("mixing beta$EQUAL_SIGN$FIXED_POINT_REAL")
# 'number of iterations used = ',I12,2X,A,' mixing'
const NUMBER_OF_ITERATIONS_USED = Regex(
    "number of iterations used$EQUAL_SIGN$INTEGER\\s+([-+\\w]+)\\s+mixing"
)
const EXCHANGE_CORRELATION = r"Exchange-correlation\s*=\s*(.*)"
# "nstep                     = ",I12
const NSTEP = Regex("nstep$EQUAL_SIGN$INTEGER")
# The following format is from https://github.com/QEF/q-e/blob/4132a64/Modules/fft_base.f90#L70-L91.
const FFT_BASE_INFO = r"""\s*(?<head>Parallelization info|G-vector sticks info)
\s*--------------------
\s*sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
(?<body>(?:\s*Min.*)?(?:\s*Max.*)?\s*Sum.*)"""m
const SYM_OPS = r"\h*(No symmetry found|(?<n>\d+)\s*Sym\. Ops\..*found)"
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L341-L381.
const K_POINTS_BLOCK = r"""
number of k points=\s*(?<nk>[0-9]+)\h*(?<metainfo>.*)
\s*(?:cart\. coord\. in units 2pi\/alat\s*(?<cartesian>\X+?)^\s*$|Number of k-points >= 100: set verbosity='high' to print them\.)
\s*(?:cryst\. coord\.\s*(?<crystal>\X+?)
\s*Dense  grid)?"""m
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L353-L354.
# '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)'
const K_POINTS_ITEM = Regex(
    "k\\(.*\\) = \\(\\s*$FIXED_POINT_REAL\\s*$FIXED_POINT_REAL\\s*$FIXED_POINT_REAL\\s*\\), wk =\\s*$FIXED_POINT_REAL",
)
const STRESS_BLOCK = r"""
^[ \t]*
total\s+stress\s*\(Ry\/bohr\*\*3\)\s+
\(kbar\)\s+P=\s*([\-|\+]? (?: [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?)
\R
(
(?:
\s*
(?:
[\-|\+]? (?: [0-9]*[\.][0-9]+ | [0-9]+[\.]?[0-9]*)
    ([E|e|d|D][+|-]?[0-9]+)?[ \t]*
){6}
){3}
)
"""mx
const SELF_CONSISTENT_CALCULATION_BLOCK =
    rs"Self-consistent Calculation" *
    capture(lazy_zero_or_more(ANY)) *
    rs"End of self-consistent calculation"
const ITERATION_BLOCK = r"(?<=iteration #)(.*?)(?=iteration #|End of self-consistent calculation)"s
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L920-L921.
# '     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2
const ITERATION_HEAD = Regex(
    "\\s*$INTEGER\\s+ecut=\\s*$FIXED_POINT_REAL\\s+Ry\\s+beta=\\s*$FIXED_POINT_REAL"
)
# These formats are from https://github.com/QEF/q-e/blob/4132a64/PW/src/c_bands.f90#L129-L130
# and https://github.com/QEF/q-e/blob/4132a64/PW/src/c_bands.f90#L65-L73.
const C_BANDS = Regex(
    """
    (?<diag>Davidson diagonalization.*|CG style diagonalization|PPCG style diagonalization)
    \\h*ethr =\\h*$GENERAL_REAL,  avg # of iterations =\\h*$FIXED_POINT_REAL""",
    "m",
)
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L917-L918.
# '     total cpu time spent up to now is ',F10.1,' secs'
const TOTAL_CPU_TIME = Regex(
    "total cpu time spent up to now is\\s*$FIXED_POINT_REAL\\s* secs"
)
const KS_ENERGIES_BLOCK = r"""
(Number\s+of\s+k-points\s+>=\s+100:\s+set\s+verbosity='high'\s+to\s+print\s+the\s+bands\.
|
\s*(------\s+SPIN\s+UP\s+------------|------\s+SPIN\s+DOWN\s+----------)?
(?:
\s*(?:k\s+=(?:\s*[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*){3}\s*\(\s*(?:[-+]?[0-9]+)\s*PWs\)\s+bands\s+\(ev\):
|
k\s+=(?:\s*[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*){3}\s+band\s+energies\s+\(ev\):
)
(?:(?:\s*[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*){1,8})+
)+
)"""mx
const KS_ENERGIES_BANDS = r"""
\s*k\s+=\s*(?<k>.*)\s*\(\s*(?<PWs>[-+]?[0-9]+)\s*PWs\)\s+bands\s+\(ev\):
(?<band>(?:\s*[-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[eE][-+]?[0-9]+)?)+)"""
const KS_ENERGIES_BAND_ENERGIES = r"""
\s*k\s+=\s*(?<k>.*)\s*band\s+energies\s+\(ev\):
(?<band>(?:\s*[-+]?(?:[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*)(?:[eE][-+]?[0-9]+)?)+)"""
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L1257-L1261.
const UNCONVERGED_ELECTRONS_ENERGY = Regex(
    """
    ^[^!]\\s+total energy\\s+=\\s*$FIXED_POINT_REAL\\s+Ry(?:\\s*Harris-Foulkes estimate\\s+=\\s*$FIXED_POINT_REAL\\s+Ry)?
    \\s*estimated scf accuracy\\s+<\\s*$GENERAL_REAL\\s+Ry""",
    "m",
)
const CONVERGED_ELECTRONS_ENERGY = Regex(
    """
    ^!\\h+total energy\\s+=\\s*$FIXED_POINT_REAL\\s+Ry(?:\\s*Harris-Foulkes estimate\\s+=\\s*$FIXED_POINT_REAL\\s+Ry)?(?:\\s*total all-electron energy\\s+=\\s*$FIXED_POINT_REAL\\s+Ry)?
    \\s*estimated scf accuracy\\s+<\\s*$GENERAL_REAL\\s+Ry
    \\s*(?<ae>total all-electron energy =.*Ry)?\\s*(?<decomp>The total energy is the sum of the following terms:
    \\s*one-electron contribution =.*Ry
    \\s*hartree contribution      =.*Ry
    \\s*xc contribution           =.*Ry
    \\s* ewald contribution        =.*Ry)?\\s*(?<one>one-center paw contrib.*Ry
    \\s*-> PAW hartree energy AE =.*Ry
    \\s*-> PAW hartree energy PS =.*Ry
    \\s*-> PAW xc energy AE      =.*Ry
    \\s*-> PAW xc energy PS      =.*Ry
    \\s*-> total E_H with PAW    =.*Ry
    \\s*-> total E_XC with PAW   =.*Ry)?\\s*(?<smearing>smearing contrib.*Ry)?""",
    "m",
)
const HOURS_MINUTES = rs"(\d+)h\s*(\d+)m"
const MINUTES_SECONDS = rs"(\d+)m\s*(\d+\.\d{2})s"
const SECONDS = rs"(\d+)\.(\d{2})s"
const TIME_FORMAT = capture(either(HOURS_MINUTES, MINUTES_SECONDS, SECONDS))
const TIMED_ITEM = Regex(
    rs"([\w:]+)\s*:\s*" *
    TIME_FORMAT *
    rs"\s*CPU\s*" *
    TIME_FORMAT *
    rs"\s*WALL" *
    maybe(rs"\s*\(\s*(\d+)\s*calls\)"),
)
const TERMINATED_DATE = r"This run was terminated on:(.+)"  # TODO: Date
const JOB_DONE = r"JOB DONE\."
# These formats are from https://github.com/QEF/q-e/blob/4132a64/UtilXlib/error_handler.f90#L48-L68.
const ERROR_BLOCK = r"%{78}(?<body>\X+?)\s*%{78}"
const ERROR_IN_ROUTINE = r"Error in routine\s+(.*)\s+\((.*)\):"
