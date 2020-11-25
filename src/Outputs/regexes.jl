# See https://gist.github.com/singularitti/e9e04c501ddfe40ba58917a754707b2e
const SIGN = char_in("-+")
const _DIGIT = char_in('0':'9')
const EQUAL = rs"\s*=\s*"
# See https://github.com/jkrumbiegel/ReadableRegex.jl/blob/771869b/README.md
const INTEGER = capture(look_for(maybe(SIGN) * one_or_more(_DIGIT), not_after = '.'))
const REAL = capture(maybe(SIGN) * maybe(zero_or_more(_DIGIT) * '.') * one_or_more(_DIGIT))
const EXP_REAL = capture(REAL * maybe(char_in("eE") * maybe(SIGN) * one_or_more(_DIGIT)))

# This format is from https://github.com/QEF/q-e/blob/4132a64/Modules/environment.f90#L215-L224.
const PARALLEL_INFO = r"(?<kind>(?:Parallel version [^,]*|Serial version))(?:, running on\s*(?<num>[0-9]+) processors)?"
const READING_INPUT_FROM = r"(?:Reading input from \s*(.*|standard input))"
const PWSCF_VERSION = r"Program PWSCF v\.(?<version>[0-9]\.[0-9]+\.?[0-9]?)"
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L374-L375.
const FFT_DIMENSIONS = r"Dense  grid:(.*)G-vectors     FFT dimensions: \((.*)\)"
# The following format is from https://github.com/QEF/q-e/blob/7357cdb/PW/src/summary.f90#L100-L119.
const SUMMARY_BLOCK = r"(bravais-lattice index\X+?)\s*celldm"  # Match between "bravais-lattice index" & the 1st of the "celldm"s, `+?` means un-greedy matching (required)
# 'bravais-lattice index     = ',I12
const BRAVAIS_LATTICE_INDEX = "bravais-lattice index" * EQUAL * INTEGER
# 'lattice parameter (alat)  = ',F12.4,'  a.u.'
const LATTICE_PARAMETER = rs"lattice parameter \(alat\)" * EQUAL * REAL
# 'unit-cell volume          = ',F12.4,' (a.u.)^3'
const UNIT_CELL_VOLUME = "unit-cell volume" * EQUAL * REAL
# 'number of atoms/cell      = ',I12
const NUMBER_OF_ATOMS_PER_CELL = "number of atoms/cell" * EQUAL * INTEGER
# 'number of atomic types    = ',I12
const NUMBER_OF_ATOMIC_TYPES = "number of atomic types" * EQUAL * INTEGER
# 'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')'
NUMBER_OF_ELECTRONS =
    "number of electrons" *
    EQUAL *
    REAL *
    maybe(rs"\(up:\s*" * REAL * rs", down:\s*" * REAL * ')')
# 'number of Kohn-Sham states= ',I12
const NUMBER_OF_KOHN_SHAM_STATES = "number of Kohn-Sham states" * EQUAL * INTEGER
# 'kinetic-energy cutoff     = ',F12.4,'  Ry'
const KINETIC_ENERGY_CUTOFF = "kinetic-energy cutoff" * EQUAL * REAL * rs"\s+Ry"
# 'charge density cutoff     = ',F12.4,'  Ry'
const CHARGE_DENSITY_CUTOFF = "charge density cutoff" * EQUAL * REAL * rs"\s+Ry"
# 'cutoff for Fock operator  = ',F12.4,'  Ry'
const CUTOFF_FOR_FOCK_OPERATOR = "cutoff for Fock operator" * EQUAL * REAL * rs"\s+Ry"
# 'convergence threshold     = ',1PE12.1
const CONVERGENCE_THRESHOLD = "convergence threshold" * EQUAL * EXP_REAL
# 'mixing beta               = ',0PF12.4
const MIXING_BETA = "mixing beta" * EQUAL * REAL
# 'number of iterations used = ',I12,2X,A,' mixing'
const NUMBER_OF_ITERATIONS_USED =
    "number of iterations used" * EQUAL * INTEGER * rs"\s+([-+\w]+)\s+mixing"
const EXCHANGE_CORRELATION = r"Exchange-correlation\s*=\s*(.*)"
# "nstep                     = ",I12
const NSTEP = "nstep" * EQUAL * INTEGER
# The following format is from https://github.com/QEF/q-e/blob/4132a64/Modules/fft_base.f90#L70-L91.
const FFT_BASE_INFO = r"""\s*(?<head>Parallelization info|G-vector sticks info)
\s*--------------------
\s*sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
(?<body>(?:\s*Min.*)?(?:\s*Max.*)?\s*Sum.*)"""m
const SYM_OPS = r"\h*(No symmetry found|(?<n>\d+)\s*Sym\. Ops\..*found)"
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L341-L381.
const K_POINTS_BLOCK = r"""
number of k points=\s*(?<nk>[0-9]+)\h*(?<metainfo>.*)
\s*(?:cart\. coord\. in units 2pi\/alat\s*(?<cart>\X+?)^\s*$|Number of k-points >= 100: set verbosity='high' to print them\.)
\s*(?:cryst\. coord\.\s*(?<cryst>\X+?)
\s*Dense  grid)?"""m
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/summary.f90#L353-L354.
# '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)'
K_POINTS_ITEM =
    rs"k\(.*\) = \(" *
    exactly(3, look_for(INTEGER; before = zero_or_more(WHITESPACE))) *
    rs", wk =\s*" *
    REAL
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/output_tau.f90#L47-L60.
const CELL_PARAMETERS_BLOCK = r"""
CELL_PARAMETERS \h+
\( (?<option>\w+) =? \s* (?<alat>[-+]?[0-9]*\.[0-9]{8})? \) \h*  # Match `alat`: `F12.8`
(?<data>
    (?: \s*
        (?:
            [-+]?[0-9]*\.[0-9]+ \s*  # Match element
        ){3}  # I need exactly 3 elements per vector
    ){3}  # I need exactly 3 vectors
)
"""x
const CELL_PARAMETERS_ITEM = r"""
\s*
([-+]?[0-9]*\.[0-9]+) \s*  # x
([-+]?[0-9]*\.[0-9]+) \s*  # y
([-+]?[0-9]*\.[0-9]+) \s*  # z
"""x
# The following format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/output_tau.f90#L64-L109.
const ATOMIC_POSITIONS_BLOCK = r"""
ATOMIC_POSITIONS \h*                   # Atomic positions start with that string
\( (?<option>\w+) \)                   # Option of the card
(?<data>
    (?:
        \s*
        [A-Za-z]+[A-Za-z0-9]{0,2} \s+  # Atom spec
        (?:
            [-+]?[0-9]*\.[0-9]+ \s*  # Match element
        ){3}                           # I need exactly 3 floats per vector.
        (?:
            [-+]?[0-9]+ \s*
        ){0,3}                     # I need exactly 3 integers in `if_pos`, if there is any.
    )+
)
"""x
const ATOMIC_POSITIONS_ITEM = r"""
\s*
([A-Za-z]+[A-Za-z0-9]{0,2}) \s+  # Atom spec
([-+]?[0-9]*\.[0-9]+) \s*  # x
([-+]?[0-9]*\.[0-9]+) \s*  # y
([-+]?[0-9]*\.[0-9]+) \s*  # z
([-+]?[0-9]+)? \s*            # if_pos(1)
([-+]?[0-9]+)? \s*            # if_pos(2)
([-+]?[0-9]+)? \s*            # if_pos(3)
"""x
const FINAL_COORDINATES_BLOCK = r"""
Begin final coordinates
(\X+?)
End final coordinates
"""
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
const SELF_CONSISTENT_CALCULATION_BLOCK = r"(Self-consistent Calculation\X+?End of self-consistent calculation)"
const ITERATION_BLOCK = r"(?<=iteration #)(.*?)(?=iteration #|End of self-consistent calculation)"s
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L920-L921.
# '     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2
const ITERATION_HEAD = INTEGER * rs"\s+ecut=\s*" * REAL * rs"\s+Ry\s+beta=\s*" * REAL
# These formats are from https://github.com/QEF/q-e/blob/4132a64/PW/src/c_bands.f90#L129-L130
# and https://github.com/QEF/q-e/blob/4132a64/PW/src/c_bands.f90#L65-L73.
C_BANDS =
    capture(
        either(
            rs"Davidson diagonalization.*",
            "CG style diagonalization",
            "PPCG style diagonalization",
        );
        as = "diag",
    ) *
    rs"\h*ethr =\h*" *
    EXP_REAL *
    rs",  avg # of iterations =\h*" *
    REAL
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/electrons.f90#L917-L918.
# '     total cpu time spent up to now is ',F10.1,' secs'
const TOTAL_CPU_TIME = rs"total cpu time spent up to now is\s*" * REAL * rs"\s* secs"
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
# const UNCONVERGED_ELECTRONS_ENERGY = Regex(
#     """
#     ^[^!]\\s+total energy\\s+=\\s*$FIXED_POINT_REAL\\s+Ry
#     \\s*Harris-Foulkes estimate\\s+=\\s*$FIXED_POINT_REAL\\s+Ry
#     \\s*estimated scf accuracy\\s+<\\s*$GENERAL_REAL\\s+Ry""",
#     "m",
# )
# const CONVERGED_ELECTRONS_ENERGY =
#     rs"!\h+total energy\s+=\s*" *
#     FIXED_POINT_REAL *
#     "\\s+Ry" *
#     maybe(rs"\s*Harris-Foulkes estimate\s+=\s*" * FIXED_POINT_REAL * "\\s+Ry") *
#     maybe(rs"\s*estimated scf accuracy\s+<\s*" * GENERAL_REAL * rs"\\s+Ry") *
#     rs"""
#     \s*(?<ae>total all-electron energy =.*Ry)?\s*(?<decomp>The total energy is the sum of the following terms:
#     \s*one-electron contribution =.*Ry
#     \s*hartree contribution      =.*Ry
#     \s*xc contribution           =.*Ry
#     \s* ewald contribution        =.*Ry)?\s*(?<one>one-center paw contrib.*Ry
#     \s*-> PAW hartree energy AE =.*Ry
#     \s*-> PAW hartree energy PS =.*Ry
#     \s*-> PAW xc energy AE      =.*Ry
#     \s*-> PAW xc energy PS      =.*Ry
#     \s*-> total E_H with PAW    =.*Ry
#     \s*-> total E_XC with PAW   =.*Ry)?\s*(?<smearing>smearing contrib.*Ry)?
#     """
# const CONVERGED_ELECTRONS_ENERGY = Regex(
#     """
#     ^!\\h+total energy\\s+=\\s*$FIXED_POINT_REAL\\s+Ry
#     \\s*Harris-Foulkes estimate\\s+=\\s*$FIXED_POINT_REAL\\s+Ry
#     \\s*estimated scf accuracy\\s+<\\s*$GENERAL_REAL\\s+Ry
#     \\s*(?<ae>total all-electron energy =.*Ry)?\\s*(?<decomp>The total energy is the sum of the following terms:
#     \\s*one-electron contribution =.*Ry
#     \\s*hartree contribution      =.*Ry
#     \\s*xc contribution           =.*Ry
#     \\s* ewald contribution        =.*Ry)?\\s*(?<one>one-center paw contrib.*Ry
#     \\s*-> PAW hartree energy AE =.*Ry
#     \\s*-> PAW hartree energy PS =.*Ry
#     \\s*-> PAW xc energy AE      =.*Ry
#     \\s*-> PAW xc energy PS      =.*Ry
#     \\s*-> total E_H with PAW    =.*Ry
#     \\s*-> total E_XC with PAW   =.*Ry)?\\s*(?<smearing>smearing contrib.*Ry)?""",
#     "m",
# )
const TIME_BLOCK = r"(init_run\X+?This run was terminated on:.*)"
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
"""mx
const TIME_ITEM =
    rs"([\w0-9:]+)\s+:\s*" *
    REAL *
    rs"s\sCPU\s*" *
    REAL *
    rs"s\sWALL\s\(\s*" *
    INTEGER *
    rs"\scalls\)"
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L35-L36.
const INIT_RUN_TIME_BLOCK = r"Called by (?<head>init_run):(?<body>\X+?)^\s*$"m
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L53-L54.
const ELECTRONS_TIME_BLOCK = r"Called by (?<head>electrons):(?<body>\X+?)^\s*$"m
# This format is from https://github.com/QEF/q-e/blob/4132a64/PW/src/print_clock_pw.f90#L78-L79.
const C_BANDS_TIME_BLOCK = r"Called by (?<head>c_bands):(?<body>\X+?)^\s*$"m
const SUM_BAND_TIME_BLOCK = r"Called by (?<head>sum_band):(?<body>\X+?)^\s*$"m
const EGTERG_TIME_BLOCK = r"Called by (?<head>\*egterg):(?<body>\X+?)^\s*$"m
const H_PSI_TIME_BLOCK = r"Called by (?<head>h_psi):(?<body>\X+?)^\s*$"m
const GENERAL_ROUTINES_TIME_BLOCK = r"(?<head>General routines)(?<body>\X+?)^\s*$"m
const PARALLEL_ROUTINES_TIME_BLOCK = r"(?<head>Parallel routines)(?<body>\X+?)^\s*$"m
const TERMINATED_DATE = r"This run was terminated on:(.+)"  # TODO: Date
const JOB_DONE = r"JOB DONE\."
# These formats are from https://github.com/QEF/q-e/blob/4132a64/UtilXlib/error_handler.f90#L48-L68.
const ERROR_BLOCK = r"%{78}(?<body>\X+?)\s*%{78}"
const ERROR_IN_ROUTINE = r"Error in routine\s+(.*)\s+\((.*)\):"
