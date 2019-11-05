using Test

using Compat: isnothing

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Namelists.CP
using QuantumESPRESSOBase.Namelists.PHonon

using QuantumESPRESSOParsers.InputParsers.Namelists

@testset "Parse empty string" begin
    @test isnothing(parse(PWscf.ControlNamelist, " "))
    @test isnothing(parse(PWscf.ControlNamelist, "&control/"))
    @test parse(PWscf.ControlNamelist, "&control\n/\n") == PWscf.ControlNamelist()
    @test_throws Meta.ParseError parse(PWscf.ControlNamelist, " ")
    @test_logs(
        (:info, "An empty Namelist found! Default values will be used!"),
        parse(PWscf.ControlNamelist, "&control\n/"),
    )
end # testset

@testset "Parse CP input" begin
    # This data is from https://github.com/QEF/q-e/blob/master/CPV/examples/example01/run_example.
    str = raw"""
    &control
        calculation='cp',
        restart_mode='from_scratch',
        nstep=20, iprint=20, isave=20,
        dt=5.0,
        ndr=90, ndw=91,
        pseudo_dir='$PSEUDO_DIR/',
        outdir='$TMP_DIR/',
    /
    &system
        ibrav=8, celldm(1)=9.28990, celldm(2)=1.73206, celldm(3)=1.09955,
        nat=18, ntyp=2, nbnd=48, nspin=1,
        ecutwfc=20.0, ecutrho=150.0,
        nr1b=16, nr2b=16, nr3b=16,
        qcutz=150., q2sigma=2.0, ecfixed=16.0,
    /
    &electrons
        electron_dynamics='damp', electron_damping=0.2,
        startingwfc='random', ampre=0.01,
        emass=700., emass_cutoff=3.,
    /
    &ions
        ion_dynamics='none',
        ion_radius(1)=1.0, ion_radius(2)=1.0,
    /
    ATOMIC_SPECIES
    O  16.00 O.pz-rrkjus.UPF
    Si  28.00 Si.pz-vbc.UPF
    ATOMIC_POSITIONS bohr
    O  3.18829368  14.83237039   1.22882961
    O  7.83231469   6.78704039   1.22882961
    O  2.07443467   5.99537992   4.73758250
    O  6.72031366  14.04231898   4.73758250
    O  3.96307134  11.26989826   7.87860582
    O  8.60802134   3.22295920   7.87860582
    O  3.96307134   4.81915267   9.14625133
    O  8.60802134  12.86448267   9.14625133
    O  3.18736469   1.25668055   5.58029607
    O  7.83324368   9.30201055   5.58029607
    O  2.07536366  10.09206195   2.07358613
    O  6.71938467   2.04673195   2.07358613
    Si  0.28891589   8.04533000   3.40456284
    Si  4.93386589   0.00000000   3.40456284
    Si  2.13389003  12.27717358  -0.04188031
    Si  6.77884003   4.23184358  -0.04188031
    Si  2.13389003   3.81348642   6.85202747
    Si  6.77884003  11.85881642   6.85202747
    """
    control = parse(CP.ControlNamelist, str)
    system = parse(CP.SystemNamelist, str)
    electrons = parse(CP.ElectronsNamelist, str)
    ions = parse(CP.IonsNamelist, str)
    # @test control == CP.ControlNamelist(
    #     calculation = "cp",
    #     restart_mode = "from_scratch",
    #     nstep = 20,
    #     iprint = 20,
    #     isave = 20,
    #     dt = 5.0,
    #     ndr = 90,
    #     ndw = 91,
    #     pseudo_dir = raw"$PSEUDO_DIR/",
    #     outdir = raw"$TMP_DIR/",
    # )
    # @test system == CP.SystemNamelist(
    #     ibrav = 8,
    #     celldm = [9.28990, 1.73206, 1.09955],
    #     nat = 18,
    #     ntyp = 2,
    #     nbnd = 48,
    #     nspin = 1,
    #     ecutwfc = 20.0,
    #     ecutrho = 150.0,
    #     nr1b = 16,
    #     nr2b = 16,
    #     nr3b = 16,
    #     qcutz = 150.0,
    #     q2sigma = 2.0,
    #     ecfixed = 16.0,
    # )
    # @test electrons = CP.ElectronsNamelist(
    #     electron_dynamics = "damp",
    #     electron_damping = 0.2,
    #     startingwfc = "random",
    #     ampre = 0.01,
    #     emass = 700.0,
    #     emass_cutoff = 3.0,
    # )
    # @test ions == CP.IonsNamelist(ion_dynamics = "none", ion_radius = [1.0, 1.0])
end # testset
