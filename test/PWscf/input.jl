using Test

using Crystallography
using QuantumESPRESSOBase.PWscf
using QuantumESPRESSOParser

@testset "Parse empty strings" begin
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
    @test tryparse(ControlNamelist, " ") === nothing
    @test tryparse(ControlNamelist, "&control/") == ControlNamelist()
    @test parse(ControlNamelist, "&control\n/\n") == ControlNamelist()
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
end

@testset "Test parsing a `PWInput`: silicon" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example01/run_example#L81-L116
    str = read("../data/si.scf.cg.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(;
        prefix = "silicon",
        tstress = true,
        tprnfor = true,
        pseudo_dir = "pseudo/",
        outdir = "./",
    )
    @test pw.system ==
          SystemNamelist(; ibrav = 2, celldm = [10.20], nat = 2, ntyp = 1, ecutwfc = 18)
    @test pw.electrons == ElectronsNamelist(; diagonalization = "cg", conv_thr = 1e-8)
    @test pw.ions == IonsNamelist()
    @test pw.cell == CellNamelist()
    @test pw.atomic_species ==
          AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
    @test pw.atomic_positions == AtomicPositionsCard([
        AtomicPosition("Si", [0.0, 0.0, 0.0]),
        AtomicPosition("Si", [0.25, 0.25, 0.25]),
    ])
    @test pw.k_points == SpecialPointsCard(
        [
            0.125 0.125 0.125 1.0
            0.125 0.125 0.375 3.0
            0.125 0.125 0.625 3.0
            0.125 0.125 0.875 3.0
            0.125 0.375 0.375 3.0
            0.125 0.375 0.625 6.0
            0.125 0.375 0.875 6.0
            0.125 0.625 0.625 3.0
            0.375 0.375 0.375 1.0
            0.375 0.375 0.625 3.0
        ],
    )
    @test pw.cell_parameters === nothing
end

@testset "Test parsing a `PWInput`: silicon bands" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example01/run_example#L125-L172
    str = read("../data/si.band.cg.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(;
        calculation = "bands",
        pseudo_dir = "pseudo/",
        outdir = "./",
        prefix = "silicon",
    )
    @test pw.system == SystemNamelist(;
        ibrav = 2,
        celldm = [10.2],
        nat = 2,
        ntyp = 1,
        ecutwfc = 18.0,
        nbnd = 8,
    )
    @test pw.electrons == ElectronsNamelist(; diagonalization = "cg")
    @test pw.atomic_species ==
          AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
    @test pw.atomic_positions == AtomicPositionsCard([
        AtomicPosition("Si", [0.0, 0.0, 0.0]),
        AtomicPosition("Si", [0.25, 0.25, 0.25]),
    ])
    @test pw.k_points == SpecialPointsCard(
        [
            0.0 0.0 0.0 1.0
            0.0 0.0 0.1 1.0
            0.0 0.0 0.2 1.0
            0.0 0.0 0.3 1.0
            0.0 0.0 0.4 1.0
            0.0 0.0 0.5 1.0
            0.0 0.0 0.6 1.0
            0.0 0.0 0.7 1.0
            0.0 0.0 0.8 1.0
            0.0 0.0 0.9 1.0
            0.0 0.0 1.0 1.0
            0.0 0.0 0.0 1.0
            0.0 0.1 0.1 1.0
            0.0 0.2 0.2 1.0
            0.0 0.3 0.3 1.0
            0.0 0.4 0.4 1.0
            0.0 0.5 0.5 1.0
            0.0 0.6 0.6 1.0
            0.0 0.7 0.7 1.0
            0.0 0.8 0.8 1.0
            0.0 0.9 0.9 1.0
            0.0 1.0 1.0 1.0
            0.0 0.0 0.0 1.0
            0.1 0.1 0.1 1.0
            0.2 0.2 0.2 1.0
            0.3 0.3 0.3 1.0
            0.4 0.4 0.4 1.0
            0.5 0.5 0.5 1.0
        ],
    )
end

@testset "Test parsing a `PWInput`: aluminium" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example01/run_example#L186-L268
    str = read("../data/al.scf.david.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(;
        calculation = "scf",
        restart_mode = "from_scratch",
        pseudo_dir = "pseudo/",
        outdir = "./",
        prefix = "al",
        tprnfor = true,
        tstress = true,
    )
    @test pw.system == SystemNamelist(;
        ibrav = 2,
        celldm = [7.50],
        nat = 1,
        ntyp = 1,
        ecutwfc = 15.0,
        occupations = "smearing",
        smearing = "marzari-vanderbilt",
        degauss = 0.05,
    )
    @test pw.electrons == ElectronsNamelist()
    @test pw.atomic_species ==
          AtomicSpeciesCard([AtomicSpecies("Al", 26.98, "Al.pz-vbc.UPF")])
    @test pw.atomic_positions ==
          AtomicPositionsCard([AtomicPosition("Al", [0.0, 0.0, 0.0])])
    @test pw.k_points == SpecialPointsCard(
        [
            0.0625000 0.0625000 0.0625000 1.00
            0.0625000 0.0625000 0.1875000 3.00
            0.0625000 0.0625000 0.3125000 3.00
            0.0625000 0.0625000 0.4375000 3.00
            0.0625000 0.0625000 0.5625000 3.00
            0.0625000 0.0625000 0.6875000 3.00
            0.0625000 0.0625000 0.8125000 3.00
            0.0625000 0.0625000 0.9375000 3.00
            0.0625000 0.1875000 0.1875000 3.00
            0.0625000 0.1875000 0.3125000 6.00
            0.0625000 0.1875000 0.4375000 6.00
            0.0625000 0.1875000 0.5625000 6.00
            0.0625000 0.1875000 0.6875000 6.00
            0.0625000 0.1875000 0.8125000 6.00
            0.0625000 0.1875000 0.9375000 6.00
            0.0625000 0.3125000 0.3125000 3.00
            0.0625000 0.3125000 0.4375000 6.00
            0.0625000 0.3125000 0.5625000 6.00
            0.0625000 0.3125000 0.6875000 6.00
            0.0625000 0.3125000 0.8125000 6.00
            0.0625000 0.3125000 0.9375000 6.00
            0.0625000 0.4375000 0.4375000 3.00
            0.0625000 0.4375000 0.5625000 6.00
            0.0625000 0.4375000 0.6875000 6.00
            0.0625000 0.4375000 0.8125000 6.00
            0.0625000 0.4375000 0.9375000 6.00
            0.0625000 0.5625000 0.5625000 3.00
            0.0625000 0.5625000 0.6875000 6.00
            0.0625000 0.5625000 0.8125000 6.00
            0.0625000 0.6875000 0.6875000 3.00
            0.0625000 0.6875000 0.8125000 6.00
            0.0625000 0.8125000 0.8125000 3.00
            0.1875000 0.1875000 0.1875000 1.00
            0.1875000 0.1875000 0.3125000 3.00
            0.1875000 0.1875000 0.4375000 3.00
            0.1875000 0.1875000 0.5625000 3.00
            0.1875000 0.1875000 0.6875000 3.00
            0.1875000 0.1875000 0.8125000 3.00
            0.1875000 0.3125000 0.3125000 3.00
            0.1875000 0.3125000 0.4375000 6.00
            0.1875000 0.3125000 0.5625000 6.00
            0.1875000 0.3125000 0.6875000 6.00
            0.1875000 0.3125000 0.8125000 6.00
            0.1875000 0.4375000 0.4375000 3.00
            0.1875000 0.4375000 0.5625000 6.00
            0.1875000 0.4375000 0.6875000 6.00
            0.1875000 0.4375000 0.8125000 6.00
            0.1875000 0.5625000 0.5625000 3.00
            0.1875000 0.5625000 0.6875000 6.00
            0.1875000 0.6875000 0.6875000 3.00
            0.3125000 0.3125000 0.3125000 1.00
            0.3125000 0.3125000 0.4375000 3.00
            0.3125000 0.3125000 0.5625000 3.00
            0.3125000 0.3125000 0.6875000 3.00
            0.3125000 0.4375000 0.4375000 3.00
            0.3125000 0.4375000 0.5625000 6.00
            0.3125000 0.4375000 0.6875000 6.00
            0.3125000 0.5625000 0.5625000 3.00
            0.4375000 0.4375000 0.4375000 1.00
            0.4375000 0.4375000 0.5625000 3.00
        ],
    )
end

@testset "Test parsing a `PWInput`: cobalt relaxation" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example02/run_example#L83-L112
    str = read("../data/co.rx.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(;
        calculation = "relax",
        pseudo_dir = "pseudo/",
        outdir = "./",
        prefix = "CO",
    )
    @test pw.system ==
          SystemNamelist(; ibrav = 0, nat = 2, ntyp = 2, ecutwfc = 24, ecutrho = 144)
    @test pw.electrons == ElectronsNamelist(; conv_thr = 1e-7)
    @test pw.ions == IonsNamelist()
    @test pw.atomic_species == AtomicSpeciesCard([
        AtomicSpecies("O", 1, "O.pz-rrkjus.UPF"),
        AtomicSpecies("C", 1.0, "C.pz-rrkjus.UPF"),
    ])
    @test pw.atomic_positions == AtomicPositionsCard(
        [AtomicPosition("C", [2.256, 0.0, 0.0]), AtomicPosition("O", [0, 0, 0], [0, 0, 0])],
        "bohr",
    )
    @test pw.k_points == GammaPointCard()
end

@testset "Test parsing a `PWInput`: FeO" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example08/run_example#L84-L124
    str = read("../data/feo_LDA.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(;
        calculation = "scf",
        restart_mode = "from_scratch",
        prefix = "feo_af",
        pseudo_dir = "pseudo/",
        outdir = "./",
        tstress = true,
        tprnfor = true,
    )
    @test pw.system == SystemNamelist(;
        ibrav = 0,
        celldm = [8.19],
        nat = 4,
        ntyp = 3,
        ecutwfc = 30.0,
        ecutrho = 240,
        nbnd = 20,
        starting_magnetization = [0, 0.5, -0.5],
        occupations = "smearing",
        smearing = "mv",
        degauss = 0.01,
        nspin = 2,
        # lda_plus_u = true,
        # Hubbard_U = [nothing, 1e-8, 1e-8],
    )
    @test pw.electrons ==
          ElectronsNamelist(; conv_thr = 1e-6, mixing_beta = 0.3, mixing_fixed_ns = 0)
    @test pw.ions == IonsNamelist()
    @test pw.cell_parameters == CellParametersCard([
        0.50 0.50 1.00
        0.50 1.00 0.50
        1.00 0.50 0.50
    ])
    @test pw.atomic_species == AtomicSpeciesCard([
        AtomicSpecies("O1", 1, "O.pz-rrkjus.UPF"),
        AtomicSpecies("Fe1", 1.0, "Fe.pz-nd-rrkjus.UPF"),
        AtomicSpecies("Fe2", 1.0, "Fe.pz-nd-rrkjus.UPF"),
    ])
    @test pw.atomic_positions == AtomicPositionsCard(
        [
            AtomicPosition("O1", [1, 1, 1] / 4),
            AtomicPosition("O1", [3, 3, 3] / 4),
            AtomicPosition("Fe1", [0, 0, 0]),
            AtomicPosition("Fe2", [0.5, 0.5, 0.5]),
        ],
        "crystal",
    )
    @test pw.k_points == KMeshCard(MonkhorstPackGrid([2, 2, 2], [0, 0, 0]))
end
