using Test

using QuantumESPRESSOBase.Inputs.PWscf
using QuantumESPRESSOParser.Inputs

@testset "Parse empty strings" begin
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
    @test tryparse(ControlNamelist, " ") === nothing
    @test tryparse(ControlNamelist, "&control/") == ControlNamelist()
    @test parse(ControlNamelist, "&control\n/\n") == ControlNamelist()
    @test_throws Meta.ParseError parse(ControlNamelist, " ")
end

@testset "Test parsing a `PWInput`: silicon" begin
    # From https://github.com/QEF/q-e/blob/121edf5/PW/examples/example01/run_example#L81-L116
    str = read("data/si.scf.cg.in", String)
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
    str = read("data/si.band.cg.in", String)
    pw = parse(PWInput, str)
    @test pw.control == ControlNamelist(
        calculation = "bands",
        pseudo_dir = raw"pseudo/",
        outdir = raw"./",
        prefix = "silicon",
    )
    @test pw.system == SystemNamelist(
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
