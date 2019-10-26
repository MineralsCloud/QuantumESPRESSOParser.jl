using Test

using Compat: isnothing
using DataFrames
using QuantumESPRESSOBase
using QuantumESPRESSOParsers.OutputParsers.PWscf

@testset "Parse scf output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/example01/reference/al.scf.ppcg.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test parse_head(str) == Dict(
        "number of electrons" => 3.0,
        "number of Kohn-Sham states" => 6,
        "lattice parameter (alat)" => 7.5,
        "number of iterations used" => 8,
        "convergence threshold" => 1.0e-6,
        "unit-cell volume" => 105.4688,
        "Exchange-correlation" => "SLA PZ NOGX NOGC ( 1  1  0  0 0 0)",
        "kinetic-energy cutoff" => 15.0,
        "bravais-lattice index" => 2,
        "number of atoms/cell" => 1,
        "number of atomic types" => 1,
        "mixing beta" => 0.7,
        "charge density cutoff" => 60.0,
    )

    parallelization_info = parse_parallelization_info(str) ## TODO : compare two tables

    @test parse_ibz(str) == (
        [
            0.0625 0.0625 0.0625 0.0078125; 
            0.0625 0.0625 0.1875 0.0234375; 
            0.0625 0.0625 0.3125 0.0234375; 
            0.0625 0.0625 0.4375 0.0234375; 
            0.0625 0.0625 0.5625 0.0234375; 
            0.0625 0.0625 0.6875 0.0234375; 
            0.0625 0.0625 0.8125 0.0234375; 
            0.0625 0.0625 0.9375 0.0234375; 
            0.0625 0.1875 0.1875 0.0234375; 
            0.0625 0.1875 0.3125 0.046875; 
            0.0625 0.1875 0.4375 0.046875; 
            0.0625 0.1875 0.5625 0.046875; 
            0.0625 0.1875 0.6875 0.046875; 
            0.0625 0.1875 0.8125 0.046875; 
            0.0625 0.1875 0.9375 0.046875; 
            0.0625 0.3125 0.3125 0.0234375; 
            0.0625 0.3125 0.4375 0.046875; 
            0.0625 0.3125 0.5625 0.046875; 
            0.0625 0.3125 0.6875 0.046875; 
            0.0625 0.3125 0.8125 0.046875; 
            0.0625 0.3125 0.9375 0.046875; 
            0.0625 0.4375 0.4375 0.0234375; 
            0.0625 0.4375 0.5625 0.046875; 
            0.0625 0.4375 0.6875 0.046875; 
            0.0625 0.4375 0.8125 0.046875; 
            0.0625 0.4375 0.9375 0.046875; 
            0.0625 0.5625 0.5625 0.0234375; 
            0.0625 0.5625 0.6875 0.046875; 
            0.0625 0.5625 0.8125 0.046875; 
            0.0625 0.6875 0.6875 0.0234375; 
            0.0625 0.6875 0.8125 0.046875; 
            0.0625 0.8125 0.8125 0.0234375; 
            0.1875 0.1875 0.1875 0.0078125; 
            0.1875 0.1875 0.3125 0.0234375; 
            0.1875 0.1875 0.4375 0.0234375; 
            0.1875 0.1875 0.5625 0.0234375; 
            0.1875 0.1875 0.6875 0.0234375; 
            0.1875 0.1875 0.8125 0.0234375; 
            0.1875 0.3125 0.3125 0.0234375; 
            0.1875 0.3125 0.4375 0.046875; 
            0.1875 0.3125 0.5625 0.046875; 
            0.1875 0.3125 0.6875 0.046875; 
            0.1875 0.3125 0.8125 0.046875; 
            0.1875 0.4375 0.4375 0.0234375; 
            0.1875 0.4375 0.5625 0.046875; 
            0.1875 0.4375 0.6875 0.046875; 
            0.1875 0.4375 0.8125 0.046875; 
            0.1875 0.5625 0.5625 0.0234375; 
            0.1875 0.5625 0.6875 0.046875; 
            0.1875 0.6875 0.6875 0.0234375; 
            0.3125 0.3125 0.3125 0.0078125; 
            0.3125 0.3125 0.4375 0.0234375; 
            0.3125 0.3125 0.5625 0.0234375; 
            0.3125 0.3125 0.6875 0.0234375; 
            0.3125 0.4375 0.4375 0.0234375; 
            0.3125 0.4375 0.5625 0.046875; 
            0.3125 0.4375 0.6875 0.046875; 
            0.3125 0.5625 0.5625 0.0234375; 
            0.4375 0.4375 0.4375 0.0078125; 
            0.4375 0.4375 0.5625 0.0234375], 
        nothing
    )


    @test parse_stress(str) == (
        [-17.35],
        [[
          -0.00011796 -0.0 -0.0
          -0.0 -0.00011796 0.0
          -0.0 0.0 -0.00011796
        ]],
        [[
          -17.35 -0.0 -0.0
          -0.0 -17.35 0.0
          -0.0 0.0 -17.35
        ]],
    )

    @test isempty(parse_cell_parameters(str))

    @test isempty(parse_atomic_positions(str))

    @test parse_scf_calculation(str) == [[
        Dict(
            "total energy" => -4.1872583,
            "time" => 0.5,
            "iteration" => 1,
            "Harris-Foulkes estimate" => -4.18806959,
            "ecut" => 15.0,
            "estimated scf accuracy" => 0.00586766,
            "beta" => 0.7,
        ),
        Dict(
            "total energy" => -4.18723304,
            "time" => 0.6,
            "iteration" => 2,
            "Harris-Foulkes estimate" => -4.18727083,
            "ecut" => 15.0,
            "estimated scf accuracy" => 0.00049876,
            "beta" => 0.7,
        ),
        Dict(
            "total energy" => -4.18725801,
            "time" => 0.7,
            "iteration" => 3,
            "Harris-Foulkes estimate" => -4.18725803,
            "ecut" => 15.0,
            "estimated scf accuracy" => 3.31e-6,
            "beta" => 0.7,
        ),
        Dict(
            "total energy" => -4.18725743,
            "time" => 1.0,
            "iteration" => 4,
            "Harris-Foulkes estimate" => -4.18725817,
            "ecut" => 15.0,
            "estimated scf accuracy" => 4.08e-6,
            "beta" => 0.7,
        ),
        Dict("time" => 1.1, "iteration" => 5, "ecut" => 15.0, "beta" => 0.7),
    ]]

    @test parse_total_energy(str) == [-4.18725747]

    @test parse_version(str) == "6.3"

    @test parse_processors_num(str) == ("Parallel version (MPI)", 4)

    @test parse_fft_dimensions(str) == [869, 15, 15, 15]

    clock = parse_clock(str) ## TODO : compare two tables

    @test isjobdone(str) == true
end

@testset "Parse vc-relax output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/VCSexample/reference/As.vcs00.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test parse_head(str) == Dict(
        "number of electrons" => 10.0,
        "number of Kohn-Sham states" => 9,
        "lattice parameter (alat)" => 7.0103,
        "number of iterations used" => 8,
        "nstep" => 55,
        "convergence threshold" => 1.0e-7,
        "unit-cell volume" => 245.3705,
        "Exchange-correlation" => "SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)",
        "kinetic-energy cutoff" => 25.0,
        "bravais-lattice index" => 0,
        "number of atoms/cell" => 2,
        "number of atomic types" => 1,
        "mixing beta" => 0.7,
        "charge density cutoff" => 100.0,
    )

# parallelization_info = parse_parallelization_info(vc) ## TODO : compare two tables

    @test parse_ibz(str) == (
        [
            0.0 0.0 0.1534638 0.0625; 
            -0.1436461 -0.2488023 0.2557731 0.1875; 
            0.2872922 0.4976046 -0.0511547 0.1875; 
            0.1436461 0.2488023 0.0511546 0.1875; 
            -0.2872922 0.0 0.3580823 0.1875; 
            0.1436461 0.746407 0.0511546 0.375; 
            0.0 0.4976046 0.1534638 0.375; 
            0.5745844 0.0 -0.2557731 0.1875; 
            0.0 0.0 0.4603915 0.0625; 
            0.4309383 0.746407 0.1534638 0.1875], 
        nothing
    )

    @test parse_stress(str) == (
        [
         217.55,
         142.98,
         33.14,
         -1.09,
         -35.67,
         -43.98,
         -23.49,
         -9.45,
         6.88,
         -3.22,
         -4.13,
         -4.46,
         -4.07,
         1.09,
         0.02,
         0.08,
         0.34,
         0.43,
         0.1,
        ],
        [
         [0.001725 -0.0 -0.0; -0.0 0.001725 -0.0; -0.0 0.0 0.00098671],
         [0.00107816 -0.0 -0.0; 0.0 0.00107816 0.0; -0.0 0.0 0.00075965],
         [0.00014848 0.0 -0.0; -0.0 0.00014848 -0.0; -0.0 0.0 0.0003788],
         [-9.215e-5 0.0 0.0; 0.0 -9.215e-5 0.0; 0.0 0.0 0.00016204],
         [-0.00029487 0.0 0.0; 0.0 -0.00029487 -0.0; 0.0 -0.0 -0.0001376],
         [-0.0003037 0.0 0.0; 0.0 -0.0003037 -0.0; 0.0 -0.0 -0.00028951],
         [-0.00017615 0.0 0.0; -0.0 -0.00017615 0.0; 0.0 0.0 -0.0001267],
         [-5.562e-5 -0.0 0.0; 0.0 -5.562e-5 0.0; 0.0 0.0 -8.156e-5],
         [8.346e-5 0.0 -0.0; -0.0 8.346e-5 0.0; -0.0 -0.0 -2.668e-5],
         [-5.63e-6 0.0 0.0; -0.0 -5.63e-6 -0.0; 0.0 0.0 -5.446e-5],
         [-1.622e-5 -0.0 0.0; 0.0 -1.622e-5 0.0; 0.0 0.0 -5.179e-5],
         [-2.303e-5 0.0 0.0; -0.0 -2.303e-5 0.0; 0.0 0.0 -4.485e-5],
         [-2.463e-5 -0.0 0.0; -0.0 -2.463e-5 0.0; 0.0 0.0 -3.382e-5],
         [4.29e-6 0.0 -0.0; -0.0 4.29e-6 -0.0; -0.0 -0.0 1.371e-5],
         [-1.15e-6 -0.0 0.0; 0.0 -1.15e-6 -0.0; 0.0 0.0 2.67e-6],
         [-5.6e-7 -0.0 0.0; 0.0 -5.6e-7 -0.0; 0.0 -0.0 2.7e-6],
         [1.8e-6 0.0 -0.0; 0.0 1.8e-6 -0.0; -0.0 -0.0 3.37e-6],
         [2.06e-6 0.0 -0.0; 0.0 2.06e-6 -0.0; -0.0 0.0 4.65e-6],
         [3.4e-7 -0.0 -0.0; -0.0 3.4e-7 -0.0; -0.0 -0.0 1.34e-6],
        ],
        [
         [253.76 -0.0 -0.0; -0.0 253.76 -0.0; -0.0 0.0 145.15],
         [158.6 -0.0 -0.0; 0.0 158.6 0.0; -0.0 0.0 111.75],
         [21.84 0.0 -0.0; -0.0 21.84 -0.0; -0.0 0.0 55.72],
         [-13.56 0.0 0.0; 0.0 -13.56 0.0; 0.0 0.0 23.84],
         [-43.38 0.0 0.0; 0.0 -43.38 -0.0; 0.0 -0.0 -20.24],
         [-44.68 0.0 0.0; 0.0 -44.68 -0.0; 0.0 -0.0 -42.59],
         [-25.91 0.0 0.0; -0.0 -25.91 0.0; 0.0 0.0 -18.64],
         [-8.18 -0.0 0.0; 0.0 -8.18 0.0; 0.0 0.0 -12.0],
         [12.28 0.0 -0.0; -0.0 12.28 0.0; -0.0 -0.0 -3.93],
         [-0.83 0.0 0.0; -0.0 -0.83 -0.0; 0.0 0.0 -8.01],
         [-2.39 -0.0 0.0; 0.0 -2.39 0.0; 0.0 0.0 -7.62],
         [-3.39 0.0 0.0; -0.0 -3.39 0.0; 0.0 0.0 -6.6],
         [-3.62 -0.0 0.0; -0.0 -3.62 0.0; 0.0 0.0 -4.97],
         [0.63 0.0 -0.0; -0.0 0.63 -0.0; -0.0 -0.0 2.02],
         [-0.17 -0.0 0.0; 0.0 -0.17 -0.0; 0.0 0.0 0.39],
         [-0.08 -0.0 0.0; 0.0 -0.08 -0.0; 0.0 -0.0 0.4],
         [0.26 0.0 -0.0; 0.0 0.26 -0.0; -0.0 -0.0 0.5],
         [0.3 0.0 -0.0; 0.0 0.3 -0.0; -0.0 0.0 0.68],
         [0.05 -0.0 -0.0; -0.0 0.05 -0.0; -0.0 -0.0 0.2],
        ],
    )

    @test parse_cell_parameters(str) == Array{Float64,2}[
        [
         4.134122414288868 -0.0 5.764073531453344
         -2.06705987518055 3.5802547047057334 5.764073545474017
         -2.06705987518055 -3.5802547047057334 5.764073545474017
        ],
        [
         4.258010239861242 -0.0 5.880745786192249
         -2.1290037879667367 3.6875446899718898 5.880745828254267
         -2.1290037879667367 -3.6875446899718898 5.880745828254267
        ],
        [
         4.257935776069806 -0.0 6.029049617385969
         -2.1289665490606833 3.6874802018889103 6.029049687489331
         -2.1289665490606833 -3.6874802018889103 6.029049687489331
        ],
        [
         4.256160990277466 -0.0 6.201972838075738
         -2.128079152659345 3.685943192680819 6.201972929210109
         -2.128079152659345 -3.685943192680819 6.201972929210109
        ],
        [
         4.243072377070925 -0.0 6.378994392459937
         -2.1215348425509064 3.674608123947842 6.378994504625316
         -2.1215348425509064 -3.674608123947842 6.378994504625316
        ],
        [
         4.218745718184832 -0.0 6.23404739066659
         -2.1093715025923556 3.6535406194474285 6.234047453759617
         -2.1093715025923556 -3.6535406194474285 6.234047453759617
        ],
        [
         4.1850966230457125 -0.0 6.224906508101251
         -2.0925469480124597 3.6243996472070634 6.22490656418394
         -2.0925469480124597 -3.6243996472070634 6.22490656418394
        ],
        [
         4.148104649075318 -0.0 6.2099855121799035
         -2.0740509540169265 3.5923636584758123 6.209985568262593
         -2.0740509540169265 -3.5923636584758123 6.209985568262593
        ],
        [
         4.168067247472185 -0.0 6.193465044954656
         -2.0840322567205276 3.6096517715389167 6.19346509402701
         -2.0840322567205276 -3.6096517715389167 6.19346509402701
        ],
        [
         4.167806336778376 -0.0 6.172953375993247
         -2.0839018013736235 3.6094258213918877 6.172953425065601
         -2.0839018013736235 -3.6094258213918877 6.172953425065601
        ],
        [
         4.1667403660822515 -0.0 6.1487371369131125
         -2.083368812520393 3.60850266526544 6.148737185985466
         -2.083368812520393 -3.60850266526544 6.148737185985466
        ],
        [
         4.164514163709053 -0.0 6.121379713385375
         -2.0822557113337936 3.6065747106368105 6.121379755447392
         -2.0822557113337936 -3.6065747106368105 6.121379755447392
        ],
        [
         4.161039210143204 -0.0 6.0916974310226735
         -2.0805182345508695 3.603565313499996 6.0916974730846905
         -2.0805182345508695 -3.603565313499996 6.0916974730846905
        ],
        [
         4.1616956159657645 -0.0 6.100511673981714
         -2.080846437462149 3.6041337816648866 6.100511709033395
         -2.080846437462149 -3.6041337816648866 6.100511709033395
        ],
        [
         4.161633953048285 -0.0 6.100722159327019
         -2.0808156060034095 3.6040803769234864 6.100722201389036
         -2.0808156060034095 -3.6040803769234864 6.100722201389036
        ],
        [
         4.161542320943423 -0.0 6.101146319720615
         -2.0807697934561467 3.6040010269276994 6.101146361782632
         -2.0807697934561467 -3.6040010269276994 6.101146361782632
        ],
        [
         4.161612171933618 -0.0 6.101837630007264
         -2.0808047189512444 3.6040615121086916 6.101837665058945
         -2.0808047189512444 -3.6040615121086916 6.101837665058945
        ],
        [
         4.161722479574197 -0.0 6.102896366036401
         -2.080859869266366 3.6041570419604976 6.102896408098418
         -2.080859869266366 -3.6041570419604976 6.102896408098418
        ],
        [
         4.161722479574197 -0.0 6.102896366036401
         -2.080859869266366 3.6041570419604976 6.102896408098418
         -2.080859869266366 -3.6041570419604976 6.102896408098418
        ],
    ]

     # @test parse_atomic_positions(vc) ==
     # QuantumESPRESSOBase.Cards.AtomicPositionsCard[
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.288384589, 0.288384588, 0.288384588], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.288384589, 0.288384588, 0.288384588], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.284847404, 0.284847401, 0.284847401], [1, 1, 1]),
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.284847404, 0.284847401, 0.284847401], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.280292903, 0.280292897, 0.280292897], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.280292903, 0.280292897, 0.280292897], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275309566, 0.275309559, 0.275309559], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275309566, 0.275309559, 0.275309559], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275518788, 0.275518781, 0.275518781], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275518788, 0.275518781, 0.275518781], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275471859, 0.275471852, 0.275471852], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275471859, 0.275471852, 0.275471852], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275290955, 0.275290948, 0.275290948], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.275290955, 0.275290948, 0.275290948], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.27492874, 0.274928733, 0.274928733], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.27492874, 0.274928733, 0.274928733], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.274339437, 0.274339429, 0.274339429], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.274339437, 0.274339429, 0.274339429], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.273610438, 0.27361043, 0.27361043], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.273610438, 0.27361043, 0.27361043], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272799426, 0.272799417, 0.272799417], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272799426, 0.272799417, 0.272799417], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.271967078, 0.271967069, 0.271967069], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.271967078, 0.271967069, 0.271967069], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.27199446, 0.271994451, 0.271994451], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.27199446, 0.271994451, 0.271994451], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272035264, 0.272035254, 0.272035254], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272035264, 0.272035254, 0.272035254], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272093332, 0.272093322, 0.272093322], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272093332, 0.272093322, 0.272093322], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272162275, 0.272162266, 0.272162266], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272162275, 0.272162266, 0.272162266], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272235517, 0.272235507, 0.272235507], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272235517, 0.272235507, 0.272235507], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272233124, 0.272233115, 0.272233115], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272233124, 0.272233115, 0.272233115], [1, 1, 1])
     #           ]
     #      ), 
     #      QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}(
     #           "crystal", QuantumESPRESSOBase.Cards.AtomicPosition[
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272233124, 0.272233115, 0.272233115], [1, 1, 1]), 
     #                QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("As", [0.272233124, 0.272233115, 0.272233115], [1, 1, 1])
     #           ]
     #      )
     # ]

    @test parse_scf_calculation(str) == [
        [
         Dict("time" => 0.3, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.4, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.4, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.5, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.5, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 0.8, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.9, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 0.9, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.0, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.0, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.1, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.1, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 1.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.5, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.6, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.6, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.7, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 1.7, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 2.0, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.1, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.1, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.2, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.3, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.3, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.3, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 2.6, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.7, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.7, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.8, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.9, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.9, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 2.9, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 3.2, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.3, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.4, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.4, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.5, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.5, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.6, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 3.8, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 3.9, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.0, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.0, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.1, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.1, "iteration" => 6, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.2, "iteration" => 7, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 4.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.5, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.6, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 4.7, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 4.9, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.0, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.0, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.1, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.1, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 5.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.5, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.6, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 5.6, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 5.9, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.0, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.0, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.1, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.1, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 6.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.6, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 6.6, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 6.9, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.0, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.0, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.1, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.1, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 7.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.5, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.6, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 7.6, "iteration" => 5, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 8.0, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 8.1, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 8.1, "iteration" => 3, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 8.1, "iteration" => 4, "ecut" => 25.0, "beta" => 0.7),
        ],
        [
         Dict("time" => 8.4, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7),
         Dict("time" => 8.5, "iteration" => 2, "ecut" => 25.0, "beta" => 0.7),
        ],
        [Dict("time" => 8.8, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7)],
        [Dict("time" => 9.1, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7)],
        [Dict("time" => 9.5, "iteration" => 1, "ecut" => 25.0, "beta" => 0.7)],
    ]

    @test parse_total_energy(str) == [
        -25.4401674,
        -25.46010129,
        -25.48262993,
        -25.49296522,
        -25.49696246,
        -25.49568749,
        -25.49793805,
        -25.49857642,
        -25.49870559,
        -25.49906146,
        -25.4992724,
        -25.49940546,
        -25.49943809,
        -25.49945971,
        -25.49946354,
        -25.4994652,
        -25.4994663,
        -25.49946664,
        -25.49946688,
    ]

    @test parse_version(str) == "6.0"

    @test parse_processors_num(str) == ("Parallel version (MPI)", 2)

    @test parse_fft_dimensions(str) == [4159, 24, 24, 24]

    clock = parse_clock(str) ## TODO : compare two tables

    @test isjobdone(str) == true
end