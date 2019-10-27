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

    @test parse_summary(str) == Dict(
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

    @test parse_fft_base_info(str) == groupby(
        DataFrame(
            [
             "sticks" "Min" 30 30 10
             "sticks" "Max" 31 31 11
             "sticks" "Sum" 121 121 43
             "gvecs" "Min" 216 216 45
             "gvecs" "Max" 218 218 46
             "gvecs" "Sum" 869 869 181
            ],
            [:kind, :stats, :dense, :smooth, :PW],
        ),
        :kind,
    )

    @test parse_ibz(str) == (
        [
         0.0625 0.0625 0.0625 0.0078125
         0.0625 0.0625 0.1875 0.0234375
         0.0625 0.0625 0.3125 0.0234375
         0.0625 0.0625 0.4375 0.0234375
         0.0625 0.0625 0.5625 0.0234375
         0.0625 0.0625 0.6875 0.0234375
         0.0625 0.0625 0.8125 0.0234375
         0.0625 0.0625 0.9375 0.0234375
         0.0625 0.1875 0.1875 0.0234375
         0.0625 0.1875 0.3125 0.046875
         0.0625 0.1875 0.4375 0.046875
         0.0625 0.1875 0.5625 0.046875
         0.0625 0.1875 0.6875 0.046875
         0.0625 0.1875 0.8125 0.046875
         0.0625 0.1875 0.9375 0.046875
         0.0625 0.3125 0.3125 0.0234375
         0.0625 0.3125 0.4375 0.046875
         0.0625 0.3125 0.5625 0.046875
         0.0625 0.3125 0.6875 0.046875
         0.0625 0.3125 0.8125 0.046875
         0.0625 0.3125 0.9375 0.046875
         0.0625 0.4375 0.4375 0.0234375
         0.0625 0.4375 0.5625 0.046875
         0.0625 0.4375 0.6875 0.046875
         0.0625 0.4375 0.8125 0.046875
         0.0625 0.4375 0.9375 0.046875
         0.0625 0.5625 0.5625 0.0234375
         0.0625 0.5625 0.6875 0.046875
         0.0625 0.5625 0.8125 0.046875
         0.0625 0.6875 0.6875 0.0234375
         0.0625 0.6875 0.8125 0.046875
         0.0625 0.8125 0.8125 0.0234375
         0.1875 0.1875 0.1875 0.0078125
         0.1875 0.1875 0.3125 0.0234375
         0.1875 0.1875 0.4375 0.0234375
         0.1875 0.1875 0.5625 0.0234375
         0.1875 0.1875 0.6875 0.0234375
         0.1875 0.1875 0.8125 0.0234375
         0.1875 0.3125 0.3125 0.0234375
         0.1875 0.3125 0.4375 0.046875
         0.1875 0.3125 0.5625 0.046875
         0.1875 0.3125 0.6875 0.046875
         0.1875 0.3125 0.8125 0.046875
         0.1875 0.4375 0.4375 0.0234375
         0.1875 0.4375 0.5625 0.046875
         0.1875 0.4375 0.6875 0.046875
         0.1875 0.4375 0.8125 0.046875
         0.1875 0.5625 0.5625 0.0234375
         0.1875 0.5625 0.6875 0.046875
         0.1875 0.6875 0.6875 0.0234375
         0.3125 0.3125 0.3125 0.0078125
         0.3125 0.3125 0.4375 0.0234375
         0.3125 0.3125 0.5625 0.0234375
         0.3125 0.3125 0.6875 0.0234375
         0.3125 0.4375 0.4375 0.0234375
         0.3125 0.4375 0.5625 0.046875
         0.3125 0.4375 0.6875 0.046875
         0.3125 0.5625 0.5625 0.0234375
         0.4375 0.4375 0.4375 0.0078125
         0.4375 0.4375 0.5625 0.0234375
        ],
        nothing,
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

    # @test parse_scf_calculation(str) == groupby(
    #     DataFrame(
    #         [
    #          1 1 15.0 "PPCG style diagonalization" 0.01 8.1 0.7 0.5 -4.1872583 -4.18806959 0.00586766; 
    #          1 2 15.0 "PPCG style diagonalization" 0.000196 2.4 0.7 0.6 -4.18723304 -4.18727083 0.00049876; 
    #          1 3 15.0 "PPCG style diagonalization" 1.66e-5 4.5 0.7 0.7 -4.18725801 -4.18725803 3.31e-6; 
    #          1 4 15.0 "PPCG style diagonalization" 1.1e-7 5.9 0.7 1.0 -4.18725743 -4.18725817 4.08e-6; 
    #          1 5 15.0 "PPCG style diagonalization" 1.1e-7 4.0 0.7 1.1 nothing nothing nothing
    #         ],
    #         [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ]
    #     ),
    #     :n,
    # )

    @test parse_total_energy(str) == [-4.18725747]

    @test parse_version(str) == "6.3"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 4)

    @test parse_fft_dimensions(str) == (869, (nr1 = 15, nr2 = 15, nr3 = 15))

    @test parse_clock(str) == groupby(
        DataFrame(
            [
             "" "init_run" 0.03 0.03 1
             "" "electrons" 0.89 1.06 1
             "" "forces" 0.0 0.0 1
             "" "stress" 0.01 0.01 1
             "init_run" "wfcinit" 0.02 0.02 1
             "init_run" "potinit" 0.0 0.0 1
             "init_run" "hinit0" 0.0 0.0 1
             "electrons" "c_bands" 0.84 0.99 6
             "electrons" "sum_band" 0.04 0.06 6
             "electrons" "v_of_rho" 0.0 0.0 6
             "electrons" "mix_rho" 0.0 0.0 6
             "c_bands" "init_us_2" 0.01 0.01 900
             "c_bands" "ppcg_k" 0.74 0.88 360
             "c_bands" "wfcrot" 0.1 0.12 300
             "h_psi" "h_psi:pot" 0.44 0.52 1930
             "h_psi" "h_psi:calbec" 0.03 0.02 1930
             "h_psi" "vloc_psi" 0.39 0.48 1930
             "h_psi" "add_vuspsi" 0.01 0.01 1930
             "General routines" "calbec" 0.02 0.02 2230
             "General routines" "fft" 0.0 0.0 24
             "General routines" "ffts" 0.0 0.0 6
             "General routines" "fftw" 0.37 0.45 22818
             "Parallel routines" "fft_scatt_xy" 0.04 0.08 22848
             "Parallel routines" "fft_scatt_yz" 0.12 0.15 22848
            ],
            [:subroutine, :item, :CPU, :wall, :calls],
        ),
        :subroutine,
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == false

    @test isjobdone(str) == true

    @test haserror(str) == false
end

@testset "Parse vc-relax output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/VCSexample/reference/As.vcs00.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test parse_summary(str) == Dict(
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

    @test parse_fft_base_info(str) == groupby(
        DataFrame(
            [
             "sticks" "Min" 174 174 60
             "sticks" "Max" 175 175 61
             "sticks" "Sum" 349 349 121
             "gvecs" "Min" 2079 2079 416
             "gvecs" "Max" 2080 2080 417
             "gvecs" "Sum" 4159 4159 833
            ],
            [:kind, :stats, :dense, :smooth, :PW],
        ),
        :kind,
    )

    @test parse_ibz(str) == (
        [
         0.0 0.0 0.1534638 0.0625
         -0.1436461 -0.2488023 0.2557731 0.1875
         0.2872922 0.4976046 -0.0511547 0.1875
         0.1436461 0.2488023 0.0511546 0.1875
         -0.2872922 0.0 0.3580823 0.1875
         0.1436461 0.746407 0.0511546 0.375
         0.0 0.4976046 0.1534638 0.375
         0.5745844 0.0 -0.2557731 0.1875
         0.0 0.0 0.4603915 0.0625
         0.4309383 0.746407 0.1534638 0.1875
        ],
        nothing,
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

    # @test parse_scf_calculation(str) == groupby(
    #     DataFrame(
    #         [
    #          1 1 25.0 "Davidson diagonalization with overlap" 0.01 4.3 0.7 0.3 nothing nothing nothing;
    #          1 2 25.0 "Davidson diagonalization with overlap" 0.000156 1.0 0.7 0.4 nothing nothing nothing;
    #          1 3 25.0 "Davidson diagonalization with overlap" 8.89e-6 1.6 0.7 0.4 nothing nothing nothing;
    #          1 4 25.0 "Davidson diagonalization with overlap" 5.01e-8 3.1 0.7 0.5 nothing nothing nothing;
    #          1 5 25.0 "Davidson diagonalization with overlap" 7.74e-9 1.4 0.7 0.5 nothing nothing nothing;
    #          2 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 4.2 0.7 0.8 nothing nothing nothing;
    #          2 2 25.0 "Davidson diagonalization with overlap" 8.24e-6 3.1 0.7 0.9 nothing nothing nothing;
    #          2 3 25.0 "Davidson diagonalization with overlap" 6.78e-6 1.0 0.7 0.9 nothing nothing nothing;
    #          2 4 25.0 "Davidson diagonalization with overlap" 1.49e-6 1.0 0.7 1.0 nothing nothing nothing;
    #          2 5 25.0 "Davidson diagonalization with overlap" 4.63e-7 2.6 0.7 1.0 nothing nothing nothing;
    #          2 6 25.0 "Davidson diagonalization with overlap" 1.04e-8 2.1 0.7 1.1 nothing nothing nothing;
    #          2 7 25.0 "Davidson diagonalization with overlap" 1.89e-9 1.1 0.7 1.1 nothing nothing nothing;
    #          3 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 4.9 0.7 1.4 nothing nothing nothing;
    #          3 2 25.0 "Davidson diagonalization with overlap" 2.69e-5 3.1 0.7 1.5 nothing nothing nothing;
    #          3 3 25.0 "Davidson diagonalization with overlap" 2.43e-5 1.0 0.7 1.5 nothing nothing nothing;
    #          3 4 25.0 "Davidson diagonalization with overlap" 5.68e-6 1.0 0.7 1.6 nothing nothing nothing;
    #          3 5 25.0 "Davidson diagonalization with overlap" 1.88e-6 2.2 0.7 1.6 nothing nothing nothing;
    #          3 6 25.0 "Davidson diagonalization with overlap" 6.49e-8 2.5 0.7 1.7 nothing nothing nothing;
    #          3 7 25.0 "Davidson diagonalization with overlap" 4.28e-9 1.8 0.7 1.7 nothing nothing nothing;
    #          4 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 5.1 0.7 2.0 nothing nothing nothing;
    #          4 2 25.0 "Davidson diagonalization with overlap" 3.7e-6 3.0 0.7 2.1 nothing nothing nothing;
    #          4 3 25.0 "Davidson diagonalization with overlap" 2.05e-6 1.0 0.7 2.1 nothing nothing nothing;
    #          4 4 25.0 "Davidson diagonalization with overlap" 3.6e-7 1.5 0.7 2.2 nothing nothing nothing;
    #          4 5 25.0 "Davidson diagonalization with overlap" 1.95e-8 3.0 0.7 2.3 nothing nothing nothing;
    #          4 6 25.0 "Davidson diagonalization with overlap" 2.0e-9 1.0 0.7 2.3 nothing nothing nothing;
    #          4 7 25.0 "Davidson diagonalization with overlap" 1.16e-9 1.5 0.7 2.3 nothing nothing nothing;
    #          5 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 5.3 0.7 2.6 nothing nothing nothing;
    #          5 2 25.0 "Davidson diagonalization with overlap" 5.22e-6 3.0 0.7 2.7 nothing nothing nothing;
    #          5 3 25.0 "Davidson diagonalization with overlap" 2.41e-6 1.0 0.7 2.7 nothing nothing nothing;
    #          5 4 25.0 "Davidson diagonalization with overlap" 3.56e-7 1.3 0.7 2.8 nothing nothing nothing;
    #          5 5 25.0 "Davidson diagonalization with overlap" 2.63e-8 3.0 0.7 2.9 nothing nothing nothing;
    #          5 6 25.0 "Davidson diagonalization with overlap" 2.83e-9 1.0 0.7 2.9 nothing nothing nothing;
    #          5 7 25.0 "Davidson diagonalization with overlap" 1.44e-9 2.0 0.7 2.9 nothing nothing nothing;
    #          6 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 4.2 0.7 3.2 nothing nothing nothing;
    #          6 2 25.0 "Davidson diagonalization with overlap" 1.41e-6 3.0 0.7 3.3 nothing nothing nothing;
    #          6 3 25.0 "Davidson diagonalization with overlap" 1.41e-6 1.0 0.7 3.4 nothing nothing nothing;
    #          6 4 25.0 "Davidson diagonalization with overlap" 6.85e-7 1.0 0.7 3.4 nothing nothing nothing;
    #          6 5 25.0 "Davidson diagonalization with overlap" 8.37e-8 3.0 0.7 3.5 nothing nothing nothing;
    #          6 6 25.0 "Davidson diagonalization with overlap" 4.29e-9 1.0 0.7 3.5 nothing nothing nothing;
    #          6 7 25.0 "Davidson diagonalization with overlap" 1.47e-9 2.0 0.7 3.6 nothing nothing nothing;
    #          7 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 4.2 0.7 3.8 nothing nothing nothing;
    #          7 2 25.0 "Davidson diagonalization with overlap" 6.2e-6 3.0 0.7 3.9 nothing nothing nothing;
    #          7 3 25.0 "Davidson diagonalization with overlap" 5.67e-6 1.0 0.7 4.0 nothing nothing nothing;
    #          7 4 25.0 "Davidson diagonalization with overlap" 1.15e-6 1.0 0.7 4.0 nothing nothing nothing;
    #          7 5 25.0 "Davidson diagonalization with overlap" 2.62e-7 3.0 0.7 4.1 nothing nothing nothing;
    #          7 6 25.0 "Davidson diagonalization with overlap" 9.51e-9 1.0 0.7 4.1 nothing nothing nothing;
    #          7 7 25.0 "Davidson diagonalization with overlap" 4.98e-9 2.0 0.7 4.2 nothing nothing nothing;
    #          8 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.7 0.7 4.4 nothing nothing nothing;
    #          8 2 25.0 "Davidson diagonalization with overlap" 1.9e-6 3.0 0.7 4.5 nothing nothing nothing;
    #          8 3 25.0 "Davidson diagonalization with overlap" 1.54e-6 1.0 0.7 4.5 nothing nothing nothing;
    #          8 4 25.0 "Davidson diagonalization with overlap" 2.67e-7 1.0 0.7 4.6 nothing nothing nothing;
    #          8 5 25.0 "Davidson diagonalization with overlap" 7.53e-8 2.7 0.7 4.7 nothing nothing nothing;
    #          9 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.8 0.7 4.9 nothing nothing nothing;
    #          9 2 25.0 "Davidson diagonalization with overlap" 2.41e-6 3.0 0.7 5.0 nothing nothing nothing;
    #          9 3 25.0 "Davidson diagonalization with overlap" 2.15e-6 1.0 0.7 5.0 nothing nothing nothing;
    #          9 4 25.0 "Davidson diagonalization with overlap" 3.58e-7 1.0 0.7 5.1 nothing nothing nothing;
    #          9 5 25.0 "Davidson diagonalization with overlap" 1.12e-7 2.7 0.7 5.1 nothing nothing nothing;
    #          10 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.6 0.7 5.4 nothing nothing nothing;
    #          10 2 25.0 "Davidson diagonalization with overlap" 3.74e-7 3.0 0.7 5.5 nothing nothing nothing;
    #          10 3 25.0 "Davidson diagonalization with overlap" 1.98e-7 1.0 0.7 5.5 nothing nothing nothing;
    #          10 4 25.0 "Davidson diagonalization with overlap" 2.5e-8 1.0 0.7 5.6 nothing nothing nothing;
    #          10 5 25.0 "Davidson diagonalization with overlap" 4.09e-9 3.0 0.7 5.6 nothing nothing nothing;
    #          11 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.4 0.7 5.9 nothing nothing nothing;
    #          11 2 25.0 "Davidson diagonalization with overlap" 1.03e-7 2.7 0.7 6.0 nothing nothing nothing;
    #          11 3 25.0 "Davidson diagonalization with overlap" 8.21e-8 1.0 0.7 6.0 nothing nothing nothing;
    #          11 4 25.0 "Davidson diagonalization with overlap" 1.21e-8 2.2 0.7 6.1 nothing nothing nothing;
    #          11 5 25.0 "Davidson diagonalization with overlap" 2.13e-9 1.7 0.7 6.1 nothing nothing nothing;
    #          12 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.5 0.7 6.4 nothing nothing nothing;
    #          12 2 25.0 "Davidson diagonalization with overlap" 1.3e-7 3.0 0.7 6.5 nothing nothing nothing;
    #          12 3 25.0 "Davidson diagonalization with overlap" 1.25e-7 1.0 0.7 6.6 nothing nothing nothing;
    #          12 4 25.0 "Davidson diagonalization with overlap" 2.09e-8 2.0 0.7 6.6 nothing nothing nothing;
    #          13 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.6 0.7 6.9 nothing nothing nothing;
    #          13 2 25.0 "Davidson diagonalization with overlap" 1.51e-7 3.0 0.7 7.0 nothing nothing nothing;
    #          13 3 25.0 "Davidson diagonalization with overlap" 1.51e-7 1.0 0.7 7.0 nothing nothing nothing;
    #          13 4 25.0 "Davidson diagonalization with overlap" 3.25e-8 1.4 0.7 7.1 nothing nothing nothing;
    #          13 5 25.0 "Davidson diagonalization with overlap" 1.36e-9 3.0 0.7 7.1 nothing nothing nothing;
    #          14 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 3.1 0.7 7.4 nothing nothing nothing;
    #          14 2 25.0 "Davidson diagonalization with overlap" 1.52e-7 3.0 0.7 7.5 nothing nothing nothing;
    #          14 3 25.0 "Davidson diagonalization with overlap" 1.52e-7 1.0 0.7 7.5 nothing nothing nothing;
    #          14 4 25.0 "Davidson diagonalization with overlap" 3.4e-8 1.0 0.7 7.6 nothing nothing nothing;
    #          14 5 25.0 "Davidson diagonalization with overlap" 9.49e-9 2.7 0.7 7.6 nothing nothing nothing;
    #          15 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 2.2 0.7 8.0 nothing nothing nothing;
    #          15 2 25.0 "Davidson diagonalization with overlap" 1.27e-8 3.0 0.7 8.1 nothing nothing nothing;
    #          15 3 25.0 "Davidson diagonalization with overlap" 1.27e-8 1.0 0.7 8.1 nothing nothing nothing;
    #          15 4 25.0 "Davidson diagonalization with overlap" 2.69e-9 1.0 0.7 8.1 nothing nothing nothing;
    #          16 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 1.2 0.7 8.4 nothing nothing nothing;
    #          16 2 25.0 "Davidson diagonalization with overlap" 1.06e-9 2.0 0.7 8.5 nothing nothing nothing;
    #          17 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 1.6 0.7 8.8 nothing nothing nothing;
    #          18 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 1.9 0.7 9.1 nothing nothing nothing;
    #          19 1 25.0 "Davidson diagonalization with overlap" 1.0e-6 1.0 0.7 9.5 nothing nothing nothing
    #         ],
    #         [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
    #     ),
    #     :n,
    # )

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

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 2)

    @test parse_fft_dimensions(str) == (4159, (nr1 = 24, nr2 = 24, nr3 = 24))

    @test parse_clock(str) == groupby(
        DataFrame(
            [
             "" "init_run" 0.12 0.15 1
             "" "electrons" 5.8 5.9 19
             "" "update_pot" 1.09 1.11 18
             "" "forces" 0.45 0.45 19
             "" "stress" 0.94 0.94 19
             "init_run" "wfcinit" 0.02 0.02 1
             "init_run" "potinit" 0.02 0.02 1
             "electrons" "c_bands" 4.94 5.0 96
             "electrons" "sum_band" 0.73 0.76 96
             "electrons" "v_of_rho" 0.07 0.08 109
             "electrons" "mix_rho" 0.02 0.02 96
             "c_bands" "init_us_2" 0.07 0.12 2310
             "c_bands" "cegterg" 4.84 4.91 960
             "*egterg" "h_psi" 3.79 3.77 3133
             "*egterg" "g_psi" 0.01 0.02 2163
             "*egterg" "cdiaghg" 0.39 0.42 2883
             "h_psi" "h_psi:pot" 3.77 3.76 3133
             "h_psi" "h_psi:calbec" 0.12 0.1 3133
             "h_psi" "vloc_psi" 3.54 3.57 3133
             "h_psi" "add_vuspsi" 0.1 0.08 3133
             "General routines" "calbec" 0.14 0.13 4083
             "General routines" "fft" 0.05 0.07 542
             "General routines" "fftw" 3.82 3.88 55066
             "General routines" "davcio" 0.0 0.0 10
             "Parallel routines" "fft_scatter" 0.55 0.57 55608
            ],
            [:subroutine, :item, :CPU, :wall, :calls],
        ),
        :subroutine,
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true

    @test haserror(str) == false
end