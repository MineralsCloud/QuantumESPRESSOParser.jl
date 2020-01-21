using Test

using Compat: isnothing
using DataFrames
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards.PWscf

using QuantumESPRESSOParsers
using QuantumESPRESSOParsers.Outputs.PWscf

@testset "Parse scf Al output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/example01/reference/al.scf.ppcg.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(
        ibrav = 1,
        alat = 7.5,
        omega = 105.4688,
        nat = 1,
        ntyp = 1,
        nelec = 3.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 6,
        ecutwfc = 15.0,
        ecutrho = 60.0,
        ecutfock = nothing,
        conv_thr = 1.0e-6,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA PZ NOGX NOGC ( 1  1  0  0 0 0)",
        nstep = nothing,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 30 30 10
         "gvecs" "Min" 216 216 45
         "sticks" "Max" 31 31 11
         "gvecs" "Max" 218 218 46
         "sticks" "Sum" 121 121 43
         "gvecs" "Sum" 869 869 181
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (KPointsCard("tpiba", SpecialKPoint[SpecialKPoint([0.0625, 0.0625, 0.0625], 0.0078125), SpecialKPoint([0.0625, 0.0625, 0.1875], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.3125], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.4375], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.5625], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.6875], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.8125], 0.0234375), SpecialKPoint([0.0625, 0.0625, 0.9375], 0.0234375), SpecialKPoint([0.0625, 0.1875, 0.1875], 0.0234375), SpecialKPoint([0.0625, 0.1875, 0.3125], 0.046875)…SpecialKPoint([0.3125, 0.3125, 0.3125], 0.0078125), SpecialKPoint([0.3125, 0.3125, 0.4375], 0.0234375), SpecialKPoint([0.3125, 0.3125, 0.5625], 0.0234375), SpecialKPoint([0.3125, 0.3125, 0.6875], 0.0234375), SpecialKPoint([0.3125, 0.4375, 0.4375], 0.0234375), SpecialKPoint([0.3125, 0.4375, 0.5625], 0.046875), SpecialKPoint([0.3125, 0.4375, 0.6875], 0.046875), SpecialKPoint([0.3125, 0.5625, 0.5625], 0.0234375), SpecialKPoint([0.4375, 0.4375, 0.4375], 0.0078125), SpecialKPoint([0.4375, 0.4375, 0.5625], 0.0234375)]), nothing)

    @test parse_stress(str) == ([-17.35],
        [[
          -0.00011796 -0.0 -0.0
          -0.0 -0.00011796 0.0
          -0.0 0.0 -0.00011796
        ]],
        [[
          -17.35 -0.0 -0.0
          -0.0 -17.35 0.0
          -0.0 0.0 -17.35
        ]],)

    @test isempty(parseall(CellParametersCard, str))

    @test isempty(parseall(AtomicPositionsCard, str))

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 15.0 0.7
             1 2 15.0 0.7
             1 3 15.0 0.7
             1 4 15.0 0.7
             1 5 15.0 0.7
            ],
            [:step, :iteration, :ecut, :β],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-4.18725747,
        -4.18725747,
        2.0e-8,
        nothing,
        [2.93900635, 0.00980673, -1.63461306, -5.50183453],
        nothing,
        0.00037704,)

    @test parse_version(str) == "6.3"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 4)

    @test parse_fft_dimensions(str) == (869, (nr1 = 15, nr2 = 15, nr3 = 15))

    @test parse_bands(str) == ([
         0.0625 0.0625 0.1875
         0.0625 0.3125 0.3125
         0.0625 0.9375 0.5625
         0.0625 0.0625 0.1875
         0.0625 0.4375 0.3125
         0.1875 0.4375 0.6875
         0.0625 0.0625 0.1875
         0.0625 0.4375 0.3125
         0.3125 0.5625 0.8125
         0.0625 0.0625 0.1875
         0.0625 0.4375 0.4375
         0.4375 0.6875 0.4375
         0.0625 0.0625 0.1875
         0.0625 0.4375 0.4375
         0.5625 0.8125 0.5625
         0.0625 0.0625 0.1875
         0.0625 0.4375 0.4375
         0.6875 0.9375 0.6875
         0.0625 0.0625 0.1875
         0.0625 0.5625 0.4375
         0.8125 0.5625 0.8125
         0.0625 0.0625 0.1875
         0.0625 0.5625 0.5625
         0.9375 0.6875 0.5625
         0.0625 0.0625 0.1875
         0.1875 0.5625 0.5625
         0.1875 0.8125 0.6875
         0.0625 0.0625 0.1875
         0.1875 0.6875 0.6875
         0.3125 0.6875 0.6875
         0.0625 0.0625 0.3125
         0.1875 0.6875 0.3125
         0.4375 0.8125 0.3125
         0.0625 0.0625 0.3125
         0.1875 0.8125 0.3125
         0.5625 0.8125 0.4375
         0.0625 0.1875 0.3125
         0.1875 0.1875 0.3125
         0.6875 0.1875 0.5625
         0.0625 0.1875 0.3125
         0.1875 0.1875 0.3125
         0.8125 0.3125 0.6875
         0.0625 0.1875 0.3125
         0.1875 0.1875 0.4375
         0.9375 0.4375 0.4375
         0.0625 0.1875 0.3125
         0.3125 0.1875 0.4375
         0.3125 0.5625 0.5625
         0.0625 0.1875 0.3125
         0.3125 0.1875 0.4375
         0.4375 0.6875 0.6875
         0.0625 0.1875 0.3125
         0.3125 0.1875 0.5625
         0.5625 0.8125 0.5625
         0.0625 0.1875 0.4375
         0.3125 0.3125 0.4375
         0.6875 0.3125 0.4375
         0.0625 0.1875 0.4375
         0.3125 0.3125 0.4375
         0.8125 0.4375 0.5625
        ],
        [
         -3.0797 -1.0176 5.603 5.8608 0.9944 -0.4378
         19.3076 13.0419 8.1191 7.0399 8.6946 9.5218
         20.7687 15.5184 8.8991 8.5733 15.1665 19.2621
         20.7687 18.6013 11.613 13.7268 16.3382 19.2621
         23.1341 19.5714 19.6079 18.1487 19.6303 22.2492
         23.1341 22.3329 22.4241 24.3659 24.4258 22.2492
         -2.7825 0.1391 0.426 5.3365 2.4087 0.4238
         17.6261 11.9822 10.1492 7.8454 7.8868 8.1833
         19.1894 14.4668 12.8402 8.7945 13.3831 18.0354
         20.2773 15.6004 19.029 14.631 14.5391 19.5087
         22.4674 18.5657 21.0994 16.1492 18.9094 19.7048
         24.4737 21.5049 24.0188 22.0439 24.8576 22.9994
         -2.19 1.5644 1.5649 -2.1898 4.064 1.5635
         16.0987 11.1708 9.0862 14.208 7.3417 7.1008
         17.7534 12.5881 11.7854 19.7035 10.8137 16.76
         18.9164 13.8692 16.7675 19.7035 13.9894 17.0583
         21.4298 17.7886 21.3158 23.0946 18.4051 18.8537
         21.8084 20.8345 23.8776 24.4169 24.5588 24.0153
         -1.3091 3.2368 2.9673 -1.6019 0.7109 2.9665
         14.7953 9.9115 8.2851 12.6564 8.4338 6.2852
         16.5529 10.7197 10.9631 18.2575 14.9986 13.9562
         17.6162 13.3154 14.0396 19.7001 19.2759 16.3266
         18.3713 17.2633 22.3899 21.8871 21.2448 18.147
         20.6643 20.3729 23.2849 24.1379 22.2836 25.2903
         -0.1498 5.0688 4.6119 -0.7292 1.8461 1.2792
         13.7348 7.6401 7.7495 11.3381 7.3552 6.8358
         15.2346 10.4177 10.3044 17.0422 13.9713 17.411
         15.5868 13.0521 11.5748 18.5896 17.0145 19.7706
         16.6095 17.0003 22.6878 18.9497 21.4663 20.782
         19.8033 20.1369 23.8328 24.1623 22.2909 21.3528
         1.2794 -1.3097 6.3931 0.4238 3.2462 2.4089
         12.3165 12.7881 7.4661 10.2733 6.5442 5.746
         12.985 15.3549 8.9448 15.7515 13.1255 16.438
         14.8647 21.4111 10.4212 16.0746 14.3435 17.5062
         15.7974 22.4619 22.4677 17.7757 21.5717 20.7542
         19.0943 23.4833 24.6836 23.3713 22.7623 21.8472
         2.9583 -0.438 2.6904 1.8457 4.8733 3.7928
         9.6949 11.4707 8.0227 9.4654 6.0143 4.9322
         12.428 14.1067 10.7233 12.9091 11.5323 14.689
         14.3846 19.0234 17.0676 15.3511 12.8296 15.7815
         15.2587 21.996 19.0231 16.9735 21.1473 20.215
         18.6079 23.4258 26.0923 22.7219 24.1573 23.0301
         4.8012 0.7108 4.0755 3.5136 2.9685 3.5179
         7.3729 10.4082 7.2235 8.8966 6.2744 4.6557
         12.159 13.0604 9.9082 10.3144 12.9436 15.4926
         14.1471 16.0838 14.9934 14.8699 17.3102 17.7991
         14.9896 21.1131 19.5279 16.4261 19.2582 19.6745
         18.3589 23.5045 26.1777 22.2723 24.3074 22.7873
         -2.4858 2.1277 5.686 -1.018 4.3403 2.1279
         15.9032 9.6067 6.6998 11.0932 5.4705 5.4737
         18.1436 12.1739 9.343 17.3994 12.1437 19.5275
         21.3064 13.3553 12.5062 20.7075 15.2504 19.5275
         23.3524 20.3734 20.9135 21.772 19.7671 20.5238
         23.7856 23.0358 25.7071 23.2785 24.5023 20.5238
         -1.8949 3.7901 5.4263 -0.1503 4.6029 3.2422
         14.3535 9.0573 6.4366 9.7669 5.7552 4.3789
         16.7501 10.4579 9.0996 16.202 11.3582 18.1679
         20.7072 11.9433 15.6027 19.2545 15.8453 18.5363
         21.6903 19.8656 17.3703 20.458 17.6337 20.0806
         23.3065 22.6358 26.6302 23.5595 26.7367 21.5258
        ],)

    @test parse_clock(str) == DataFrame([
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
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == false

    @test isjobdone(str) == true
end

@testset "Parse scf Si output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/example01/reference/si.scf.cg.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 2,
        alat = 10.2,
        omega = 265.302,
        nat = 2,
        ntyp = 1,
        nelec = 8.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 4,
        ecutwfc = 18.0,
        ecutrho = 72.0,
        ecutfock = nothing,
        conv_thr = 1.0e-8,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)",
        nstep = nothing,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 63 63 21
         "gvecs" "Min" 682 682 132
         "sticks" "Max" 64 64 22
         "gvecs" "Max" 686 686 135
         "sticks" "Sum" 253 253 85
         "gvecs" "Sum" 2733 2733 531
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
            0.125 0.125 0.125 0.0625
            0.125 0.125 0.375 0.1875
            0.125 0.125 0.625 0.1875
            0.125 0.125 0.875 0.1875
            0.125 0.375 0.375 0.1875
            0.125 0.375 0.625 0.375
            0.125 0.375 0.875 0.375
            0.125 0.625 0.625 0.1875
            0.375 0.375 0.375 0.0625
            0.375 0.375 0.625 0.1875
        ],
        cryst = nothing,)

    @test parse_stress(str) == ([-10.24],
        [[
          -6.961e-5 0.0 0.0
          0.0 -6.961e-5 0.0
          0.0 -0.0 -6.961e-5
        ]],
        [[
          -10.24 0.0 0.0
          0.0 -10.24 0.0
          0.0 -0.0 -10.24
        ]],)

    @test isempty(parse_cell_parameters(str))

    @test isempty(parse_atomic_positions(str))

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 18.0 CGDiagonalization() 0.01 3.0 0.7 0.1 -15.8409163 -15.86196769 0.06153621
             1 2 18.0 CGDiagonalization() 0.000769 3.0 0.7 0.1 -15.84402044 -15.84433077 0.0021625
             1 3 18.0 CGDiagonalization() 2.7e-5 3.8 0.7 0.1 -15.84450651 -15.84454243 7.508e-5
             1 4 18.0 CGDiagonalization() 9.39e-7 4.1 0.7 0.1 -15.84452598 -15.8445296 8.06e-6
             1 5 18.0 CGDiagonalization() 1.01e-7 3.9 0.7 0.1 -15.84452722 -15.84452725 7.0e-8
             1 6 18.0 CGDiagonalization() 8.22e-10 4.5 0.7 0.2 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-15.84452726,
        -15.84452726,
        1.1e-9,
        nothing,
        [4.79352741, 1.07664069, -4.8149367, -16.89975867],
        nothing,
        nothing,)

    @test parse_version(str) == "6.0"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 4)

    @test parse_fft_dimensions(str) == (2733, (nr1 = 20, nr2 = 20, nr3 = 20))

    @test parse_bands(str) == ([
         0.125 0.125 0.875
         0.125 0.875 0.125
         0.125 0.125 0.625
         0.125 0.375 0.625
         0.125 0.375 0.375
         0.375 0.125 0.375
         0.125 0.375 0.375
         0.125 0.625 0.375
         0.625 0.125 0.375
         0.125 0.375 0.625
        ],
        [
         -5.6039 3.5165 -3.549 2.1614
         4.6467 3.9919 0.3751 4.323
         5.9568 -2.4615 2.8565 -4.0849
         5.9568 -0.5936 4.2745 0.2304
         -5.0584 2.7226 -2.2719 5.1432
         3.0175 3.5069 -0.7033 5.1432
         4.9012 -4.5395 2.0784 -3.3347
         4.9909 1.5909 3.2106 -0.5842
         -3.9883 3.8905 -2.822 3.934
         1.3106 5.4637 -0.439 4.6556
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 0.01 0.02 1
         "" "electrons" 0.11 0.12 1
         "" "forces" 0.0 0.0 1
         "" "stress" 0.01 0.01 1
         "init_run" "wfcinit" 0.01 0.01 1
         "init_run" "potinit" 0.0 0.0 1
         "electrons" "c_bands" 0.09 0.09 7
         "electrons" "sum_band" 0.02 0.02 7
         "electrons" "v_of_rho" 0.0 0.0 7
         "electrons" "mix_rho" 0.0 0.0 7
         "c_bands" "init_us_2" 0.0 0.0 170
         "c_bands" "ccgdiagg" 0.08 0.07 70
         "c_bands" "wfcrot" 0.01 0.02 60
         "h_psi" "h_psi:pot" 0.06 0.06 828
         "h_psi" "h_psi:calbec" 0.0 0.01 828
         "h_psi" "vloc_psi" 0.05 0.05 828
         "h_psi" "add_vuspsi" 0.0 0.0 828
         "h_psi" "h_1psi" 0.05 0.05 768
         "General routines" "calbec" 0.01 0.01 1646
         "General routines" "fft" 0.0 0.0 34
         "General routines" "fftw" 0.05 0.05 2376
         "General routines" "davcio" 0.0 0.0 10
         "Parallel routines" "fft_scatter" 0.02 0.02 2410
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == false

    @test isjobdone(str) == true
end

@testset "Parse scf NaCl output" begin
    url = "https://raw.githubusercontent.com/atztogo/phonopy/master/example/NaCl-pwscf/NaCl-001.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 64,
        alat = 21.5062,
        omega = 9947.007,
        nat = 64,
        ntyp = 2,
        nelec = 512.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 256,
        ecutwfc = 70.0,
        ecutrho = 280.0,
        ecutfock = nothing,
        conv_thr = 1.0e-9,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PW   PBX  PBC ( 1  4  3  4 0 0)",
        nstep = nothing,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 644 644 165
         "gvecs" "Min" 49140 49140 6426
         "sticks" "Max" 645 645 167
         "gvecs" "Max" 49142 49142 6427
         "sticks" "Sum" 10309 10309 2661
         "gvecs" "Sum" 786247 786247 102831
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [0.25 0.25 0.25 2.0], cryst = nothing)

    @test parse_stress(str) == ([0.63],
        [[
          4.28e-6 0.0 0.0
          0.0 4.27e-6 0.0
          0.0 0.0 4.27e-6
        ]],
        [[
          0.63 0.0 0.0
          0.0 0.63 0.0
          0.0 0.0 0.63
        ]],)

    @test isempty(parse_cell_parameters(str))

    @test isempty(parse_atomic_positions(str))

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 70.0 DavidsonDiagonalization() 0.01 1.0 0.7 10.7 -6041.20603805 -6045.69337535 5.2864571
             1 2 70.0 DavidsonDiagonalization() 0.00103 4.0 0.7 22.7 -6043.01689773 -6043.59224951 0.77895678
             1 3 70.0 DavidsonDiagonalization() 0.000152 2.0 0.7 32.3 -6043.23571556 -6043.25444434 0.0214773
             1 4 70.0 DavidsonDiagonalization() 4.19e-6 4.0 0.7 44.8 -6043.24515439 -6043.24634462 0.00257225
             1 5 70.0 DavidsonDiagonalization() 5.02e-7 2.0 0.7 52.9 -6043.24549219 -6043.24548845 2.587e-5
             1 6 70.0 DavidsonDiagonalization() 5.05e-9 4.0 0.7 68.3 -6043.24556318 -6043.24556762 9.92e-6
             1 7 70.0 DavidsonDiagonalization() 1.94e-9 2.0 0.7 75.9 -6043.24556342 -6043.24556419 1.4e-6
             1 8 70.0 DavidsonDiagonalization() 2.73e-10 3.0 0.7 85.6 -6043.24556395 -6043.24556397 7.0e-8
             1 9 70.0 DavidsonDiagonalization() 1.29e-11 2.0 0.7 94.5 -6043.24556394 -6043.24556396 2.0e-8
             1 10 70.0 DavidsonDiagonalization() 3.22e-12 3.0 0.7 103.6 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-6043.24556394,
        -6043.24556394,
        9.5e-10,
        -39940.429322,
        [-2623.52269882, 1363.88574749, -789.46164384, -2182.32845246],
        nothing,
        nothing,)

    @test parse_version(str) == "5.4.0"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 16)

    @test parse_fft_dimensions(str) == (786247, (nr1 = 120, nr2 = 120, nr3 = 120))

    @test parse_bands(str) == ([0.25 0.25 0.25],
        [-48.1202 -48.1201 -48.1201 -48.12 -48.1187 -48.1185 -48.1185 -48.1184 -48.1183 -48.1181 -48.1181 -48.1181 -48.1181 -48.1181 -48.1181 -48.118 -48.118 -48.118 -48.118 -48.1179 -48.1179 -48.1179 -48.1178 -48.1178 -48.1176 -48.1175 -48.1175 -48.1174 -48.1161 -48.116 -48.116 -48.1159 -20.0531 -20.0529 -20.0529 -20.0515 -20.0488 -20.0486 -20.0484 -20.0483 -20.0483 -20.0482 -20.0481 -20.048 -20.0478 -20.047 -20.0466 -20.0465 -20.0461 -20.0459 -20.0455 -20.0455 -20.0454 -20.0452 -20.0451 -20.0449 -20.0448 -20.0447 -20.0447 -20.0444 -20.0403 -20.0402 -20.0402 -20.0401 -20.0401 -20.04 -20.0395 -20.0394 -20.0393 -20.039 -20.039 -20.0388 -20.0381 -20.0379 -20.0378 -20.0367 -20.0365 -20.0365 -20.0364 -20.0363 -20.0361 -20.0361 -20.036 -20.0359 -20.0357 -20.0357 -20.0356 -20.0355 -20.0353 -20.035 -20.0346 -20.0345 -20.0344 -20.0341 -20.034 -20.0339 -20.0337 -20.0337 -20.0334 -20.0334 -20.0331 -20.0331 -20.033 -20.0327 -20.0326 -20.0326 -20.0326 -20.0325 -20.0325 -20.0324 -20.0322 -20.032 -20.0319 -20.0319 -20.0316 -20.0314 -20.0314 -20.0312 -20.0311 -20.031 -20.0308 -20.0305 -20.0304 -20.0302 -20.0301 -20.0296 -20.0294 -20.0293 -11.714 -11.6188 -11.6179 -11.6172 -11.5458 -11.545 -11.5442 -11.4969 -11.4925 -11.4917 -11.4909 -11.456 -11.4543 -11.4542 -11.4542 -11.4542 -11.4529 -11.4431 -11.443 -11.4428 -11.4169 -11.416 -11.4154 -11.4102 -11.4101 -11.4099 -11.4008 -11.3996 -11.3996 -11.3996 -11.3996 -11.3976 -0.8874 -0.8874 -0.8874 -0.783 -0.6658 -0.6658 -0.6658 -0.6358 -0.6349 -0.6341 -0.6261 -0.6257 -0.6254 -0.6253 -0.6249 -0.6246 -0.5201 -0.5201 -0.5195 -0.5195 -0.519 -0.5189 -0.483 -0.4824 -0.4818 -0.4217 -0.4216 -0.4215 -0.1146 -0.1146 -0.1143 -0.1143 -0.114 -0.114 -0.0591 -0.0588 -0.0585 -0.0383 -0.0381 -0.0378 0.1489 0.1489 0.1489 0.1899 0.1901 0.1905 0.1905 0.1909 0.1911 0.2862 0.2864 0.2869 0.2869 0.2874 0.2876 0.3165 0.3165 0.3165 0.3739 0.3739 0.3739 0.4176 0.4176 0.4176 0.449 0.4494 0.4499 0.469 0.4696 0.47 0.4839 0.4843 0.4847 0.4898 0.502 0.5028 0.503 0.503 0.5031 0.5041 0.5436 0.5436 0.5436 0.6285 0.6285 0.6504 0.6509 0.6513 0.6599 0.6608 0.6616 0.6888 0.6892 0.6896 0.8041 0.8041],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 4.14 4.3 1
         "" "electrons" 98.15 99.0 1
         "" "forces" 1.6 1.65 1
         "" "stress" 11.58 11.58 1
         "init_run" "wfcinit" 3.2 3.26 1
         "init_run" "potinit" 0.31 0.32 1
         "electrons" "c_bands" 84.89 85.71 10
         "electrons" "sum_band" 9.46 9.49 10
         "electrons" "v_of_rho" 1.01 1.02 11
         "electrons" "newd" 1.14 1.19 11
         "electrons" "PAW_pot" 1.73 1.73 11
         "electrons" "mix_rho" 0.2 0.2 10
         "c_bands" "init_us_2" 0.4 0.41 21
         "c_bands" "cegterg" 83.49 84.29 10
         "sum_band" "sum_band:bec" 0.01 0.01 10
         "sum_band" "addusdens" 1.5 1.51 10
         "*egterg" "h_psi" 41.33 41.44 38
         "*egterg" "s_psi" 5.69 5.69 38
         "*egterg" "g_psi" 0.24 0.24 27
         "*egterg" "cdiaghg" 15.26 15.27 37
         "h_psi" "add_vuspsi" 5.67 5.68 38
         "General routines" "calbec" 9.12 9.13 53
         "General routines" "fft" 0.73 0.77 167
         "General routines" "fftw" 30.1 30.19 17922
         "Parallel routines" "fft_scatter" 10.96 10.94 18089
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "NaCl-001.in"

    @test isrelaxed(str) == false

    @test isjobdone(str) == true
end

@testset "Parse scf SiO2 output" begin
    url = "https://raw.githubusercontent.com/maxhutch/deprecated-quantum-espresso/master/XSpectra/examples/reference/SiO2.scf.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 9,
        alat = 9.2863,
        omega = 762.9417,
        nat = 9,
        ntyp = 3,
        nelec = 48.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 30,
        ecutwfc = 20.0,
        ecutrho = 150.0,
        ecutfock = nothing,
        conv_thr = 1.0e-9,
        mixing_beta = 0.3,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PW   PBE  PBE ( 1  4  3  4 0 0)",
        nstep = nothing,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Sum" 889 475 151
         "gvecs" "Sum" 23595 9203 1559
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
            0.0 0.0 0.0 0.25
            0.0 0.0 -0.4545041 0.25
            0.0 -0.5773503 0.0 0.25
            0.0 -0.5773503 -0.4545041 0.25
            0.5 -0.2886751 0.0 0.5
            0.5 -0.2886751 -0.4545041 0.5
        ],
        cryst = [
            0.0 0.0 0.0 0.25
            0.0 0.0 -0.5 0.25
            0.0 -0.5 0.0 0.25
            0.0 -0.5 -0.5 0.25
            0.5 -0.5 0.0 0.5
            0.5 -0.5 -0.5 0.5
        ],)

    @test all(isempty, parse_stress(str))

    @test isempty(parse_cell_parameters(str))

    @test isempty(parse_atomic_positions(str))

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 20.0 DavidsonDiagonalization() 0.01 3.3 0.3 4.8 -215.47818098 -215.53465529 0.3369274
             1 2 20.0 DavidsonDiagonalization() 0.000702 1.0 0.3 5.8 -215.48170145 -215.48997319 0.1173383
             1 3 20.0 DavidsonDiagonalization() 0.000244 2.7 0.3 6.9 -215.48819658 -215.48809974 0.01949303
             1 4 20.0 DavidsonDiagonalization() 4.06e-5 4.8 0.3 8.1 -215.4897598 -215.48959757 0.0004444
             1 5 20.0 DavidsonDiagonalization() 9.26e-7 5.3 0.3 9.8 -215.49021498 -215.49010809 3.11e-5
             1 6 20.0 DavidsonDiagonalization() 6.48e-8 3.0 0.3 11.1 -215.49027155 -215.49022257 9.02e-6
             1 7 20.0 DavidsonDiagonalization() 1.88e-8 3.0 0.3 12.5 -215.49029996 -215.49027536 3.07e-6
             1 8 20.0 DavidsonDiagonalization() 6.39e-9 2.2 0.3 13.9 -215.49031266 -215.49030052 4.0e-8
             1 9 20.0 DavidsonDiagonalization() 7.44e-11 3.0 0.3 15.3 -215.49031863 -215.49031268 3.5e-9
             1 10 20.0 DavidsonDiagonalization() 7.39e-12 3.0 0.3 16.8 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-215.49032431,
        -215.49031863,
        5.2e-10,
        nothing,
        [-102.06258058, 76.26219017, -50.59625207, -139.09368184],
        nothing,
        0.0,)

    @test parse_version(str) == "5.2.0"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 1)

    @test parse_fft_dimensions(str) == (23595, (nr1 = 40, nr2 = 40, nr3 = 40))

    @test parse_bands(str) == ([
         0.0 0.0 0.5
         0.0 -0.5774 -0.2887
         0.0 0.0 0.0
         0.0 0.0 0.5
         0.0 -0.5774 -0.2887
         -0.4545 -0.4545 -0.4545
        ],
        [
         -17.0334 -5.7053 -0.3559 1.0458 8.5471 -16.5782 -6.7838 -0.1516 0.6782 10.5342 -16.5545 -5.8544 -0.4926 0.8894 10.8541 -16.4867 -5.6872 -0.617 1.2993 11.5725 -16.554 -5.8545 -0.4923 0.8903 10.8543 -16.4862 -5.6873 -0.617 1.2994 11.5727
         -15.5439 -5.7051 -0.3556 1.561 12.0141 -16.5775 -3.7867 -0.151 0.9906 10.535 -16.0878 -4.8542 -0.3446 1.4566 11.6301 -16.0555 -4.7586 -0.591 1.4355 12.0567 -16.0884 -4.8539 -0.3447 1.4561 11.6302 -16.056 -4.7583 -0.5911 1.4353 12.057
         -15.5432 -2.9612 -0.1432 1.5616 12.0141 -14.9317 -3.7863 -0.0524 1.3044 11.7234 -15.4934 -3.3991 -0.1052 1.5187 12.3931 -15.6172 -3.9004 0.1213 1.6749 12.3288 -15.4938 -3.3993 -0.1058 1.5191 12.3931 -15.6177 -3.9004 0.1214 1.6744 12.3284
         -14.9187 -2.9296 0.1125 1.8715 12.4598 -14.9147 -3.113 0.5108 2.3399 12.4599 -14.9799 -2.9537 0.2737 1.7938 12.4139 -15.0016 -2.8272 0.4585 1.7248 12.5029 -14.9796 -2.9541 0.273 1.7939 12.4142 -15.002 -2.8276 0.4579 1.7248 12.5031
         -14.8904 -2.742 1.0266 2.2461 12.46 -14.9139 -3.1126 0.5114 2.34 12.4602 -14.9429 -2.916 0.4648 2.3038 13.2561 -14.9318 -2.5641 0.4729 1.9562 12.9247 -14.9439 -2.9156 0.4649 2.3042 13.2561 -14.9318 -2.564 0.4729 1.9573 12.9248
         -14.8896 -2.7419 1.027 2.2465 13.768 -14.8691 -2.2079 0.678 2.4292 13.3749 -14.9202 -2.1394 0.8078 2.5349 13.3438 -14.9111 -2.4471 1.0148 2.3604 13.1935 -14.919 -2.1395 0.8082 2.5344 13.344 -14.91 -2.447 1.015 2.3604 13.1935
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 1.89 2.01 1
         "" "electrons" 14.18 14.52 1
         "init_run" "wfcinit" 0.69 0.69 1
         "init_run" "potinit" 0.12 0.16 1
         "electrons" "c_bands" 11.08 11.18 11
         "electrons" "sum_band" 2.15 2.22 11
         "electrons" "v_of_rho" 0.5 0.52 11
         "electrons" "v_h" 0.03 0.03 11
         "electrons" "v_xc" 0.47 0.5 11
         "electrons" "newd" 0.45 0.56 11
         "electrons" "mix_rho" 0.05 0.08 11
         "c_bands" "init_us_2" 0.12 0.13 138
         "c_bands" "cegterg" 10.79 10.88 66
         "sum_band" "sum_band:bec" 0.0 0.0 66
         "sum_band" "addusdens" 0.58 0.65 11
         "*egterg" "h_psi" 6.87 6.89 270
         "*egterg" "s_psi" 0.91 0.91 270
         "*egterg" "g_psi" 0.07 0.07 198
         "*egterg" "cdiaghg" 0.85 0.86 258
         "*egterg" "cegterg:over" 1.14 1.14 198
         "*egterg" "cegterg:upda" 0.63 0.63 198
         "*egterg" "cegterg:last" 0.35 0.35 66
         "h_psi" "h_psi:vloc" 4.93 4.95 270
         "h_psi" "h_psi:vnl" 1.92 1.92 270
         "h_psi" "add_vuspsi" 0.91 0.91 270
         "General routines" "calbec" 1.36 1.36 336
         "General routines" "fft" 0.26 0.29 175
         "General routines" "ffts" 0.01 0.01 22
         "General routines" "fftw" 4.87 4.89 13186
         "General routines" "interpolate" 0.05 0.05 22
         "Parallel routines" "fft_scatter" 0.33 0.34 13383
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == false

    @test isjobdone(str) == true
end

@testset "Parse vc-relax As output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/VCSexample/reference/As.vcs00.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 2,
        alat = 7.0103,
        omega = 245.3705,
        nat = 2,
        ntyp = 1,
        nelec = 10.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 9,
        ecutwfc = 25.0,
        ecutrho = 100.0,
        ecutfock = nothing,
        conv_thr = 1.0e-7,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)",
        nstep = 55,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 174 174 60
         "gvecs" "Min" 2079 2079 416
         "sticks" "Max" 175 175 61
         "gvecs" "Max" 2080 2080 417
         "sticks" "Sum" 349 349 121
         "gvecs" "Sum" 4159 4159 833
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
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
        cryst = nothing,)

    @test parse_stress(str) == ([
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
        ],)

    @test parse_cell_parameters(str) == [
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

     # @test parse_atomic_positions(str) ==
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

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 25.0 DavidsonDiagonalization() 0.01 4.3 0.7 0.3 -25.43995582 -25.4437112 0.01555792
             1 2 25.0 DavidsonDiagonalization() 0.000156 1.0 0.7 0.4 -25.44012469 -25.44030751 0.00088879
             1 3 25.0 DavidsonDiagonalization() 8.89e-6 1.6 0.7 0.4 -25.44015903 -25.44016035 5.01e-6
             1 4 25.0 DavidsonDiagonalization() 5.01e-8 3.1 0.7 0.5 -25.44016731 -25.44016767 7.7e-7
             1 5 25.0 DavidsonDiagonalization() 7.74e-9 1.4 0.7 0.5 nothing nothing nothing
             2 1 25.0 DavidsonDiagonalization() 1.0e-6 4.2 0.7 0.8 -25.45861709 -25.70454741 0.00082435
             2 2 25.0 DavidsonDiagonalization() 8.24e-6 3.1 0.7 0.9 -25.46013239 -25.46040686 0.00067849
             2 3 25.0 DavidsonDiagonalization() 6.78e-6 1.0 0.7 0.9 -25.46011089 -25.46016215 0.00014927
             2 4 25.0 DavidsonDiagonalization() 1.49e-6 1.0 0.7 1.0 -25.46009315 -25.46011707 4.631e-5
             2 5 25.0 DavidsonDiagonalization() 4.63e-7 2.6 0.7 1.0 -25.46010105 -25.46010159 1.04e-6
             2 6 25.0 DavidsonDiagonalization() 1.04e-8 2.1 0.7 1.1 -25.46010129 -25.46010136 1.9e-7
             2 7 25.0 DavidsonDiagonalization() 1.89e-9 1.1 0.7 1.1 nothing nothing nothing
             3 1 25.0 DavidsonDiagonalization() 1.0e-6 4.9 0.7 1.4 -25.47745488 -25.91214701 0.00268934
             3 2 25.0 DavidsonDiagonalization() 2.69e-5 3.1 0.7 1.5 -25.48276145 -25.48371495 0.00243353
             3 3 25.0 DavidsonDiagonalization() 2.43e-5 1.0 0.7 1.5 -25.4826746 -25.48286068 0.0005678
             3 4 25.0 DavidsonDiagonalization() 5.68e-6 1.0 0.7 1.6 -25.48260137 -25.48269578 0.00018831
             3 5 25.0 DavidsonDiagonalization() 1.88e-6 2.2 0.7 1.6 -25.48262653 -25.48262997 6.49e-6
             3 6 25.0 DavidsonDiagonalization() 6.49e-8 2.5 0.7 1.7 -25.48262992 -25.48263004 4.3e-7
             3 7 25.0 DavidsonDiagonalization() 4.28e-9 1.8 0.7 1.7 nothing nothing nothing
             4 1 25.0 DavidsonDiagonalization() 1.0e-6 5.1 0.7 2.0 -25.49248333 -25.62939439 0.00037033
             4 2 25.0 DavidsonDiagonalization() 3.7e-6 3.0 0.7 2.1 -25.49296778 -25.49305279 0.00020532
             4 3 25.0 DavidsonDiagonalization() 2.05e-6 1.0 0.7 2.1 -25.4929614 -25.49297819 3.596e-5
             4 4 25.0 DavidsonDiagonalization() 3.6e-7 1.5 0.7 2.2 -25.49296383 -25.49296488 1.95e-6
             4 5 25.0 DavidsonDiagonalization() 1.95e-8 3.0 0.7 2.3 -25.49296527 -25.4929653 2.0e-7
             4 6 25.0 DavidsonDiagonalization() 2.0e-9 1.0 0.7 2.3 -25.49296521 -25.49296527 1.2e-7
             4 7 25.0 DavidsonDiagonalization() 1.16e-9 1.5 0.7 2.3 nothing nothing nothing
             5 1 25.0 DavidsonDiagonalization() 1.0e-6 5.3 0.7 2.6 -25.49636723 -25.64573313 0.00052245
             5 2 25.0 DavidsonDiagonalization() 5.22e-6 3.0 0.7 2.7 -25.49696036 -25.49706357 0.00024087
             5 3 25.0 DavidsonDiagonalization() 2.41e-6 1.0 0.7 2.7 -25.49695895 -25.49697471 3.557e-5
             5 4 25.0 DavidsonDiagonalization() 3.56e-7 1.3 0.7 2.8 -25.49696079 -25.4969622 2.63e-6
             5 5 25.0 DavidsonDiagonalization() 2.63e-8 3.0 0.7 2.9 -25.49696251 -25.49696257 2.8e-7
             5 6 25.0 DavidsonDiagonalization() 2.83e-9 1.0 0.7 2.9 -25.49696244 -25.49696252 1.4e-7
             5 7 25.0 DavidsonDiagonalization() 1.44e-9 2.0 0.7 2.9 nothing nothing nothing
             6 1 25.0 DavidsonDiagonalization() 1.0e-6 4.2 0.7 3.2 -25.49530736 -25.61150725 0.00014094
             6 2 25.0 DavidsonDiagonalization() 1.41e-6 3.0 0.7 3.3 -25.49570743 -25.49578663 0.00023062
             6 3 25.0 DavidsonDiagonalization() 1.41e-6 1.0 0.7 3.4 -25.49568509 -25.49571515 6.854e-5
             6 4 25.0 DavidsonDiagonalization() 6.85e-7 1.0 0.7 3.4 -25.49568478 -25.49568972 8.37e-6
             6 5 25.0 DavidsonDiagonalization() 8.37e-8 3.0 0.7 3.5 -25.49568754 -25.49568767 4.3e-7
             6 6 25.0 DavidsonDiagonalization() 4.29e-9 1.0 0.7 3.5 -25.49568746 -25.49568755 1.5e-7
             6 7 25.0 DavidsonDiagonalization() 1.47e-9 2.0 0.7 3.6 nothing nothing nothing
             7 1 25.0 DavidsonDiagonalization() 1.0e-6 4.2 0.7 3.8 -25.49684745 -25.30507535 0.00061972
             7 2 25.0 DavidsonDiagonalization() 6.2e-6 3.0 0.7 3.9 -25.49795324 -25.49818001 0.0005669
             7 3 25.0 DavidsonDiagonalization() 5.67e-6 1.0 0.7 4.0 -25.49793986 -25.49798289 0.0001145
             7 4 25.0 DavidsonDiagonalization() 1.15e-6 1.0 0.7 4.0 -25.49793038 -25.49794575 2.62e-5
             7 5 25.0 DavidsonDiagonalization() 2.62e-7 3.0 0.7 4.1 -25.49793823 -25.49793843 9.5e-7
             7 6 25.0 DavidsonDiagonalization() 9.51e-9 1.0 0.7 4.1 -25.49793797 -25.49793825 5.0e-7
             7 7 25.0 DavidsonDiagonalization() 4.98e-9 2.0 0.7 4.2 nothing nothing nothing
             8 1 25.0 DavidsonDiagonalization() 1.0e-6 3.7 0.7 4.4 -25.49828297 -25.40005836 0.00019031
             8 2 25.0 DavidsonDiagonalization() 1.9e-6 3.0 0.7 4.5 -25.49857725 -25.49864174 0.00015415
             8 3 25.0 DavidsonDiagonalization() 1.54e-6 1.0 0.7 4.5 -25.49857749 -25.49858675 2.667e-5
             8 4 25.0 DavidsonDiagonalization() 2.67e-7 1.0 0.7 4.6 -25.49857458 -25.49857875 7.53e-6
             8 5 25.0 DavidsonDiagonalization() 7.53e-8 2.7 0.7 4.7 nothing nothing nothing
             9 1 25.0 DavidsonDiagonalization() 1.0e-6 3.8 0.7 4.9 -25.49832741 -25.38301822 0.00024105
             9 2 25.0 DavidsonDiagonalization() 2.41e-6 3.0 0.7 5.0 -25.49870513 -25.49879485 0.00021533
             9 3 25.0 DavidsonDiagonalization() 2.15e-6 1.0 0.7 5.0 -25.49870766 -25.49871944 3.584e-5
             9 4 25.0 DavidsonDiagonalization() 3.58e-7 1.0 0.7 5.1 -25.49870296 -25.4987092 1.123e-5
             9 5 25.0 DavidsonDiagonalization() 1.12e-7 2.7 0.7 5.1 nothing nothing nothing
             10 1 25.0 DavidsonDiagonalization() 1.0e-6 3.6 0.7 5.4 -25.49901709 -25.53840213 3.743e-5
             10 2 25.0 DavidsonDiagonalization() 3.74e-7 3.0 0.7 5.5 -25.49906083 -25.49906974 1.984e-5
             10 3 25.0 DavidsonDiagonalization() 1.98e-7 1.0 0.7 5.5 -25.49906136 -25.49906237 2.5e-6
             10 4 25.0 DavidsonDiagonalization() 2.5e-8 1.0 0.7 5.6 -25.49906132 -25.49906154 4.1e-7
             10 5 25.0 DavidsonDiagonalization() 4.09e-9 3.0 0.7 5.6 nothing nothing nothing
             11 1 25.0 DavidsonDiagonalization() 1.0e-6 3.4 0.7 5.9 -25.49926138 -25.47952947 1.028e-5
             11 2 25.0 DavidsonDiagonalization() 1.03e-7 2.7 0.7 6.0 -25.49927224 -25.49927575 8.21e-6
             11 3 25.0 DavidsonDiagonalization() 8.21e-8 1.0 0.7 6.0 -25.4992721 -25.49927285 1.21e-6
             11 4 25.0 DavidsonDiagonalization() 1.21e-8 2.2 0.7 6.1 -25.49927238 -25.49927249 2.1e-7
             11 5 25.0 DavidsonDiagonalization() 2.13e-9 1.7 0.7 6.1 nothing nothing nothing
             12 1 25.0 DavidsonDiagonalization() 1.0e-6 3.5 0.7 6.4 -25.49938883 -25.47370616 1.302e-5
             12 2 25.0 DavidsonDiagonalization() 1.3e-7 3.0 0.7 6.5 -25.49940541 -25.49941055 1.25e-5
             12 3 25.0 DavidsonDiagonalization() 1.25e-7 1.0 0.7 6.6 -25.49940503 -25.49940626 2.09e-6
             12 4 25.0 DavidsonDiagonalization() 2.09e-8 2.0 0.7 6.6 nothing nothing nothing
             13 1 25.0 DavidsonDiagonalization() 1.0e-6 3.6 0.7 6.9 -25.4994147 -25.46715169 1.51e-5
             13 2 25.0 DavidsonDiagonalization() 1.51e-7 3.0 0.7 7.0 -25.49943827 -25.49944506 1.723e-5
             13 3 25.0 DavidsonDiagonalization() 1.51e-7 1.0 0.7 7.0 -25.49943774 -25.49943933 3.25e-6
             13 4 25.0 DavidsonDiagonalization() 3.25e-8 1.4 0.7 7.1 -25.499438 -25.49943808 1.4e-7
             13 5 25.0 DavidsonDiagonalization() 1.36e-9 3.0 0.7 7.1 nothing nothing nothing
             14 1 25.0 DavidsonDiagonalization() 1.0e-6 3.1 0.7 7.4 -25.49942959 -25.46106847 1.52e-5
             14 2 25.0 DavidsonDiagonalization() 1.52e-7 3.0 0.7 7.5 -25.49946014 -25.49946654 1.647e-5
             14 3 25.0 DavidsonDiagonalization() 1.52e-7 1.0 0.7 7.5 -25.49945989 -25.49946105 3.4e-6
             14 4 25.0 DavidsonDiagonalization() 3.4e-8 1.0 0.7 7.6 -25.49945947 -25.49946003 9.5e-7
             14 5 25.0 DavidsonDiagonalization() 9.49e-9 2.7 0.7 7.6 nothing nothing nothing
             15 1 25.0 DavidsonDiagonalization() 1.0e-6 2.2 0.7 8.0 -25.49946121 -25.5097648 1.27e-6
             15 2 25.0 DavidsonDiagonalization() 1.27e-8 3.0 0.7 8.1 -25.49946359 -25.49946413 1.37e-6
             15 3 25.0 DavidsonDiagonalization() 1.27e-8 1.0 0.7 8.1 -25.49946358 -25.49946367 2.7e-7
             15 4 25.0 DavidsonDiagonalization() 2.69e-9 1.0 0.7 8.1 nothing nothing nothing
             16 1 25.0 DavidsonDiagonalization() 1.0e-6 1.2 0.7 8.4 -25.49946518 -25.49949373 1.1e-7
             16 2 25.0 DavidsonDiagonalization() 1.06e-9 2.0 0.7 8.5 nothing nothing nothing
             17 1 25.0 DavidsonDiagonalization() 1.0e-6 1.6 0.7 8.8 nothing nothing nothing
             18 1 25.0 DavidsonDiagonalization() 1.0e-6 1.9 0.7 9.1 nothing nothing nothing
             19 1 25.0 DavidsonDiagonalization() 1.0e-6 1.0 0.7 9.5 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-25.4401674,
        -25.44016741,
        2.0e-8,
        nothing,
        nothing,
        nothing,
        nothing,)

    @test parse_version(str) == "6.0"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 2)

    @test parse_fft_dimensions(str) == (4159, (nr1 = 24, nr2 = 24, nr3 = 24))

    @test parse_bands(str) == ([
         0.0 0.2398 0.1432
         0.0 0.0469 0.5611
         0.1535 -0.277 -0.0
         -0.1436 0.0 -0.2386
         -0.2488 0.328 0.0
         0.2558 0.1385 0.0
         0.2873 0.7195 0.4295
         0.4976 0.0469 0.4208
         -0.0512 0.0 0.7289
         0.1436 0.4797 0.1432
         0.2488 0.1406 0.0
         0.0512 0.5539 -0.0
         -0.2873 0.0 0.1439
         0.0 -0.2343 -0.1404
         0.3581 0.0 -0.2432
         0.1436 0.0 0.2398
         0.7464 0.4217 0.2808
         0.0512 0.4154 0.4863
         0.0 0.7195 -0.048
         0.4976 0.1406 0.1404
         0.1535 0.0 0.2432
         0.5746 -0.0 0.048
         0.0 0.1408 -0.2808
         -0.2558 -0.1396 0.0
         0.0 -0.2418 0.3357
         0.0 0.2346 0.1404
         0.4604 0.2792 0.7295
         0.4309 0.4836 0.048
         0.7464 -0.0469 0.0
         0.1535 0.1396 0.4863
         0.0 0.2418 0.1439
         0.0 0.0469 0.5616
         0.152 -0.2792 -0.0
         -0.1413 0.0 -0.2398
         -0.2448 0.3285 0.0
         0.2534 0.1396 0.0
         0.2826 0.7253 0.4316
         0.4895 0.0469 0.4212
         -0.0507 0.0 0.7295
         0.1413 0.4836 0.1439
         0.2448 0.1408 0.0
         0.0507 0.5584 0.0
         -0.2826 0.0 0.1436
         0.0 -0.2346 -0.1404
         0.3547 0.0 -0.2431
         0.1413 0.0 0.2394
         0.7343 0.4223 0.2807
         0.0507 0.4188 0.4863
         0.0 0.7253 -0.0479
         0.4895 0.1408 0.1404
         0.152 0.0 0.2431
         0.5652 -0.0 0.0479
         0.0 0.1411 -0.2807
         -0.2534 -0.1408 0.0
         0.0 -0.2439 0.3352
         0.0 0.2352 0.1404
         0.4561 0.2817 0.7294
         0.4239 0.4879 0.0479
         0.7343 -0.047 0.0
         0.152 0.1408 0.4863
         0.0 0.2439 0.1436
         0.0 0.047 0.5615
         0.149 -0.2817 -0.0
         -0.1372 0.0 -0.2394
         -0.2376 0.3293 0.0
         0.2484 0.1408 0.0
         0.2744 0.7318 0.4309
         0.4753 0.047 0.4211
         -0.0497 0.0 0.7294
         0.1372 0.4879 0.1436
         0.2376 0.1411 0.0
         0.0497 0.5633 -0.0
         -0.2744 0.0 0.1436
         0.0 -0.2352 -0.1404
         0.3477 0.0 -0.2431
         0.1372 0.0 0.2394
         0.7129 0.4233 0.2808
         0.0497 0.4225 0.4863
         0.0 0.7318 -0.0479
         0.4753 0.1411 0.1404
         0.149 0.0 0.2431
         0.5488 -0.0 0.0479
         0.0 0.1415 -0.2808
         -0.2484 -0.1402 -0.0
         0.0 -0.2428 0.3352
         0.0 0.2358 0.1404
         0.447 0.2803 0.7294
         0.4116 0.4855 0.0479
         0.7129 -0.0472 0.0
         0.149 0.1402 0.4863
         0.0 0.2428 0.1436
         0.0 0.0472 0.5615
         0.1453 -0.2803 0.0
         -0.1372 0.0 -0.2394
         -0.2376 0.3301 0.0
         0.2422 0.1402 0.0
         0.2744 0.7283 0.4309
         0.4753 0.0472 0.4211
         -0.0484 0.0 0.7294
         0.1372 0.4855 0.1436
         0.2376 0.1415 0.0
         0.0484 0.5606 -0.0
         -0.2744 -0.0 0.1436
         0.0 -0.2358 -0.1404
         0.3391 0.0 -0.2431
         0.1372 0.0 0.2394
         0.7129 0.4245 0.2808
         0.0484 0.4205 0.4863
         0.0 0.7283 -0.0479
         0.4753 0.1415 0.1404
         0.1453 0.0 0.2431
         0.5488 0.0 0.0479
         0.0 0.142 -0.2808
         -0.2422 -0.1402 -0.0
         0.0 -0.2428 0.3351
         0.0 0.2366 0.1404
         0.436 0.2803 0.7294
         0.4116 0.4856 0.0479
         0.7129 -0.0473 0.0
         0.1453 0.1402 0.4863
         0.0 0.2428 0.1436
         0.0 0.0473 0.5615
         0.1413 -0.2803 0.0
         -0.1373 0.0 -0.2394
         -0.2377 0.3312 0.0
         0.2355 0.1402 0.0
         0.2745 0.7283 0.4309
         0.4755 0.0473 0.4211
         -0.0471 0.0 0.7294
         0.1373 0.4856 0.1436
         0.2377 0.142 0.0
         0.0471 0.5607 -0.0
         -0.2745 -0.0 0.1436
         0.0 -0.2366 -0.1404
         0.3297 0.0 -0.2431
         0.1373 0.0 0.2394
         0.7132 0.4259 0.2808
         0.0471 0.4205 0.4863
         0.0 0.7283 -0.0479
         0.4755 0.142 0.1404
         0.1413 0.0 0.2431
         0.549 0.0 0.0479
         0.0 0.1425 -0.2808
         -0.2355 -0.1402 0.0
         0.0 -0.2428 0.3351
         0.0 0.2375 0.1404
         0.4239 0.2804 0.7294
         0.4118 0.4857 0.0479
         0.7132 -0.0475 0.0
         0.1413 0.1402 0.4863
         0.0 0.2428 0.1436
         0.0 0.0475 0.5615
         0.1374 -0.2804 0.0
         -0.1377 0.0 -0.2394
         -0.2385 0.3325 0.0
         0.229 0.1402 0.0
         0.2754 0.7285 0.4308
         0.4769 0.0475 0.4211
         -0.0458 0.0 0.7294
         0.1377 0.4857 0.1436
         0.2385 0.1425 0.0
         0.0458 0.5608 -0.0
         -0.2754 -0.0 0.1436
         0.0 -0.2375 -0.1404
         0.3205 0.0 -0.2431
         0.1377 0.0 0.2393
         0.7154 0.4275 0.2807
         0.0458 0.4206 0.4863
         0.0 0.7285 -0.0479
         0.4769 0.1425 0.1404
         0.1374 0.0 0.2431
         0.5507 0.0 0.0479
         0.0 0.1432 -0.2807
         -0.229 -0.1403 -0.0
         0.0 -0.243 0.335
         0.0 0.2386 0.1404
         0.4121 0.2806 0.7294
         0.413 0.4859 0.0479
         0.7154 -0.0477 0.0
         0.1374 0.1403 0.4863
         0.0 0.243 0.1436
         -0.0 0.0477 0.5615
         0.1406 -0.2806 0.0
         -0.1385 0.0 -0.2393
         -0.2398 0.334 0.0
         0.2343 0.1403 -0.0
         0.277 0.7289 0.4308
         0.4797 0.0477 0.4211
         -0.0469 0.0 0.7294
         0.1385 0.4859 0.1436
        ],
        [
         -6.996 -0.9951 3.9837 4.7763 5.4039 9.7649 11.0995 14.3053 15.7304
         4.5198 4.2049 4.0324 6.0774 9.9174 12.2869 12.1495 14.9467 -4.899
         5.9668 4.6912 5.5262 8.8711 11.8325 13.4241 13.6575 -5.9081 -2.0664
         5.9668 7.4545 8.3361 9.3404 11.9263 13.8207 -4.4938 -1.5546 2.1328
         8.4358 8.4354 9.0633 10.8065 13.8731 -5.1281 -1.9411 5.8115 4.6465
         11.0403 9.6017 9.7281 12.9339 -4.3781 -2.3575 1.829 5.8115 5.9556
         11.7602 11.6648 14.8956 -6.1544 -2.6717 2.6479 3.5028 7.0297 10.071
         11.7603 13.8033 -6.9208 -1.1704 1.7851 4.6972 4.0923 8.5024 10.4058
         16.5646 -5.0422 -0.681 2.4553 2.7624 5.8583 9.7647 8.5024 13.2094
         -5.925 -3.4005 4.0063 3.5458 6.0707 9.275 12.8843 9.6236 15.2427
         0.3918 3.9116 4.6905 5.0126 9.8877 10.9504 14.213 15.7323 -7.1154
         5.3512 4.8912 5.9857 9.3681 12.3194 11.9438 14.8322 -4.8915 1.7766
         5.6502 6.1714 8.6976 11.0897 13.4459 13.5275 -5.9434 -2.0734 5.6249
         9.2996 8.9981 9.1861 11.2791 13.9087 -4.6234 -1.6135 2.1362 5.6249
         10.5304 10.417 10.6243 13.2743 -5.1403 -2.1188 5.7742 4.6512 6.5386
         11.7007 11.5751 12.7756 -4.6857 -2.3404 1.8198 5.7742 5.9524 9.9963
         13.5632 15.8399 -6.1959 -3.063 2.6107 3.389 6.9792 10.0668 10.5594
         15.7167 -6.8441 -1.3231 1.5249 4.7811 4.0354 8.4545 10.4092 10.5594
         -4.349 -0.1812 2.3847 2.5047 5.818 9.5304 8.4545 13.2126 14.5365
         -2.4704 3.9899 3.3903 5.6475 9.3853 12.7505 9.5573 15.2394 -6.094
         4.7884 5.5436 4.8602 9.2662 11.0325 14.1736 15.6533 -7.1108 -0.8459
         6.1554 6.9369 9.1833 11.7281 11.971 14.6501 -4.9332 1.7757 3.9954
         7.8796 9.3856 10.8202 12.7816 13.5812 -6.1106 -2.1258 5.6317 5.6858
         10.8149 10.0527 11.0477 13.3547 -4.6419 -1.6292 2.108 5.6317 8.0603
         12.585 11.7987 13.0239 -5.4095 -2.1485 5.6038 4.6174 6.5384 8.314
         13.8262 13.7047 -4.7575 -2.7484 1.9033 5.6038 5.9048 10.0096 9.0599
         17.7262 -6.0763 -3.1967 2.3161 3.4069 6.787 10.0064 10.5641 11.8948
         -6.3695 -0.5006 1.4397 4.4239 4.0957 8.4074 10.3403 10.5641 13.9393
         1.3043 2.6581 2.4379 5.4463 9.4627 8.4074 13.1355 14.5342 -4.5664
         4.986 3.3347 5.5196 8.9335 12.8455 9.4957 15.1626 -6.0886 -3.1902
         7.1721 5.8635 9.02 10.3905 14.3489 15.4679 -7.0989 -0.845 4.5911
         8.5435 9.0445 11.5236 11.3409 14.7342 -5.1359 1.8057 3.9993 4.7638
         10.8048 11.4134 12.5778 13.0688 -6.1966 -2.1312 5.6498 5.6985 6.2495
         12.4703 11.7446 13.1806 -4.9468 -1.5356 1.9967 5.6498 8.0637 9.3232
         13.9613 14.3726 -5.4877 -2.544 5.607 4.4734 6.565 8.315 9.6674
         15.351 -4.5948 -2.8812 1.6386 5.607 5.8554 10.0433 9.0512 10.4271
         -5.5427 -2.8756 2.2417 3.0281 6.8472 9.9436 10.5907 11.8957 15.6461
         1.1265 1.6092 4.2991 3.7829 8.5351 10.1472 10.5907 13.9457 -6.5362
         3.5658 3.116 5.3834 9.0261 8.5351 12.928 14.5476 -4.5604 0.1901
         4.2978 6.3908 8.7844 12.2207 9.6588 15.0609 -6.0745 -3.1879 4.7466
         7.516 9.5726 10.1635 13.7255 15.4887 -7.1843 -0.8265 4.5938 5.3226
         10.4216 11.763 11.143 14.0208 -5.2344 1.6069 4.0225 4.7741 6.711
         13.7076 13.5185 12.9049 -6.399 -2.0181 5.5142 5.7245 6.2499 9.4311
         13.7747 14.7116 -5.0299 -1.8369 1.9887 5.5142 8.0847 9.3292 10.2417
         16.9047 -5.3726 -2.6528 5.1218 4.477 6.3899 8.3414 9.6598 11.4785
         -3.8393 -2.5754 1.5419 5.1218 5.9747 9.8025 9.0623 10.4239 13.4721
         -1.8099 2.2982 2.8887 6.1359 10.0835 10.4059 11.9275 15.6473 -5.7224
         2.3271 5.1805 3.6882 8.0192 10.2109 10.4059 13.9679 -6.5319 -0.6125
         4.2466 6.2807 8.9188 8.0192 12.9827 14.4545 -4.5446 0.1938 2.9736
         8.0539 10.1401 12.0165 9.005 15.1491 -6.1759 -3.1732 4.7537 4.0583
         11.6204 10.8753 13.4916 14.9399 -7.2695 -0.9554 4.618 5.3237 5.3454
         13.3234 11.7936 13.7918 -5.4862 1.5088 3.8683 4.7962 6.7154 10.2064
         15.7203 13.9318 -6.4142 -2.372 5.3895 5.5351 6.2689 9.4255 11.9703
         17.349 -4.8394 -1.9398 1.7084 5.3895 7.934 9.3598 10.2543 12.0595
         -4.7124 -2.4289 4.951 4.0591 6.2986 8.1613 9.6756 11.4771 13.7749
         -1.4722 2.3698 4.951 5.4878 9.5708 8.9904 10.4575 13.4689 -4.1417
         3.0016 2.8294 5.8619 9.4907 10.2722 11.7161 15.6841 -5.7161 -2.5595
         6.6926 4.596 7.8002 9.5017 10.2722 13.8083 -6.5202 -0.6133 1.8734
         7.7778 8.9677 7.8002 12.253 14.4335 -4.6599 0.2171 2.9814 2.8679
         12.3034 13.3663 8.7364 14.637 -6.2756 -3.2771 4.7713 4.0637 6.2068
         13.0675 14.9485 14.7778 -7.3588 -1.0283 4.459 5.348 5.3394 9.9288
         13.4305 15.4149 -5.5128 1.3222 3.7527 4.6346 6.7436 10.2189 12.5149
         16.0962 -6.3146 -2.5096 5.2062 5.3358 6.1379 9.4501 11.9732 13.7327
         -4.0542 -0.676 1.623 5.2062 7.8289 9.1486 10.2873 12.0619 14.0356
         -1.5061 4.8049 3.9141 6.1388 8.0622 9.576 11.5041 13.7668 -4.9945
         3.7085 4.8049 5.2851 9.3093 9.0257 10.2467 13.4898 -4.1343 -2.1873
         3.7297 5.6083 9.2247 10.0667 11.5982 15.4507 -5.6995 -2.5568 2.822
         6.0243 8.3783 9.2626 10.0667 13.6775 -6.6049 -0.5914 1.8761 4.7961
         10.0592 8.3783 11.9975 14.2992 -4.7705 0.0519 3.0011 2.8688 6.1168
         15.9113 9.7417 14.4722 -6.3821 -3.3512 4.6402 4.0769 6.2082 9.4199
         17.7151 15.4923 -7.6125 -1.1762 4.3516 5.1874 5.3514 9.928 11.1847
         18.4776 -5.4147 0.4718 3.624 4.4664 6.5542 10.2454 12.5223 12.2121
         -5.8586 -1.5595 4.6841 5.1214 6.0687 9.2928 12.0016 13.7382 13.7219
         0.8362 1.7266 4.6841 7.6299 8.9884 10.0555 12.0881 14.0337 -4.4627
         5.884 3.8322 5.4406 7.8773 9.5899 11.3253 13.7829 -4.9898 -1.8995
         5.884 5.8418 8.5567 8.9246 10.1596 13.3541 -4.114 -2.183 1.8753
         7.4111 9.5554 9.3412 11.3918 15.3139 -5.8178 -2.5402 2.8294 3.5244
         10.0628 10.0386 9.3412 13.478 -6.6865 -0.7383 1.8876 4.793 4.1502
         10.0628 12.468 13.7101 -4.9023 -0.0642 2.8617 2.8847 6.1226 9.8007
         12.1193 15.5951 -6.689 -3.4709 4.5145 3.9753 6.2306 9.4136 12.9819
         17.3945 -7.436 -1.7663 4.221 5.0875 5.2745 9.9489 11.1858 14.3191
         -4.8493 1.5395 3.1355 4.2778 6.4226 10.0515 12.5469 12.2196 14.9494
         -0.0497 4.8048 4.5745 5.9178 9.2498 11.801 13.7749 13.724 -5.919
         2.4337 4.8048 6.9412 8.7657 9.8413 11.9049 14.0633 -4.4573 -1.5446
         4.7831 6.3259 7.1207 9.4839 11.2437 13.6772 -4.9757 -1.8919 5.7999
         7.5089 9.035 8.3087 10.0124 13.3043 -4.2587 -2.1641 1.8712 5.7999
         11.6828 9.9139 10.529 15.1419 -5.9354 -2.6585 2.8477 3.5301 7.0178
         12.0643 9.9139 12.8404 -6.7802 -0.8098 1.8057 4.8075 4.148 8.5064
         14.476 14.2088 -5.2853 -0.2543 2.7294 2.7795 6.1537 9.8131 8.5064
         17.7702 -6.478 -3.8985 4.3437 3.8812 6.0792 9.4365 12.9805 9.6297
         -7.1388 -1.2449 3.6785 4.9435 5.2843 9.8037 11.2153 14.3105 15.7264
         3.6956 3.8265 3.7857 6.2586 9.8553 12.3689 12.2497 14.949 -4.9053
         5.5401 4.724 5.3076 9.1045 11.6756 13.5239 13.7477 -5.9105 -2.0632
         5.5401 7.3009 7.9915 9.5958 11.7927 13.8701 -4.4401 -1.5508 2.1293
         7.8027 7.9573 8.831 11.0663 13.6908 -5.0776 -1.868 5.809 4.6419
         10.4 9.0674 9.1957 13.1578 -4.3973 -2.2988 1.8807 5.809 5.9558
         11.1877 11.395 14.339 -6.0549 -2.7428 2.7173 3.5434 7.0272 10.0715
         11.1878 13.3404 -7.0536 -0.9576 1.7417 4.7127 4.1645 8.5044 10.4012
         15.8503 -5.0835 -0.9368 2.5877 2.7179 5.9458 9.8415 8.5044 13.2047
         -6.1037 -3.5576 3.8678 3.7321 5.9903 9.293 13.0114 9.6264 15.2423
         -0.0926 3.9312 4.3382 5.1961 9.7385 11.0148 14.334 15.7316 -7.1164
         4.9491 4.473 5.607 9.6198 12.2206 12.0372 14.9881 -4.8947 1.7725
         5.2926 5.8681 8.3399 11.431 13.3441 13.5853 -5.8905 -2.0696 5.6233
         8.6212 8.6727 8.8753 11.578 13.7871 -4.5633 -1.5358 2.1348 5.6233
         9.7784 9.7664 10.2625 13.5553 -5.1709 -2.0366 5.8281 4.6493 6.5349
         10.9547 10.6825 12.4893 -4.5476 -2.4054 1.8204 5.8281 5.9542 9.9933
         12.8747 15.5321 -6.3961 -2.883 2.5891 3.4417 7.0503 10.0691 10.5563
         15.0275 -6.8848 -1.6089 1.6437 4.6889 4.0579 8.5173 10.408 10.5563
         -4.5986 -0.4112 2.178 2.623 5.7876 9.6405 8.5173 13.2115 14.5339
         -2.7947 3.9874 3.312 5.8421 9.2678 12.8067 9.6441 15.2412 -6.0951
         4.4624 5.1343 4.676 9.5559 10.9025 14.184 15.763 -7.1126 -0.8483
         5.7105 6.4719 8.9763 12.0 11.8686 14.728 -4.8706 1.7767 3.9926
         7.2601 9.2608 10.5361 13.0876 13.4827 -6.0312 -2.0576 5.6292 5.6839
         10.1666 9.4341 10.7628 13.6109 -4.6742 -1.6275 2.1486 5.6292 8.0577
         11.8235 11.2463 12.8153 -5.2874 -2.1874 5.6831 4.6659 6.5393 8.3103
         13.0623 13.2718 -4.986 -2.5616 1.8239 5.6831 5.9674 10.0043 9.0571
         17.0368 -6.1416 -3.3787 2.4482 3.3447 6.8748 10.0862 10.5628 11.8906
         -6.5393 -0.8719 1.32 4.5905 4.0215 8.425 10.4357 10.5628 13.9366
         0.7863 2.4926 2.2161 5.6137 9.4371 8.425 13.2425 14.5361 -4.5677
         4.6177 3.3541 5.2311 9.1465 12.711 9.5187 15.2663 -6.0907 -3.1919
         6.583 5.3896 8.8126 10.6894 14.1741 15.5513 -7.1098 -0.8449 4.588
         7.9753 9.0956 11.2671 11.6274 14.5936 -5.0402 1.7735 3.9979 4.7621
         10.2973 11.1517 12.1522 13.3037 -6.1796 -2.1336 5.6333 5.6934 6.2468
         11.59 11.4229 12.7936 -4.8089 -1.6241 2.0483 5.6333 8.0632 9.3197
         13.1956 13.7248 -5.6398 -2.3643 5.5372 4.5405 6.5369 8.315 9.6642
         14.7697 -4.6815 -3.0833 1.7627 5.5372 5.874 10.0135 9.0559 10.422
         -5.7323 -3.0196 2.0565 3.2007 6.7154 9.9672 10.5647 11.8963 15.6409
         0.5606 1.5298 4.086 3.9278 8.3977 10.2343 10.5647 13.9436 -6.5372
         3.2351 2.7767 4.9615 9.2229 8.3977 13.0217 14.5323 -4.5627 0.1874
         3.9403 5.9601 8.4258 12.511 9.4837 15.1029 -6.0874 -3.1885 4.7451
         6.9283 9.3485 9.831 14.0188 15.4022 -7.1332 -0.846 4.593 5.3194
         9.9143 11.6508 10.8208 14.3524 -5.2181 1.6841 3.9996 4.77 6.7076
         12.8591 13.0513 12.6196 -6.3108 -2.1239 5.5937 5.7025 6.2504 9.427
         13.0476 13.9653 -5.2229 -1.693 1.9534 5.5937 8.0639 9.3273 10.2387
         16.0006 -5.4306 -2.9101 5.342 4.4168 6.4597 8.3138 9.664 11.4746
         -4.1082 -2.7216 1.3876 5.342 5.8452 9.9473 9.0464 10.4261 13.4684
         -2.186 2.2628 2.8013 6.4606 9.9306 10.4972 11.8944 15.6477 -5.7238
         2.0667 4.7589 3.4427 8.2602 10.0777 10.4972 13.9472 -6.5336 -0.6156
         3.845 5.8349 8.6291 8.2602 12.8531 14.4822 -4.559 0.1927 2.9718
         7.4613 9.5421 11.6029 9.3097 15.0327 -6.1156 -3.188 4.7511 4.0571
         10.8724 10.5695 13.1358 15.1917 -7.23 -0.9003 4.5937 5.3237 5.3427
         12.7555 11.4492 13.2865 -5.3759 1.5499 3.9438 4.7772 6.7138 10.2041
         14.9539 13.38 -6.6619 -2.2039 5.4462 5.6573 6.2487 9.429 11.9667
         16.4086 -4.9359 -2.202 1.8347 5.4462 8.0071 9.3302 10.2495 12.0562
         -4.9483 -2.5328 4.8234 4.2488 6.3374 8.2372 9.6554 11.4786 13.7713
         -1.8627 1.9629 4.8234 5.7147 9.6763 8.9876 10.4208 13.4713 -4.1434
         2.7438 2.85 5.7482 9.7734 10.3314 11.8041 15.6462 -5.7186 -2.5614
         6.1567 4.1652 7.6864 9.8218 10.3314 13.8948 -6.531 -0.6125 1.872
         7.2477 8.9254 7.6864 12.5879 14.4396 -4.5915 0.1938 2.9784 2.8657
         11.5344 12.7248 8.5388 14.8745 -6.2294 -3.2266 4.7556 4.0618 6.2038
         12.2785 14.3057 14.4177 -7.264 -0.9972 4.5321 5.3228 5.3426 9.9258
         12.861 14.6319 -5.7926 1.6217 3.8047 4.7369 6.7157 10.2141 12.5122
         15.3162 -6.3754 -2.7079 5.4677 5.4272 6.1914 9.4217 11.973 13.7284
         -4.3343 -1.2789 1.4961 5.4677 7.8753 9.2552 10.258 12.0617 14.0312
         -1.8251 4.8634 3.792 6.39 8.104 9.5881 11.4749 13.7714 -4.9958
         3.2422 4.8634 5.1633 9.6238 9.0054 10.3192 13.4664 -4.1372 -2.1892
         3.3821 5.7349 8.9887 10.3683 11.6487 15.547 -5.7145 -2.5576 2.8203
         5.5124 8.1057 9.0906 10.3683 13.7359 -6.5552 -0.6151 1.8753 4.7934
         9.6601 8.1057 11.7007 14.578 -4.7197 0.1323 2.9835 2.8688 6.1132
         15.0078 9.2552 14.1384 -6.2658 -3.3181 4.7189 4.0649 6.208 9.4159
         16.7382 15.1548 -7.435 -0.9324 4.3992 5.2561 5.3358 9.9293 11.1808
         17.3744 -5.4784 1.0583 3.7731 4.5438 6.6456 10.2228 12.5199 12.2088
         -6.0209 -1.9928 4.9875 5.3684 6.0979 9.336 11.9727 13.7367 13.7185
         0.3364 1.6598 4.9875 7.9481 9.06 10.193 12.0612 14.0355 -4.4643
         5.4804 3.8624 5.9294 8.1784 9.5788 11.392 13.762 -4.9915 -1.9016
         5.4804 5.5866 9.0505 9.1896 10.1951 13.396 -4.1325 -2.1844 1.8731
         6.706 9.3935 9.8112 11.7134 15.3737 -5.7469 -2.557 2.8266 3.5232
         9.4593 9.6757 9.8112 13.7731 -6.6488 -0.6817 1.8764 4.7947 4.1475
         9.4593 12.2441 14.0536 -4.7413 -0.0129 2.9452 2.8681 6.1204 9.7981
         11.2677 15.0567 -6.4764 -3.2936 4.572 4.0354 6.2074 9.4173 12.9775
         16.7047 -7.4695 -1.3873 4.3857 5.1318 5.2833 9.9263 11.1863 14.3146
         -5.0508 0.875 3.5017 4.5004 6.4802 10.1737 12.5241 12.2171 14.9441
         -0.5731 4.8152 4.9215 6.174 9.2661 11.8929 13.7386 13.7239 -5.9205
         2.1761 4.8152 7.3549 9.0738 9.9376 11.9873 14.0312 -4.4593 -1.5476
         4.4291 5.8062 7.6221 9.7383 11.2773 13.6897 -4.9889 -1.8947 5.7985
         6.9023 8.8912 8.6863 10.2382 13.3248 -4.1734 -2.1824 1.8733 5.7985
         10.9014 9.6303 11.1145 15.362 -5.8809 -2.6011 2.8314 3.528 7.0158
         11.3372 9.6303 13.2324 -6.6706 -0.7798 1.8476 4.7906 4.1494 8.5042
         13.7576 13.8295 -5.0352 0.0324 2.7895 2.824 6.1235 9.8082 8.5042
         16.9829 -6.522 -3.6224 4.5757 3.9239 6.1449 9.4095 12.9821 9.6263
         -7.3956 -1.5504 4.0845 5.1413 5.2765 9.8592 11.1846 14.3155 15.7223
         2.1406 3.4605 4.0969 6.4761 9.9447 12.4666 12.2212 14.9503 -4.9071
         4.8134 4.7994 5.6973 9.3666 11.7297 13.6447 13.7235 -5.914 -2.0658
         4.8134 7.1287 8.5065 9.9048 11.8411 13.9405 -4.4562 -1.5478 2.128
         6.7362 7.4391 9.2561 11.3662 13.6799 -5.0194 -1.8901 5.8055 4.6405
         9.2817 8.4662 9.8166 13.4192 -4.3336 -2.2285 1.8686 5.8055 5.9535
         10.156 10.9266 14.9562 -5.9309 -2.7052 2.7959 3.5314 7.0242 10.0688
         10.156 13.052 -6.8699 -0.7135 1.7705 4.7389 4.146 8.5061 10.3982
         14.5874 -5.1183 -0.4992 2.7516 2.7446 6.0503 9.817 8.5061 13.2013
         -6.4261 -3.7324 4.1509 3.9535 6.0285 9.3275 12.9782 9.6283 15.2383
        ],)

    @test parse_clock(str) == DataFrame([
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
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse vc-relax graphene output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/7d5cebcf1250114756b88c6064ebe82e6f8fd835/PW/examples/ESM_example/reference/graphene_bc1_vc-relax.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 2,
        alat = 5.0,
        omega = 675.0316,
        nat = 2,
        ntyp = 1,
        nelec = 8.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 8,
        ecutwfc = 36.0,
        ecutrho = 144.0,
        ecutfock = nothing,
        conv_thr = 1.0e-8,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)",
        nstep = 50,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 63 63 22
         "gvecs" "Min" 4920 4920 1053
         "sticks" "Max" 64 64 23
         "gvecs" "Max" 4924 4924 1065
         "sticks" "Sum" 255 255 91
         "gvecs" "Sum" 19689 19689 4231
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
            0.0 0.0 0.0 0.0138889
            0.0 0.0933154 0.0 0.0277778
            0.0 0.1866307 0.0 0.0277778
            0.0 0.2799461 0.0 0.0277778
            0.0 0.3732615 0.0 0.0277778
            0.0 0.4665769 0.0 0.0277778
            0.0 -0.5598922 0.0 0.0138889
            0.0833333 0.0419919 0.0 0.0277778
            0.0833333 0.1353073 0.0 0.0277778
            0.0833333 0.2286227 0.0 0.0277778
            0.0833333 0.321938 0.0 0.0277778
            0.0833333 0.4152534 0.0 0.0277778
            0.0833333 0.5085688 0.0 0.0277778
            0.0833333 -0.5179003 0.0 0.0277778
            0.0833333 -0.4245849 0.0 0.0277778
            0.0833333 -0.3312696 0.0 0.0277778
            0.0833333 -0.2379542 0.0 0.0277778
            0.0833333 -0.1446388 0.0 0.0277778
            0.0833333 -0.0513235 0.0 0.0277778
            0.1666667 0.0839838 0.0 0.0277778
            0.1666667 0.1772992 0.0 0.0277778
            0.1666667 0.2706146 0.0 0.0277778
            0.1666667 0.36393 0.0 0.0277778
            0.1666667 0.4572453 0.0 0.0277778
            0.1666667 0.5505607 0.0 0.0277778
            0.1666667 -0.4759084 0.0 0.0277778
            0.1666667 -0.382593 0.0 0.0277778
            0.1666667 -0.2892777 0.0 0.0277778
            0.1666667 -0.1959623 0.0 0.0277778
            0.1666667 -0.1026469 0.0 0.0277778
            0.1666667 -0.0093315 0.0 0.0277778
            0.25 0.1259758 0.0 0.0277778
            0.25 0.2192911 0.0 0.0277778
            0.25 0.3126065 0.0 0.0277778
            0.25 0.4059219 0.0 0.0277778
            0.25 0.4992372 0.0 0.0277778
            0.25 0.5925526 0.0 0.0277778
            0.25 -0.4339165 0.0 0.0277778
            0.25 -0.3406011 0.0 0.0277778
            0.25 -0.2472857 0.0 0.0277778
            0.25 -0.1539704 0.0 0.0277778
            0.25 -0.060655 0.0 0.0277778
            0.25 0.0326604 0.0 0.0277778
            0.3333333 0.1679677 0.0 0.0277778
            0.3333333 0.261283 0.0 0.0277778
            0.3333333 0.3545984 0.0 0.0277778
            0.3333333 0.4479138 0.0 0.0277778
            0.3333333 0.5412292 0.0 0.0277778
            0.3333333 0.6345445 0.0 0.0277778
            0.3333333 -0.3919246 0.0 0.0277778
            0.3333333 -0.2986092 0.0 0.0277778
            0.3333333 -0.2052938 0.0 0.0277778
            0.3333333 -0.1119784 0.0 0.0277778
            0.3333333 -0.0186631 0.0 0.0277778
            0.3333333 0.0746523 0.0 0.0277778
            0.4166667 0.2099596 0.0 0.0277778
            0.4166667 0.303275 0.0 0.0277778
            0.4166667 0.3965903 0.0 0.0277778
            0.4166667 0.4899057 0.0 0.0277778
            0.4166667 0.5832211 0.0 0.0277778
            0.4166667 0.6765364 0.0 0.0277778
            0.4166667 -0.3499326 0.0 0.0277778
            0.4166667 -0.2566173 0.0 0.0277778
            0.4166667 -0.1633019 0.0 0.0277778
            0.4166667 -0.0699865 0.0 0.0277778
            0.4166667 0.0233288 0.0 0.0277778
            0.4166667 0.1166442 0.0 0.0277778
            -0.5 -0.2519515 0.0 0.0138889
            -0.5 -0.1586361 0.0 0.0277778
            -0.5 -0.0653208 0.0 0.0277778
            -0.5 0.0279946 0.0 0.0277778
            -0.5 0.12131 0.0 0.0277778
            -0.5 0.2146254 0.0 0.0277778
            -0.5 -0.8118437 0.0 0.0138889
        ],
        cryst = nothing,)

    @test parse_stress(str) == ([-102.58, -88.9, -56.37, 18.73, -0.55, -0.49, 1.27, 0.49, 0.01, 1.64],
        [
         [-0.00087326 -0.00016767 0.0; -0.00016767 -0.00121864 0.0; 0.0 0.0 0.0],
         [-0.00081252 -0.00021559 0.0; -0.00021559 -0.00100054 0.0; 0.0 0.0 0.0],
         [-0.00051278 -0.00017739 0.0; -0.00017739 -0.00063676 0.0; 0.0 0.0 0.0],
         [0.00018175 -3.601e-5 0.0; -3.601e-5 0.00020018 0.0; 0.0 0.0 0.0],
         [1.24e-6 -8.081e-5 0.0; -8.081e-5 -1.236e-5 0.0; 0.0 0.0 0.0],
         [-7.36e-6 -4.911e-5 0.0; -4.911e-5 -2.71e-6 0.0; 0.0 0.0 0.0],
         [1.192e-5 -2.031e-5 0.0; -2.031e-5 1.407e-5 0.0; 0.0 0.0 0.0],
         [4.21e-6 -1.97e-6 0.0; -1.97e-6 5.87e-6 0.0; 0.0 0.0 0.0],
         [-2.9e-7 1.0e-7 0.0; 1.0e-7 4.5e-7 0.0; 0.0 0.0 0.0],
         [1.602e-5 8.4e-7 0.0; 8.4e-7 1.751e-5 0.0; 0.0 0.0 0.0],
        ],
        [
         [-128.46 -24.67 0.0; -24.67 -179.27 0.0; 0.0 0.0 0.0],
         [-119.53 -31.71 0.0; -31.71 -147.18 0.0; 0.0 0.0 0.0],
         [-75.43 -26.09 0.0; -26.09 -93.67 0.0; 0.0 0.0 0.0],
         [26.74 -5.3 0.0; -5.3 29.45 0.0; 0.0 0.0 0.0],
         [0.18 -11.89 0.0; -11.89 -1.82 0.0; 0.0 0.0 0.0],
         [-1.08 -7.22 0.0; -7.22 -0.4 0.0; 0.0 0.0 0.0],
         [1.75 -2.99 0.0; -2.99 2.07 0.0; 0.0 0.0 0.0],
         [0.62 -0.29 0.0; -0.29 0.86 0.0; 0.0 0.0 0.0],
         [-0.04 0.01 0.0; 0.01 0.07 0.0; 0.0 0.0 0.0],
         [2.36 0.12 0.0; 0.12 2.58 0.0; 0.0 0.0 0.0],
        ],)

    @test parse_cell_parameters(str) == []

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.038428924, 0.020283963, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.032025555, 1.501346494, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.052503797, 0.028271169, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.03617185, 1.458627758, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.047848423, 0.02095099, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.014559869, 1.37569809, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.039029671, 0.02192862, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.007333179, 1.392506815, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.041681224, 0.018998183, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.007340697, 1.389283028, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.041394146, 0.020079401, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.003730145, 1.388883645, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.040390405, 0.020188548, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.000341993, 1.389654118, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.040179553, 0.020179171, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.000179347, 1.390083827, 0.0], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.040179553, 0.020179171, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [0.000179347, 1.390083827, 0.0], [1, 1, 1])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 36.0 DavidsonDiagonalization() 0.01 9.7 0.7 4.4 -22.7419296 -22.79397166 0.09226467
             1 2 36.0 DavidsonDiagonalization() 0.00115 4.3 0.7 6.6 -22.75741266 -22.75819232 0.00220438
             1 3 36.0 nothing nothing nothing 0.7 13.2 -22.76081155 -22.7611125 0.00052542
             1 4 36.0 nothing nothing nothing 0.7 18.0 -22.76092694 -22.76104118 0.00018529
             1 5 36.0 nothing nothing nothing 0.7 21.7 -22.76094767 -22.76095229 6.36e-6
             1 6 36.0 nothing nothing nothing 0.7 25.1 -22.76095745 -22.76095813 2.04e-6
             1 7 36.0 nothing nothing nothing 0.7 27.1 -22.76095719 -22.76095763 8.3e-7
             1 8 36.0 nothing nothing nothing 0.7 29.1 -22.76095732 -22.76095737 7.0e-8
             1 9 36.0 nothing nothing nothing 0.7 31.4 nothing nothing nothing
             2 1 36.0 nothing nothing nothing 0.7 42.5 -22.74101355 -20.23788558 0.03490288
             2 2 36.0 nothing nothing nothing 0.7 48.0 -22.80746148 -22.82028254 0.02982682
             2 3 36.0 nothing nothing nothing 0.7 50.3 -22.80293946 -22.80857334 0.00900702
             2 4 36.0 DavidsonDiagonalization() 0.000113 4.6 0.7 52.7 -22.80556447 -22.80557135 3.118e-5
             2 5 36.0 nothing nothing nothing 0.7 59.5 -22.80582575 -22.80582612 3.41e-6
             2 6 36.0 nothing nothing nothing 0.7 62.3 -22.80582609 -22.80582691 1.69e-6
             2 7 36.0 nothing nothing nothing 0.7 64.5 -22.80582667 -22.80582667 2.5e-7
             2 8 36.0 DavidsonDiagonalization() 3.14e-9 3.3 0.7 66.5 -22.80582673 -22.80582674 3.1e-7
             2 9 36.0 DavidsonDiagonalization() 3.14e-9 1.0 0.7 68.1 -22.80582676 -22.80582673 3.0e-7
             2 10 36.0 DavidsonDiagonalization() 3.14e-9 1.0 0.7 69.6 -22.80582686 -22.80582676 3.4e-7
             2 11 36.0 DavidsonDiagonalization() 3.14e-9 1.0 0.7 71.2 -22.80582681 -22.80582686 5.3e-7
             2 12 36.0 DavidsonDiagonalization() 3.14e-9 1.0 0.7 72.9 -22.80582659 -22.80582681 4.6e-7
             2 13 36.0 DavidsonDiagonalization() 3.14e-9 1.2 0.7 74.4 -22.80582659 -22.80582661 1.0e-7
             2 14 36.0 DavidsonDiagonalization() 1.27e-9 1.6 0.7 76.1 -22.8058266 -22.80582659 5.0e-8
             2 15 36.0 DavidsonDiagonalization() 6.28e-10 1.0 0.7 77.6 -22.80582661 -22.80582661 6.0e-8
             2 16 36.0 DavidsonDiagonalization() 6.28e-10 1.0 0.7 79.2 -22.8058266 -22.80582661 6.0e-8
             2 17 36.0 DavidsonDiagonalization() 6.28e-10 1.0 0.7 80.8 -22.80582658 -22.8058266 3.0e-8
             2 18 36.0 DavidsonDiagonalization() 3.16e-10 2.6 0.7 82.5 nothing nothing nothing
             3 1 36.0 nothing nothing nothing 0.7 93.2 -22.7468979 -19.70006135 0.04847438
             3 2 36.0 nothing nothing nothing 0.7 98.7 -22.83816887 -22.85561547 0.04013339
             3 3 36.0 nothing nothing nothing 0.7 101.2 -22.83310771 -22.83973618 0.01191896
             3 4 36.0 nothing nothing nothing 0.7 103.6 -22.8352954 -22.83554913 0.00040865
             3 5 36.0 nothing nothing nothing 0.7 109.0 -22.83604946 -22.83608996 5.904e-5
             3 6 36.0 nothing nothing nothing 0.7 112.7 -22.83605935 -22.83608475 3.967e-5
             3 7 36.0 nothing nothing nothing 0.7 114.9 -22.83606655 -22.83606653 1.1e-7
             3 8 36.0 nothing nothing nothing 0.7 118.3 -22.83606848 -22.83606858 6.1e-7
             3 9 36.0 nothing nothing nothing 0.7 119.8 -22.83606822 -22.83606848 4.3e-7
             3 10 36.0 DavidsonDiagonalization() 1.34e-9 4.1 0.7 122.1 -22.83606833 -22.83606833 2.0e-8
             3 11 36.0 DavidsonDiagonalization() 1.88e-10 2.1 0.7 123.8 -22.83606832 -22.83606833 9.7e-9
             3 12 36.0 DavidsonDiagonalization() 1.21e-10 1.0 0.7 125.3 -22.83606832 -22.83606832 3.8e-9
             3 13 36.0 DavidsonDiagonalization() 4.8e-11 3.0 0.7 127.3 nothing nothing nothing
             4 1 36.0 nothing nothing nothing 0.7 138.0 -22.65852563 -17.71493516 0.10235542
             4 2 36.0 nothing nothing nothing 0.7 143.7 -22.85396198 -22.89387522 0.0930354
             4 3 36.0 DavidsonDiagonalization() 0.00116 6.0 0.7 146.0 -22.84206699 -22.85806784 0.02537891
             4 4 36.0 nothing nothing nothing 0.7 149.0 -22.85024885 -22.85019128 9.07e-5
             4 5 36.0 nothing nothing nothing 0.7 157.0 -22.85071173 -22.85102108 0.00042757
             4 6 36.0 nothing nothing nothing 0.7 161.1 -22.85081769 -22.85082454 2.86e-5
             4 7 36.0 nothing nothing nothing 0.7 163.0 -22.85080973 -22.850819 1.641e-5
             4 8 36.0 nothing nothing nothing 0.7 164.9 -22.85081206 -22.85081234 7.5e-7
             4 9 36.0 nothing nothing nothing 0.7 167.4 -22.85081298 -22.85081298 4.0e-8
             4 10 36.0 nothing nothing nothing 0.7 169.8 -22.850813 -22.850813 4.0e-8
             4 11 36.0 nothing nothing nothing 0.7 171.4 -22.85081305 -22.850813 5.0e-8
             4 12 36.0 DavidsonDiagonalization() 5.34e-10 4.5 0.7 173.8 -22.85081304 -22.85081311 2.0e-7
             4 13 36.0 DavidsonDiagonalization() 5.34e-10 1.1 0.7 175.5 -22.850813 -22.85081305 9.0e-8
             4 14 36.0 DavidsonDiagonalization() 5.34e-10 3.0 0.7 177.7 nothing nothing nothing
             5 1 36.0 nothing nothing nothing 0.7 189.0 2.75183783 -25.01230609 26.97120926
             5 2 36.0 DavidsonDiagonalization() 0.01 11.6 0.7 194.0 -27.42833959 -36.70187397 14.61760936
             5 3 36.0 DavidsonDiagonalization() 0.01 1.9 0.7 195.7 -24.34775088 -27.50787812 4.95822197
             5 4 36.0 DavidsonDiagonalization() 0.01 1.8 0.7 197.4 -23.50565075 -24.39481148 1.69673074
             5 5 36.0 DavidsonDiagonalization() 0.01 1.2 0.7 199.1 -23.10587731 -23.52899921 0.77985757
             5 6 36.0 DavidsonDiagonalization() 0.00975 1.6 0.7 200.8 -23.08144068 -23.12608834 0.3441246
             5 7 36.0 DavidsonDiagonalization() 0.0043 3.8 0.7 202.8 -23.22956281 -23.1018355 0.28523064
             5 8 36.0 DavidsonDiagonalization() 0.00357 1.0 0.7 204.4 -23.06844039 -23.23529512 0.41107317
             5 9 36.0 DavidsonDiagonalization() 0.00357 1.2 0.7 206.2 -23.12987732 -23.07537965 0.24007957
             5 10 36.0 DavidsonDiagonalization() 0.003 1.2 0.7 207.8 -22.97003006 -23.13400087 0.29389488
             5 11 36.0 DavidsonDiagonalization() 0.003 1.8 0.7 209.5 -22.79914779 -22.98954102 0.26734624
             5 12 36.0 nothing nothing nothing 0.7 211.6 -22.85358071 -22.84524303 0.00240209
             5 13 36.0 nothing nothing nothing 0.7 220.4 -22.87236221 -22.87055919 0.04042041
             5 14 36.0 nothing nothing nothing 0.7 225.9 -22.85596002 -22.87242031 0.0451086
             5 15 36.0 nothing nothing nothing 0.7 229.2 -22.84875066 -22.85749748 0.01242644
             5 16 36.0 nothing nothing nothing 0.7 234.1 -22.85324898 -22.85557035 0.00747411
             5 17 36.0 nothing nothing nothing 0.7 236.1 -22.85451365 -22.85337607 0.00355046
             5 18 36.0 DavidsonDiagonalization() 3.0e-5 1.2 0.7 237.7 -22.85162468 -22.85455341 0.005188
             5 19 36.0 DavidsonDiagonalization() 3.0e-5 2.8 0.7 239.6 -22.85215577 -22.85209607 0.00030514
             5 20 36.0 nothing nothing nothing 0.7 242.5 -22.8522129 -22.85223046 0.00032305
             5 21 36.0 nothing nothing nothing 0.7 244.1 -22.85210114 -22.85221703 0.00029516
             5 22 36.0 DavidsonDiagonalization() 3.69e-6 1.1 0.7 245.9 -22.85209229 -22.85210961 0.00011835
             5 23 36.0 DavidsonDiagonalization() 1.48e-6 1.2 0.7 247.5 -22.85206792 -22.85209783 6.805e-5
             5 24 36.0 nothing nothing nothing 0.7 250.0 -22.85210062 -22.8521004 6.0e-7
             5 25 36.0 nothing nothing nothing 0.7 253.1 -22.85210234 -22.85210327 9.1e-7
             5 26 36.0 nothing nothing nothing 0.7 254.7 -22.85210149 -22.85210236 3.2e-7
             5 27 36.0 DavidsonDiagonalization() 4.0e-9 3.3 0.7 256.9 -22.85210189 -22.85210162 4.0e-8
             5 28 36.0 DavidsonDiagonalization() 5.06e-10 2.7 0.7 258.7 -22.85210202 -22.85210189 1.0e-8
             5 29 36.0 DavidsonDiagonalization() 1.71e-10 2.3 0.7 260.5 nothing nothing nothing
             6 1 36.0 nothing nothing nothing 0.7 265.6 -22.85228019 -22.84145323 5.113e-5
             6 2 36.0 nothing nothing nothing 0.7 267.4 -22.85228434 -22.85228576 1.73e-6
             6 3 36.0 DavidsonDiagonalization() 2.16e-8 4.9 0.7 270.1 -22.85228748 -22.85228792 6.1e-7
             6 4 36.0 DavidsonDiagonalization() 7.62e-9 3.0 0.7 272.3 -22.85228763 -22.85228772 1.4e-7
             6 5 36.0 DavidsonDiagonalization() 1.76e-9 2.1 0.7 274.2 -22.85228765 -22.85228765 2.2e-9
             6 6 36.0 DavidsonDiagonalization() 2.7e-11 5.7 0.7 277.2 -22.85228766 -22.85228766 3.6e-9
             6 7 36.0 DavidsonDiagonalization() 2.7e-11 2.3 0.7 279.1 nothing nothing nothing
             7 1 36.0 nothing nothing nothing 0.7 285.6 -22.85224559 -22.74355249 7.98e-5
             7 2 36.0 nothing nothing nothing 0.7 291.4 -22.8523917 -22.852421 7.006e-5
             7 3 36.0 nothing nothing nothing 0.7 293.8 -22.8523814 -22.85239483 2.014e-5
             7 4 36.0 DavidsonDiagonalization() 2.52e-7 3.2 0.7 296.1 -22.85238924 -22.85239056 1.53e-6
             7 5 36.0 DavidsonDiagonalization() 1.91e-8 3.3 0.7 298.2 -22.85238993 -22.85238989 1.0e-8
             7 6 36.0 DavidsonDiagonalization() 1.87e-10 6.5 0.7 301.3 -22.85238992 -22.85239 1.3e-7
             7 7 36.0 DavidsonDiagonalization() 1.87e-10 4.2 0.7 303.8 -22.85238994 -22.85238994 3.3e-9
             7 8 36.0 DavidsonDiagonalization() 4.08e-11 3.8 0.7 305.9 nothing nothing nothing
             8 1 36.0 nothing nothing nothing 0.7 311.7 -22.85239662 -22.90076596 1.404e-5
             8 2 36.0 nothing nothing nothing 0.7 316.6 -22.85242478 -22.85243122 1.31e-5
             8 3 36.0 nothing nothing nothing 0.7 318.3 -22.85242251 -22.85242538 3.92e-6
             8 4 36.0 nothing nothing nothing 0.7 320.5 -22.85242394 -22.85242396 5.0e-8
             8 5 36.0 DavidsonDiagonalization() 6.21e-10 6.0 0.7 323.3 -22.85242407 -22.85242408 2.0e-8
             8 6 36.0 DavidsonDiagonalization() 2.13e-10 4.1 0.7 325.8 -22.85242407 -22.85242409 2.0e-8
             8 7 36.0 DavidsonDiagonalization() 2.13e-10 3.0 0.7 327.7 -22.85242408 -22.85242408 3.1e-10
             8 8 36.0 DavidsonDiagonalization() 3.94e-12 5.1 0.7 330.6 -22.85242408 -22.85242408 5.8e-10
             8 9 36.0 DavidsonDiagonalization() 3.94e-12 1.0 0.7 332.2 -22.85242408 -22.85242408 2.2e-10
             8 10 36.0 DavidsonDiagonalization() 2.77e-12 3.7 0.7 334.5 -22.85242408 -22.85242408 1.5e-10
             8 11 36.0 DavidsonDiagonalization() 1.81e-12 1.5 0.7 336.2 nothing nothing nothing
             9 1 36.0 nothing nothing nothing 0.7 342.6 -22.85241577 -22.8824481 5.48e-6
             9 2 36.0 nothing nothing nothing 0.7 346.4 -22.85242599 -22.85242902 5.34e-6
             9 3 36.0 nothing nothing nothing 0.7 348.1 -22.8524256 -22.85242622 1.6e-6
             9 4 36.0 DavidsonDiagonalization() 2.0e-8 2.8 0.7 349.8 -22.85242569 -22.85242581 4.2e-7
             9 5 36.0 DavidsonDiagonalization() 5.26e-9 3.0 0.7 351.8 -22.85242584 -22.85242582 1.0e-7
             9 6 36.0 DavidsonDiagonalization() 1.29e-9 2.8 0.7 353.6 -22.85242578 -22.85242585 1.2e-7
             9 7 36.0 DavidsonDiagonalization() 1.29e-9 3.0 0.7 355.4 -22.85242581 -22.8524258 6.3e-9
             9 8 36.0 DavidsonDiagonalization() 7.86e-11 3.9 0.7 357.8 -22.85242581 -22.85242581 1.0e-8
             9 9 36.0 DavidsonDiagonalization() 7.86e-11 1.0 0.7 359.3 -22.85242581 -22.85242581 9.3e-9
             9 10 36.0 DavidsonDiagonalization() 7.86e-11 1.0 0.7 360.9 -22.85242581 -22.85242581 8.9e-9
             9 11 36.0 DavidsonDiagonalization() 7.86e-11 1.0 0.7 362.5 -22.85242581 -22.85242581 6.6e-9
             9 12 36.0 DavidsonDiagonalization() 7.86e-11 1.0 0.7 364.0 -22.85242584 -22.85242581 1.0e-8
             9 13 36.0 DavidsonDiagonalization() 7.86e-11 4.2 0.7 366.2 -22.8524258 -22.85242585 1.0e-7
             9 14 36.0 DavidsonDiagonalization() 7.86e-11 4.8 0.7 368.5 -22.85242581 -22.85242581 5.8e-9
             9 15 36.0 DavidsonDiagonalization() 7.25e-11 1.0 0.7 370.1 -22.85242581 -22.85242581 6.4e-10
             9 16 36.0 DavidsonDiagonalization() 7.95e-12 4.0 0.7 372.4 nothing nothing nothing
             10 1 36.0 nothing nothing nothing 0.7 385.2 -22.81064615 -22.95951827 0.22994284
             10 2 36.0 DavidsonDiagonalization() 0.00287 3.2 0.7 387.1 -22.8436412 -22.84401242 0.00233278
             10 3 36.0 nothing nothing nothing 0.7 394.4 -22.85030808 -22.85100307 0.00106604
             10 4 36.0 nothing nothing nothing 0.7 398.6 -22.85052625 -22.8505562 5.107e-5
             10 5 36.0 nothing nothing nothing 0.7 402.7 -22.85055296 -22.85055539 4.46e-6
             10 6 36.0 nothing nothing nothing 0.7 404.9 -22.85055387 -22.85055429 5.4e-7
             10 7 36.0 nothing nothing nothing 0.7 407.3 -22.85055419 -22.85055423 4.0e-8
             10 8 36.0 DavidsonDiagonalization() 4.81e-10 4.0 0.7 409.6 -22.85055421 -22.85055422 4.8e-10
             10 9 36.0 DavidsonDiagonalization() 5.98e-12 5.6 0.7 412.6 -22.85055421 -22.85055421 5.4e-10
             10 10 36.0 DavidsonDiagonalization() 5.98e-12 2.9 0.7 414.3 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-22.76095737,
        -22.76095737,
        5.9e-10,
        nothing,
        [39.24275678, -19.23546571, -6.59178126, -36.17669235],
        nothing,
        0.00022517,)

    @test parse_version(str) == "6.1"

    @test parse_parallel_info(str) == ("Parallel version (MPI & OpenMP)", 1)

    @test parse_fft_dimensions(str) == (19689, (nr1 = 20, nr2 = 20, nr3 = 120))

    @test parse_bands(str) == ([
         0.0 0.0 -0.4278
         0.0 0.1695 0.0
         0.0 -0.5325 0.3527
         0.0 0.0 -0.3235
         0.0933 0.172 0.0
         0.0 -0.4272 0.3555
         0.0 0.0 -0.2193
         0.1866 0.1745 0.0
         0.0 -0.3219 0.3582
         0.0 0.0 -0.1151
         0.2799 0.177 0.0
         0.0 -0.2166 0.361
         0.0 0.0 -0.0109
         0.3733 0.1795 0.0
         0.0 -0.1113 0.3637
         0.0 0.0 0.0933
         0.4666 0.182 0.0
         0.0 -0.006 0.4581
         0.0 0.0 0.2469
         -0.5599 0.2768 0.0
         0.0 0.1489 0.4609
         0.0833 0.0 0.3511
         0.042 0.2793 0.0
         0.0 0.2541 0.4637
         0.0833 0.0 0.4554
         0.1353 0.2818 0.0
         0.0 0.3594 0.4664
         0.0833 0.0 0.5596
         0.2286 0.2843 0.0
         0.0 0.4647 0.4692
         0.0833 0.0 0.6638
         0.3219 0.2868 0.0
         0.0 0.57 0.4719
         0.0833 0.0 0.768
         0.4153 0.2893 0.0
         0.0 0.6753 0.4416
         0.0833 0.0 -0.3784
         0.5086 0.2618 0.0
         0.0 -0.4828 0.4443
         0.0833 0.0 -0.2742
         -0.5179 0.2643 0.0
         0.0 -0.3776 0.4471
         0.0833 0.0 -0.1699
         -0.4246 0.2668 0.0
         0.0 -0.2723 0.4499
         0.0833 0.0 -0.0657
         -0.3313 0.2693 0.0
         0.0 -0.167 0.4526
         0.0833 0.0 0.0385
         -0.238 0.2718 0.0
         0.0 -0.0617 0.4554
         0.0833 0.0 0.1427
         -0.1446 0.2743 0.0
         0.0 0.0436 -0.5498
         0.0833 0.0 -0.2963
         -0.0513 0.369 0.0
         0.0 0.1985 -0.547
         0.1667 0.0 -0.1921
         0.084 0.3715 0.0
         0.0 0.3038 -0.5442
         0.1667 0.0 -0.0879
         0.1773 0.374 0.0
         0.0 0.409 -0.5415
         0.1667 0.0 0.0163
         0.2706 0.3765 0.0
         0.0 0.5143 -0.5387
         0.1667 0.0 0.1206
         0.3639 0.379 0.0
         0.0 0.6196 -0.536
         0.1667 0.0 0.2248
         0.4572 0.3815 0.0
         0.0 0.7249 -0.5663
         0.1667 0.0 -0.9216
         0.5506 0.3541 0.0
         0.0 -0.4332 0.0
         0.1667 0.0 0.0
         -0.4759 0.3566 0.0
         0.0 -0.3279 0.0029
         0.1667 0.0 0.1042
         -0.3826 0.3591 0.0
         0.0 -0.2227 0.0059
         0.1667 0.0 0.2084
         -0.2893 0.3616 0.0
         0.0 -0.1174 0.0088
         0.1667 0.0 0.3125
         -0.196 0.3641 0.0
         0.0 -0.0121 0.0117
         0.1667 0.0 0.4167
         -0.1026 0.3666 0.0
         0.0 0.0932 0.0146
         0.1667 0.0 0.5209
         -0.0093 0.4613 0.0
         0.0 0.2481 -0.0176
         0.25 0.0 -0.6251
         0.126 0.4638 0.0
         0.0 0.3534 0.0917
         0.25 0.0 0.0495
         0.2193 0.4663 0.0
         0.0 0.4587 0.0946
         0.25 0.0 0.1537
         0.3126 0.4688 0.0
         0.0 0.5639 0.0975
         0.25 0.0 0.2579
         0.4059 0.4713 0.0
         0.0 0.6692 0.1005
         0.25 0.0 0.3621
         0.4992 0.4738 0.0
         0.0 0.7745 0.1034
         0.25 0.0 0.4662
         0.5926 0.4463 0.0
         0.0 -0.3836 0.1063
         0.25 0.0 0.5704
         -0.4339 0.4488 0.0
         0.0 -0.2783 0.0741
         0.25 0.0 -0.5755
         -0.3406 0.4513 0.0
         0.0 -0.173 0.077
         0.25 0.0 -0.4714
         -0.2473 0.4538 0.0
         0.0 -0.0678 0.08
         0.25 0.0 -0.3672
         -0.154 0.4563 0.0
         0.0 0.0375 0.0829
         0.25 0.0 -0.263
         -0.0607 0.4588 0.0
         0.0 0.1428 0.0858
         0.25 0.0 -0.1588
         0.0327 -0.5536 0.0
         0.0 -0.2977 0.0887
         0.3333 0.0 -0.0547
         0.168 -0.5511 0.0
         0.0 -0.1924 0.1833
         0.3333 0.0 0.099
         0.2613 -0.5486 0.0
         0.0 -0.0871 0.1863
         0.3333 0.0 0.2032
         0.3546 -0.5461 0.0
         0.0 0.0181 0.1892
         0.3333 0.0 0.3074
         0.4479 -0.5436 0.0
         0.0 0.1234 0.1921
         0.3333 0.0 0.4116
         0.5412 -0.5411 0.0
         0.0 0.2287 0.1951
         0.3333 0.0 0.5157
         0.6345 -0.5686 0.0
         0.0 -0.9294 0.198
         0.3333 0.0 0.6199
         -0.3919 0.0 0.0
         0.0 0.0 0.1658
         0.3333 0.0 -0.526
         -0.2986 0.0023 0.0
         0.0 0.1041 0.1687
         0.3333 0.0 -0.4218
         -0.2053 0.0047 0.0
         0.0 0.2081 0.1716
         0.3333 0.0 -0.3177
         -0.112 0.007 0.0
         0.0 0.3122 0.1745
         0.3333 0.0 -0.2135
         -0.0187 0.0093 0.0
         0.0 0.4162 0.1775
         0.3333 0.0 -0.1093
         0.0747 0.0116 0.0
         0.0 0.5203 0.1804
         0.4167 0.0 -0.0051
         0.21 -0.014 0.0
         0.0 -0.6244 0.275
         0.4167 0.0 0.1486
         0.3033 0.0914 0.0
         0.0 0.0489 0.2779
         0.4167 0.0 0.2527
         0.3966 0.0937 0.0
         0.0 0.153 0.2809
         0.4167 0.0 0.3569
         0.4899 0.096 0.0
         0.0 0.257 0.2838
         0.4167 0.0 0.4611
         0.5832 0.0983 0.0
         0.0 0.3611 0.2867
         0.4167 0.0 0.5653
         0.6765 0.1007 0.0
         0.0 0.4652 0.2896
         0.4167 0.0 0.6694
         -0.3499 0.103 0.0
         0.0 0.5692 0.2574
         0.4167 0.0 -0.4765
         -0.2566 0.0774 0.0
         0.0 -0.5754 0.2604
         0.4167 0.0 -0.3723
         -0.1633 0.0797 0.0
         0.0 -0.4714 0.2633
         0.4167 0.0 -0.2681
         -0.07 0.082 0.0
         0.0 -0.3673 0.2662
         0.4167 0.0 -0.164
         0.0233 0.0844 0.0
         0.0 -0.2633 0.2691
         0.4167 0.0 -0.0598
         0.1166 0.0867 0.0
         0.0 -0.1592 0.2721
         -0.5 0.0 0.0444
         -0.252 0.089 0.0
         0.0 -0.0551 0.3667
         -0.5 0.0 0.1981
         -0.1586 0.1827 0.0
         0.0 0.0978 0.3696
         -0.5 0.0 0.3023
         -0.0653 0.185 0.0
         0.0 0.2019 0.3725
         -0.5 0.0 0.4064
         0.028 0.1874 0.0
         0.0 0.306 0.3755
         -0.5 0.0 0.5106
         0.1213 0.1897 0.0
         0.0 0.41 0.3784
         -0.5 0.0 0.6148
         0.2146 0.192 0.0
         0.0 0.5141 0.3813
         -0.5 0.0 0.719
         -0.8118 0.1943 0.0
         0.0 0.6181 0.3491
         0.0 0.0 -0.427
         0.0 0.1687 0.0
         0.0 -0.5265 0.352
         0.0004 0.0 -0.3228
         0.0962 0.1711 0.0
         0.0 -0.4225 0.355
         0.0008 0.0 -0.2186
         0.1925 0.1734 0.0
         0.0 -0.3184 0.3579
         0.0012 0.0 -0.1144
         0.2887 0.1757 0.0
         0.0 -0.2143 0.3608
         0.0016 0.0 -0.0103
         0.385 0.1781 0.0
         0.0 -0.1103 0.3637
         0.0021 0.0 0.0939
         0.4812 0.1804 0.0
         0.0 -0.0062 0.4583
         -0.0025 0.0 0.2476
         -0.5775 0.2741 0.0
         0.0 0.1467 0.4613
         0.0854 0.0 0.3518
         0.0437 0.2764 0.0
         0.0 0.2508 0.4642
         0.0858 0.0 0.456
         0.1399 0.2787 0.0
         0.0 0.3549 0.4671
         0.0862 0.0 0.5601
         0.2362 0.281 0.0
         0.0 0.4589 0.4701
         0.0866 0.0 0.6643
         0.3324 0.2834 0.0
         0.0 0.563 0.473
         0.087 0.0 0.7685
         0.4287 0.2857 0.0
         0.0 0.667 0.4408
         0.0874 0.0 -0.3775
         0.5249 0.2601 0.0
         0.0 -0.4776 0.4437
         0.0829 0.0 -0.2733
         -0.5338 0.2624 0.0
         0.0 -0.3736 0.4466
         0.0833 0.0 -0.1691
         -0.4376 0.2648 0.0
         0.0 -0.2695 0.4495
         0.0837 0.0 -0.0649
         -0.3413 0.2671 0.0
         0.0 -0.1654 0.4525
         0.0841 0.0 0.0393
         -0.2451 0.2694 0.0
         0.0 -0.0614 0.4554
         0.0846 0.0 0.1434
         -0.1488 0.2717 0.0
         0.0 0.0427 -0.55
         0.085 0.0 -0.2971
         -0.0526 0.3654 0.0
         0.0 0.1957 -0.5471
         0.1708 0.0 -0.193
         0.0874 0.3677 0.0
         0.0 0.2997 -0.5441
         0.1712 0.0 -0.0888
         0.1836 0.3701 0.0
         0.0 0.4038 -0.5412
         0.1716 0.0 0.0154
         0.2799 0.3724 0.0
         0.0 0.5078 -0.5383
         0.172 0.0 0.1196
         0.3761 0.3747 0.0
         0.0 0.6119 -0.5354
         0.1724 0.0 0.2238
         0.4724 0.377 0.0
         0.0 0.716 -0.5676
         0.1728 0.0 -0.9222
         0.5686 0.3515 0.0
         0.0 -0.4287 0.0
         0.1683 0.0 0.0
         -0.4901 0.3538 0.0
         0.0 -0.3246 0.0029
         0.1687 0.0 0.1041
         -0.3939 0.3561 0.0
         0.0 -0.2206 0.0059
         0.1691 0.0 0.2083
         -0.2976 0.3584 0.0
         0.0 -0.1165 0.0088
         0.1695 0.0 0.3124
         -0.2014 0.3608 0.0
         0.0 -0.0125 0.0118
         0.1699 0.0 0.4166
         -0.1051 0.3631 0.0
         0.0 0.0916 0.0147
         0.1703 0.0 0.5207
         -0.0089 0.4568 0.0
         0.0 0.2446 -0.0177
         0.2561 0.0 -0.6249
         0.131 0.4591 0.0
         0.0 0.3486 0.0917
         0.2566 0.0 0.0495
         0.2273 0.4614 0.0
         0.0 0.4527 0.0946
         0.257 0.0 0.1537
         0.3235 0.4637 0.0
         0.0 0.5568 0.0976
         0.2574 0.0 0.2578
         0.4198 0.4661 0.0
         0.0 0.6608 0.1005
         0.2578 0.0 0.362
         0.516 0.4684 0.0
         0.0 0.7649 0.1035
         0.2582 0.0 0.4661
         0.6123 0.4428 0.0
         0.0 -0.3798 0.1064
         0.2537 0.0 0.5702
         -0.4465 0.4451 0.0
         0.0 -0.2757 0.074
         0.2541 0.0 -0.5753
         -0.3502 0.4475 0.0
         0.0 -0.1717 0.0769
         0.2545 0.0 -0.4712
         -0.254 0.4498 0.0
         0.0 -0.0676 0.0799
         0.2549 0.0 -0.3671
         -0.1577 0.4521 0.0
         0.0 0.0365 0.0828
         0.2553 0.0 -0.2629
         -0.0615 0.4544 0.0
         0.0 0.1405 0.0858
         0.2557 0.0 -0.1588
         0.0348 -0.5481 0.0
         0.0 -0.2935 0.0887
         0.3415 0.0 -0.0546
         0.1747 -0.5458 0.0
         0.0 -0.1894 0.1833
         0.3419 0.0 0.0991
         0.271 -0.5435 0.0
         0.0 -0.0854 0.1863
         0.3423 0.0 0.2032
         0.3672 -0.5411 0.0
         0.0 0.0187 0.1892
         0.3428 0.0 0.3073
         0.4635 -0.5388 0.0
         0.0 0.1227 0.1922
         0.3432 0.0 0.4115
         0.5597 -0.5365 0.0
         0.0 0.2268 0.1951
         0.3436 0.0 0.5156
         0.656 -0.5621 0.0
         0.0 -0.9179 0.1981
         0.339 0.0 0.6198
         -0.4028 0.0 0.0
         0.0 0.0 0.1656
         0.3395 0.0 -0.5258
         -0.3065 0.0025 0.0
         0.0 0.1041 0.1686
         0.3399 0.0 -0.4217
         -0.2103 0.005 0.0
         0.0 0.2082 0.1715
         0.3403 0.0 -0.3175
         -0.114 0.0075 0.0
         0.0 0.3123 0.1745
         0.3407 0.0 -0.2134
         -0.0178 0.0101 0.0
         0.0 0.4164 0.1774
         0.3411 0.0 -0.1092
         0.0785 0.0126 0.0
         0.0 0.5205 0.1804
         0.4269 0.0 -0.0051
         0.2184 -0.0151 0.0
         0.0 -0.6245 0.275
         0.4273 0.0 0.1486
         0.3147 0.0914 0.0
         0.0 0.0491 0.2779
         0.4277 0.0 0.2527
         0.4109 0.094 0.0
         0.0 0.1532 0.2809
         0.4281 0.0 0.3569
         0.5072 0.0965 0.0
         0.0 0.2573 0.2838
         0.4286 0.0 0.461
         0.6034 0.099 0.0
         0.0 0.3614 0.2868
         0.429 0.0 0.5652
         0.6997 0.1015 0.0
         0.0 0.4655 0.2897
         0.4244 0.0 0.6693
         -0.3591 0.104 0.0
         0.0 0.5696 0.2573
         0.4248 0.0 -0.4763
         -0.2628 0.0763 0.0
         0.0 -0.5754 0.2602
         0.4253 0.0 -0.3721
         -0.1666 0.0789 0.0
         0.0 -0.4714 0.2632
         0.4257 0.0 -0.268
         -0.0703 0.0814 0.0
         0.0 -0.3673 0.2661
         0.4261 0.0 -0.1639
         0.0259 0.0839 0.0
         0.0 -0.2632 0.2691
         0.4265 0.0 -0.0597
         0.1222 0.0864 0.0
         0.0 -0.1591 0.272
         -0.5123 0.0 0.0444
         -0.2621 0.0889 0.0
         0.0 -0.055 0.3666
         -0.5119 0.0 0.1981
         -0.1658 0.1829 0.0
         0.0 0.0982 0.3696
         -0.5115 0.0 0.3022
         -0.0696 0.1854 0.0
         0.0 0.2023 0.3725
         -0.511 0.0 0.4064
         0.0267 0.1879 0.0
         0.0 0.3064 0.3755
         -0.5106 0.0 0.5105
         0.1229 0.1904 0.0
         0.0 0.4105 0.3784
         -0.5102 0.0 0.6147
         0.2192 0.1929 0.0
         0.0 0.5146 0.3814
         -0.5148 0.0 0.7188
         -0.8396 0.1955 0.0
         0.0 0.6187 0.3489
         0.0 0.0 -0.4268
         0.0 0.1678 0.0
         0.0 -0.5263 0.3519
         0.0011 0.0 -0.3226
         0.0997 0.1703 0.0
         0.0 -0.4222 0.3548
         0.0022 0.0 -0.2185
         0.1993 0.1728 0.0
         0.0 -0.3182 0.3578
         0.0033 0.0 -0.1143
         0.299 0.1753 0.0
         0.0 -0.2141 0.3607
         0.0045 0.0 -0.0102
         0.3987 0.1778 0.0
         0.0 -0.11 0.3637
         0.0056 0.0 0.094
         0.4984 0.1804 0.0
         0.0 -0.0059 0.4583
         -0.0067 0.0 0.2476
         -0.598 0.2743 0.0
         0.0 0.1473 0.4612
         0.088 0.0 0.3518
         0.0459 0.2768 0.0
         0.0 0.2514 0.4642
         0.0891 0.0 0.4559
         0.1455 0.2793 0.0
         0.0 0.3555 0.4671
         0.0903 0.0 0.5601
         0.2452 0.2819 0.0
         0.0 0.4596 0.4701
         0.0914 0.0 0.6642
         0.3449 0.2844 0.0
         0.0 0.5637 0.473
         0.0925 0.0 0.7684
         0.4445 0.2869 0.0
         0.0 0.6678 0.4406
         0.0936 0.0 -0.3772
         0.5442 0.2592 0.0
         0.0 -0.4772 0.4435
         0.0813 0.0 -0.2731
         -0.5522 0.2617 0.0
         0.0 -0.3731 0.4465
         0.0824 0.0 -0.1689
         -0.4525 0.2642 0.0
         0.0 -0.2691 0.4494
         0.0836 0.0 -0.0648
         -0.3528 0.2668 0.0
         0.0 -0.165 0.4524
         0.0847 0.0 0.0393
         -0.2532 0.2693 0.0
         0.0 -0.0609 0.4553
         0.0858 0.0 0.1435
         -0.1535 0.2718 0.0
         0.0 0.0432 -0.5499
         0.0869 0.0 -0.2972
         -0.0538 0.3657 0.0
         0.0 0.1964 -0.547
         0.176 0.0 -0.193
         0.0917 0.3683 0.0
         0.0 0.3005 -0.544
         0.1772 0.0 -0.0889
         0.1914 0.3708 0.0
         0.0 0.4046 -0.5411
         0.1783 0.0 0.0153
         0.291 0.3733 0.0
         0.0 0.5087 -0.5381
         0.1794 0.0 0.1194
         0.3907 0.3758 0.0
         0.0 0.6128 -0.5352
         0.1805 0.0 0.2236
         0.4904 0.3783 0.0
         0.0 0.7169 -0.5676
         0.1816 0.0 -0.922
         0.5901 0.3506 0.0
         0.0 -0.4281 0.0
         0.1693 0.0 0.0
         -0.5063 0.3532 0.0
         0.0 -0.324 0.0029
         0.1705 0.0 0.1041
         -0.4067 0.3557 0.0
         0.0 -0.2199 0.0059
         0.1716 0.0 0.2083
         -0.307 0.3582 0.0
         0.0 -0.1159 0.0088
         0.1727 0.0 0.3124
         -0.2073 0.3607 0.0
         0.0 -0.0118 0.0118
         0.1738 0.0 0.4166
         -0.1076 0.3632 0.0
         0.0 0.0923 0.0147
         0.1749 0.0 0.5207
         -0.008 0.4572 0.0
         0.0 0.2455 -0.0177
         0.2641 0.0 -0.6249
         0.1376 0.4597 0.0
         0.0 0.3496 0.0917
         0.2652 0.0 0.0495
         0.2372 0.4622 0.0
         0.0 0.4537 0.0946
         0.2663 0.0 0.1537
         0.3369 0.4647 0.0
         0.0 0.5578 0.0976
         0.2674 0.0 0.2578
         0.4366 0.4672 0.0
         0.0 0.6619 0.1005
         0.2685 0.0 0.362
         0.5362 0.4698 0.0
         0.0 0.766 0.1035
         0.2697 0.0 0.4661
         0.6359 0.4421 0.0
         0.0 -0.379 0.1064
         0.2574 0.0 0.5702
         -0.4605 0.4446 0.0
         0.0 -0.2749 0.074
         0.2585 0.0 -0.5753
         -0.3608 0.4471 0.0
         0.0 -0.1708 0.0769
         0.2596 0.0 -0.4712
         -0.2611 0.4496 0.0
         0.0 -0.0668 0.0799
         0.2607 0.0 -0.3671
         -0.1615 0.4521 0.0
         0.0 0.0373 0.0828
         0.2618 0.0 -0.2629
         -0.0618 0.4547 0.0
         0.0 0.1414 0.0858
         0.263 0.0 -0.1588
         0.0379 -0.5486 0.0
         0.0 -0.2946 0.0887
         0.3521 0.0 -0.0546
         0.1834 -0.5461 0.0
         0.0 -0.1905 0.1833
         0.3532 0.0 0.0991
         0.2831 -0.5436 0.0
         0.0 -0.0864 0.1863
         0.3543 0.0 0.2032
         0.3828 -0.5411 0.0
         0.0 0.0176 0.1892
         0.3554 0.0 0.3073
         0.4824 -0.5386 0.0
         0.0 0.1217 0.1922
         0.3566 0.0 0.4115
         0.5821 -0.536 0.0
         0.0 0.2258 0.1951
         0.3577 0.0 0.5156
         0.6818 -0.5637 0.0
         0.0 -0.9192 0.1981
         0.3454 0.0 0.6198
         -0.4146 0.0 0.0
         0.0 0.0 0.1656
         0.3465 0.0 -0.5258
         -0.315 0.0028 0.0
         0.0 0.1042 0.1686
         0.3476 0.0 -0.4217
         -0.2153 0.0055 0.0
         0.0 0.2084 0.1715
         0.3487 0.0 -0.3175
         -0.1156 0.0083 0.0
         0.0 0.3127 0.1745
         0.3499 0.0 -0.2134
         -0.0159 0.011 0.0
         0.0 0.4169 0.1774
         0.351 0.0 -0.1092
         0.0837 0.0138 0.0
         0.0 0.5211 0.1804
         0.4401 0.0 -0.0051
         0.2293 -0.0166 0.0
         0.0 -0.6253 0.275
         0.4412 0.0 0.1486
         0.3289 0.0916 0.0
         0.0 0.0494 0.2779
         0.4423 0.0 0.2527
         0.4286 0.0944 0.0
         0.0 0.1536 0.2809
         0.4435 0.0 0.3569
         0.5283 0.0971 0.0
         0.0 0.2578 0.2838
         0.4446 0.0 0.461
         0.6279 0.0999 0.0
         0.0 0.362 0.2868
         0.4457 0.0 0.5652
         0.7276 0.1027 0.0
         0.0 0.4663 0.2897
         0.4334 0.0 0.6693
         -0.3688 0.1054 0.0
         0.0 0.5705 0.2573
         0.4345 0.0 -0.4763
         -0.2691 0.0751 0.0
         0.0 -0.5759 0.2602
         0.4357 0.0 -0.3721
         -0.1694 0.0778 0.0
         0.0 -0.4717 0.2632
         0.4368 0.0 -0.268
         -0.0698 0.0806 0.0
         0.0 -0.3675 0.2661
         0.4379 0.0 -0.1639
         0.0299 0.0833 0.0
         0.0 -0.2633 0.2691
         0.439 0.0 -0.0597
         0.1296 0.0861 0.0
         0.0 -0.159 0.272
         -0.5281 0.0 0.0444
         -0.2751 0.0889 0.0
         0.0 -0.0548 0.3666
         -0.527 0.0 0.1981
         -0.1754 0.1833 0.0
         0.0 0.0988 0.3696
         -0.5259 0.0 0.3022
         -0.0758 0.186 0.0
         0.0 0.203 0.3725
         -0.5248 0.0 0.4064
         0.0239 0.1888 0.0
         0.0 0.3072 0.3755
         -0.5237 0.0 0.5105
         0.1236 0.1915 0.0
         0.0 0.4114 0.3784
         -0.5226 0.0 0.6147
         0.2233 0.1943 0.0
         0.0 0.5156 0.3814
         -0.5348 0.0 0.7188
         -0.8731 0.1971 0.0
         0.0 0.6199 0.3489
         0.0 0.0 -0.4268
         0.0 0.1667 0.0
         0.0 -0.5265 0.3519
         0.0025 0.0 -0.3226
         0.1053 0.1695 0.0
         0.0 -0.4223 0.3548
         0.005 0.0 -0.2185
         0.2106 0.1722 0.0
         0.0 -0.3181 0.3578
         0.0075 0.0 -0.1143
         0.3158 0.175 0.0
         0.0 -0.2139 0.3607
         0.01 0.0 -0.0102
         0.4211 0.1777 0.0
         0.0 -0.1097 0.3637
         0.0125 0.0 0.094
         0.5264 0.1805 0.0
         0.0 -0.0054 0.4583
         -0.015 0.0 0.2476
         -0.6317 0.2749 0.0
         0.0 0.1482 0.4612
         0.0923 0.0 0.3518
         0.0496 0.2776 0.0
         0.0 0.2524 0.4642
         0.0948 0.0 0.4559
         0.1549 0.2804 0.0
         0.0 0.3566 0.4671
         0.0973 0.0 0.5601
         0.2602 0.2832 0.0
         0.0 0.4608 0.4701
         0.0998 0.0 0.6642
         0.3655 0.2859 0.0
         0.0 0.565 0.473
         0.1022 0.0 0.7684
         0.4707 0.2887 0.0
         0.0 0.6692 0.4406
         0.1047 0.0 -0.3772
         0.576 0.2583 0.0
         0.0 -0.4771 0.4435
         0.0773 0.0 -0.2731
         -0.5821 0.2611 0.0
         0.0 -0.3729 0.4465
         0.0798 0.0 -0.1689
         -0.4768 0.2638 0.0
         0.0 -0.2687 0.4494
         0.0823 0.0 -0.0648
         -0.3715 0.2666 0.0
         0.0 -0.1645 0.4524
         0.0848 0.0 0.0393
         -0.2662 0.2694 0.0
         0.0 -0.0603 0.4553
         0.0873 0.0 0.1435
         -0.1609 0.2721 0.0
         0.0 0.0439 -0.5499
         0.0898 0.0 -0.2972
         -0.0557 0.3665 0.0
         0.0 0.1975 -0.547
         0.1845 0.0 -0.193
         0.0992 0.3693 0.0
         0.0 0.3018 -0.544
         0.187 0.0 -0.0889
         0.2045 0.372 0.0
         0.0 0.406 -0.5411
         0.1895 0.0 0.0153
         0.3098 0.3748 0.0
         0.0 0.5102 -0.5381
         0.192 0.0 0.1194
         0.4151 0.3776 0.0
         0.0 0.6144 -0.5352
         0.1945 0.0 0.2236
         0.5204 0.3803 0.0
         0.0 0.7186 -0.5676
         0.197 0.0 -0.922
         0.6256 0.3499 0.0
        ],
        [
         -22.5973 -0.9516 -19.3138 -1.6904 -24.0212 -0.9793 -19.6948 -1.6876
         -11.2683 -0.0565 -14.9535 3.9088 -12.2977 -0.032 -14.9339 3.2785
         -9.3879 0.15 -13.4338 4.4947 -7.6032 0.1662 -13.5928 4.424
         -6.9496 1.8177 -7.2705 5.7637 -7.447 1.8338 -7.2837 5.4305
         -1.1366 -22.4122 -2.1074 -19.2148 -1.2049 -23.3777 -1.585 -19.1011
         -0.6766 -10.8808 4.558 -17.2455 -0.2606 -11.5121 4.84 -17.0959
         -0.2849 -10.1929 5.1864 -11.2714 -0.0628 -9.8294 5.5827 -11.4642
         -0.0886 -8.9937 5.5399 -6.955 1.608 -8.6356 6.1809 -6.9917
         -22.4643 -0.3952 -20.4726 -2.5672 -23.856 -0.305 -20.9414 -2.5314
         -11.1122 0.5018 -13.7831 4.6167 -12.0969 0.6512 -13.9331 3.9402
         -9.87 0.7199 -12.2365 5.7948 -8.2146 0.8573 -12.0549 5.7308
         -7.3174 2.0609 -8.4885 6.7868 -7.8264 2.5104 -8.6655 7.1032
         -0.9597 -22.0142 -0.5804 -18.2684 -0.9794 -22.8877 0.0604 -18.2211
         -0.3655 -11.0859 2.9071 -17.5984 -0.0316 -10.924 3.0722 -17.4056
         -0.1023 -10.4283 3.5916 -12.6588 0.1682 -10.9197 4.0272 -12.8033
         0.0981 -9.8688 4.1801 -6.1943 1.8351 -9.5463 4.3326 -6.266
         -22.067 0.1478 -21.3825 -3.0538 -23.362 0.3704 -21.9009 -3.0263
         -11.1159 1.042 -12.7634 6.4539 -11.4998 1.3323 -12.9655 5.7368
         -10.6484 1.2767 -10.8779 7.3389 -9.8273 1.5496 -10.493 7.3304
         -8.1153 1.7942 -9.4872 7.4403 -8.6508 3.1791 -9.7698 7.4629
         -0.4219 -21.3423 0.8355 -18.2924 -0.303 -22.0764 1.5756 -18.2354
         0.1971 -12.2082 1.6378 -16.5185 0.6535 -12.207 1.726 -16.3999
         0.4447 -11.0969 2.1661 -14.2349 0.8613 -10.9676 2.727 -14.2904
         0.6566 -9.6697 2.8712 -5.7892 2.5138 -9.9626 2.9427 -5.856
         -21.4121 0.447 -22.0099 -3.0934 -22.5449 1.4958 -22.5499 -3.1494
         -12.7545 1.0727 -11.7495 6.454 -11.9168 2.0577 -11.8855 6.2671
         -9.8922 2.2279 -10.1898 7.0831 -10.525 2.4627 -10.53 6.9431
         -8.933 3.6199 -10.0126 8.5742 -9.5555 2.7054 -9.5904 7.9623
         0.3296 -20.4136 0.7828 -18.2734 0.8221 -20.9571 0.8248 -18.2358
         0.3886 -13.7915 1.7025 -16.3335 1.7754 -13.9253 1.7768 -16.3978
         0.5583 -11.9041 1.7974 -14.4576 2.0159 -12.0527 2.0165 -14.2921
         1.3695 -8.6452 1.9922 -5.8227 2.9019 -8.6742 2.9314 -5.8574
         -20.5151 -1.1134 -22.3413 -2.9795 -21.4155 0.0652 -22.877 -3.1473
         -14.4962 2.3663 -11.0764 6.3932 -14.1235 3.0619 -10.9161 6.267
         -9.5809 3.5578 -10.5627 7.1125 -10.3071 4.0189 -10.8833 6.9432
         -8.8828 3.9119 -9.7001 8.1391 -9.2216 4.3229 -9.5921 7.9628
         -1.4629 -19.2755 0.3352 -18.2025 0.4944 -19.5587 0.3751 -18.2236
         0.4101 -15.5007 1.276 -17.4078 2.3839 -15.868 1.3362 -17.4011
         1.7561 -12.2853 1.5199 -12.9704 3.4072 -12.6249 1.5541 -12.8051
         2.632 -7.4708 2.647 -6.2521 3.5082 -7.1998 3.1837 -6.2697
         -19.4236 -2.6812 -22.3715 -2.7976 -19.9985 -1.9043 -22.8774 -3.0223
         -16.1412 3.6102 -11.4234 5.9024 -16.242 5.0424 -10.9163 5.7385
         -9.9881 4.0941 -10.5918 7.239 -10.7912 5.8673 -10.8939 7.3278
         -7.7453 4.7884 -9.2861 7.554 -7.7501 5.8931 -9.5803 7.4631
         -3.1408 -18.2819 0.2965 -19.1018 -1.713 -18.2217 0.3747 -19.1018
         0.3012 -16.7504 1.2353 -17.0873 3.2904 -17.4129 1.3359 -17.0934
         3.3535 -12.342 1.4815 -11.6298 4.4178 -12.8036 1.5537 -11.4661
         4.2404 -6.656 2.5228 -6.9851 5.4254 -6.2621 3.1833 -6.995
         -18.5782 -3.6539 -21.0723 -2.2843 -18.6342 -3.02 -21.4201 -2.5276
         -17.1781 3.5286 -14.323 4.0959 -17.8483 5.7683 -14.1115 3.9427
         -10.1264 5.8594 -9.9879 5.9352 -10.9572 7.339 -10.3254 5.7298
         -7.0601 5.9742 -9.105 6.8812 -6.8639 7.4779 -9.227 7.1021
         -4.0175 -19.0328 -0.4113 -20.0375 -2.8611 -19.535 0.5179 -19.9976
         0.2522 -15.8901 2.082 -16.1778 3.2225 -15.8993 2.3875 -16.228
         5.2179 -12.0898 2.7319 -10.9851 6.7653 -12.6146 3.4109 -10.8182
         5.2634 -7.2361 3.3401 -7.7621 6.8507 -7.1857 3.5139 -7.7591
         -22.4643 -2.9454 -20.3852 -1.426 -23.857 -1.9118 -20.6266 -1.6853
         -11.1221 3.5514 -15.169 3.4422 -12.0984 5.0728 -15.1158 3.2813
         -9.5467 4.3743 -10.8044 4.5776 -8.2216 5.848 -11.1236 4.4233
         -7.7765 5.1224 -8.353 5.5848 -7.8203 5.9167 -8.3501 5.4296
         -0.9629 -20.2082 -1.4479 -20.6945 -0.981 -20.9377 -0.6984 -20.6228
         -0.1045 -14.2781 3.0215 -15.0142 -0.0332 -13.9505 3.4966 -15.1097
         -0.0959 -11.533 3.3674 -11.2977 0.1666 -12.0425 4.136 -11.1329
         0.1065 -8.4051 4.3376 -8.3686 1.8336 -8.6593 4.5149 -8.3506
         -22.2133 -1.3772 -19.4396 -0.428 -23.53 0.0663 -19.5456 -0.6976
         -10.8384 2.6501 -15.8056 3.6079 -11.7023 3.0862 -15.8849 3.4955
         -10.2439 3.5285 -12.2525 4.2952 -9.2406 4.0411 -12.6143 4.1174
         -8.5598 3.6968 -7.3422 4.6301 -8.5183 4.3468 -7.1942 4.5132
         -0.6236 -21.1886 -2.4861 -21.037 -0.5332 -22.062 -1.9168 -20.9406
         0.2405 -12.6461 4.3813 -13.8258 0.4207 -12.1883 5.0502 -13.9356
         0.3782 -10.7011 5.012 -12.1775 0.6254 -10.9963 5.8508 -12.0505
         0.5086 -9.476 5.2493 -8.6947 2.2836 -9.9503 5.8979 -8.666
         -21.702 0.285 -18.2983 0.2947 -22.8772 1.5139 -18.2378 0.0528
         -11.6659 1.2801 -16.2344 3.1495 -10.9179 2.0737 -16.4036 3.0678
         -10.2537 2.4333 -13.8176 4.1103 -10.9066 2.4797 -14.2872 4.023
         -9.3607 2.7978 -6.1473 4.4017 -9.5711 2.7234 -5.8535 4.3284
         0.0748 -21.9123 -3.5399 -21.0547 0.3651 -22.8781 -3.1493 -20.9406
         0.9489 -11.1931 5.7518 -14.0427 1.3277 -10.9152 6.2733 -13.937
         0.9849 -10.2951 5.9815 -11.9127 1.5464 -10.8531 6.946 -12.0495
         1.1899 -9.7689 7.0269 -8.6969 3.1768 -9.6213 7.9848 -8.6656
         -20.9411 0.2858 -17.2042 0.2052 -21.9066 0.3825 -16.98 0.0525
         -13.3891 1.413 -16.6957 3.1293 -12.9708 1.3439 -16.9766 3.068
         -10.0135 1.5591 -14.7682 4.0941 -10.4801 1.5615 -15.2934 4.0232
         -9.3907 1.7547 -4.9595 4.382 -9.7749 3.1913 -4.4669 4.3286
         -0.461 -22.3615 -4.5813 -20.7473 1.5525 -23.3729 -4.4623 -20.6227
         1.0142 -10.8133 5.4672 -15.2529 1.7091 -11.5076 5.933 -15.1113
         1.3112 -10.0318 6.8292 -10.9246 2.7097 -9.7465 8.566 -11.1311
         2.0139 -9.1649 7.7758 -8.3715 2.9282 -8.7259 8.5706 -8.3491
         -19.9557 -0.3265 -17.9179 -0.6396 -20.6347 -0.2989 -18.2335 -0.6987
         -15.1282 0.572 -16.2664 3.5552 -15.1053 0.6571 -16.3985 3.4951
         -10.4497 0.7879 -14.1956 4.7199 -11.1201 0.8633 -14.2975 4.1148
         -8.3083 1.2871 -5.9199 4.8769 -8.355 2.5163 -5.8551 4.5136
         -2.2079 -22.5287 -3.5635 -18.6971 -0.7118 -23.5401 -3.143 -18.6374
         1.0364 -11.0105 5.3765 -18.0443 3.4761 -11.7083 6.2702 -17.8325
         2.5907 -9.6041 5.6828 -10.7703 4.1434 -9.1819 6.9521 -10.9822
         3.4216 -8.9438 6.6245 -6.8237 4.4949 -8.5624 7.9629 -6.8759
         -18.8504 -0.5532 -19.118 -2.8989 -19.1157 -0.5277 -19.5414 -2.83
         -16.6612 0.3436 -15.6203 3.8659 -17.0911 0.4256 -15.8784 3.211
         -10.6293 0.557 -12.8212 6.8602 -11.4453 0.6288 -12.6269 6.7739
         -7.227 2.1351 -7.1054 6.9334 -6.989 2.2875 -7.1947 6.8624
         -3.6857 -21.7144 -2.3523 -19.2014 -2.5523 -22.5612 -1.9089 -19.1015
         0.8931 -12.1283 3.339 -17.2606 3.9582 -11.9262 5.0543 -17.0956
         4.3137 -10.0943 4.84 -11.2684 5.703 -10.5382 5.8278 -11.4644
         5.1194 -9.6537 5.4964 -6.9534 7.0788 -9.5355 6.4181 -6.9922
         -18.7186 0.5593 -20.1424 -2.5657 -19.0885 0.8194 -20.6234 -2.5309
         -16.8534 1.0117 -14.8195 4.6054 -17.1243 1.7717 -15.1075 3.94
         -10.5305 1.5367 -11.495 5.8134 -11.4356 2.0112 -11.1381 5.7302
         -7.0837 1.6917 -8.1498 6.7972 -6.9756 2.9218 -8.35 7.1027
         -3.8795 -21.1982 -1.1805 -19.6409 -2.566 -21.9159 -0.6899 -19.5422
         0.8124 -12.9849 1.7966 -16.0212 3.9476 -12.9974 3.4991 -15.8804
         4.4958 -10.5264 3.3796 -12.4726 5.7425 -10.4456 4.1168 -12.6202
         5.3732 -9.5208 4.6516 -7.1764 7.1171 -9.7803 4.5179 -7.1947
         -19.8211 0.0742 -20.9106 -1.8671 -20.611 1.5645 -21.4179 -1.9166
         -15.4503 1.2721 -13.8456 5.1377 -15.1461 1.7164 -14.1012 5.0489
         -10.1469 2.3981 -10.7509 6.0145 -11.0926 2.7192 -10.3409 5.8283
         -8.0993 2.5237 -8.9594 6.4489 -8.3341 2.9335 -9.2267 5.8942
         -2.5582 -20.4214 -0.0838 -19.7785 -0.7309 -20.96 0.526 -19.6954
         0.799 -13.7554 1.3375 -14.9358 3.5071 -13.9518 2.3895 -14.9343
         2.7308 -11.9365 2.3089 -13.6225 4.1128 -12.0201 3.4131 -13.5885
         3.6231 -8.6637 3.258 -7.2778 4.5253 -8.677 3.5008 -7.2853
         -20.8386 -1.092 -21.3951 -1.4354 -21.8882 0.0558 -21.9017 -1.5904
         -13.8226 2.3612 -12.8588 4.9399 -13.0122 3.0583 -12.9645 4.8337
         -9.4928 3.5562 -10.7661 5.7016 -10.4411 4.0162 -10.4939 5.5775
         -9.2131 4.0007 -9.4793 6.2732 -9.7563 4.3192 -9.7704 6.1748
         -0.8396 -19.414 0.7574 -19.6002 1.5416 -19.7165 1.5753 -19.5422
         0.7472 -14.6314 1.6317 -15.7987 1.7336 -14.9289 1.7252 -15.8789
         1.2844 -13.2043 2.2465 -12.7862 2.7341 -13.5803 2.7261 -12.6218
         2.1681 -7.5754 2.8571 -7.1909 2.9522 -7.2942 2.942 -7.1963
         -21.6307 -2.3842 -21.5845 -1.6913 -22.8641 -1.5817 -22.0643 -1.915
         -12.2098 3.7855 -12.5818 5.1873 -10.928 4.8219 -12.1968 5.0484
         -10.1238 4.5023 -10.7036 6.0129 -10.9038 5.5745 -10.9856 5.8291
         -8.638 5.0596 -9.6813 6.0324 -9.5461 6.1638 -9.9546 5.894
         0.1801 -18.2724 1.0027 -19.1242 0.3826 -18.2557 1.5021 -19.1019
         0.4425 -15.9106 1.3786 -17.0624 1.3442 -16.3837 2.0549 -17.0931
         0.8222 -13.6786 2.3208 -11.635 1.5635 -14.2949 2.4678 -11.4661
         1.0974 -6.4158 2.5919 -6.9888 3.1916 -5.8585 2.7115 -6.9951
         -22.1714 -3.6532 -21.4758 -2.286 -23.5223 -3.1389 -21.9027 -2.5275
         -10.7895 5.1532 -13.3582 4.1137 -11.6936 6.2778 -12.9759 3.9428
         -10.76 5.2707 -10.0932 5.9034 -9.1811 6.9288 -10.4797 5.7295
         -7.8015 5.869 -9.5548 6.8627 -8.5752 7.9954 -9.7707 7.1018
         -0.5575 -17.9283 0.5066 -18.6133 -0.5227 -18.2193 1.5699 -18.6397
         -0.1118 -16.3727 1.5309 -17.8515 0.4307 -16.4226 1.724 -17.8279
         0.3005 -13.5474 2.7502 -11.1542 0.6357 -14.2925 2.7249 -10.9835
         0.5039 -6.1179 3.5245 -6.8734 2.2936 -5.8432 2.9407 -6.8798
         -22.4505 -3.9481 -19.7952 -2.5829 -23.8544 -3.1464 -20.0014 -2.8261
         -11.0959 5.1205 -16.2498 3.3801 -12.0955 6.2796 -16.2342 3.2144
         -9.795 5.4292 -10.4215 6.9928 -8.0872 6.9827 -10.8091 6.7711
         -7.422 6.0626 -7.7837 7.093 -7.9546 7.9893 -7.7576 6.8634
         -0.9391 -19.0568 -2.3903 -24.0202 -0.9774 -19.684 -1.6877 -24.0279
         -0.3415 -15.2218 2.4847 -12.297 -0.0298 -14.9417 3.2994 -12.3014
         -0.0837 -12.9342 3.92 -7.5711 0.1701 -13.6006 4.4246 -7.5225
         0.1154 -7.1746 5.2234 -7.4816 1.837 -7.2724 5.4314 -7.5172
         -22.0686 -2.7492 -19.0456 -1.2044 -23.3662 -1.5789 -19.1051 -1.2259
         -10.6912 4.0328 -16.9562 -0.2598 -11.5059 4.8636 -17.1001 -0.3092
         -9.955 4.2797 -11.0709 -0.0615 -9.8504 5.6008 -11.4576 -0.1433
         -9.5512 4.9054 -7.1159 1.6087 -8.6267 6.2043 -6.992 1.5496
         -0.4216 -20.1446 -3.116 -23.855 -0.3092 -20.9346 -2.532 -23.8625
         0.4302 -14.0362 3.1346 -12.0967 0.6473 -13.9191 3.9624 -12.1006
         0.588 -11.8827 5.0475 -8.1352 0.8552 -12.0805 5.7313 -8.1636
         0.6992 -8.32 5.9895 -7.9103 2.5073 -8.6569 7.1038 -7.8675
         -21.7046 -1.3162 -18.2721 -0.9787 -22.8803 0.0774 -18.2239 -0.9999
         -10.7865 2.2915 -17.2048 -0.0311 -10.9361 3.09 -17.4103 -0.0798
         -10.37 2.7252 -12.3604 0.1693 -10.9227 4.0439 -12.7983 0.0877
         -10.2897 3.9237 -6.4248 1.8358 -9.5376 4.3506 -6.2655 1.7772
         0.0783 -21.0032 -3.5774 -23.361 0.3606 -21.898 -3.0263 -23.3677
         0.9038 -12.9352 4.8455 -11.5013 1.3227 -12.9395 5.7603 -11.5032
         1.0572 -10.7329 6.4426 -9.7745 1.5419 -10.5275 7.3323 -9.7917
         1.2189 -9.2656 6.5243 -8.7118 3.1689 -9.7647 7.4639 -8.681
         -21.0887 0.0384 -18.0422 -0.3033 -22.074 1.6022 -18.2373 -0.3219
         -12.2706 0.9697 -16.5615 0.6533 -12.2093 1.7387 -16.4031 0.6069
         -11.0465 1.5316 -13.7384 0.8619 -10.9644 2.7394 -14.2884 0.7812
         -9.5998 2.6912 -5.9417 2.5133 -9.9665 2.9559 -5.8549 2.4568
         0.1131 -21.5986 -3.7369 -22.5437 1.4792 -22.5505 -3.1473 -22.5491
         0.9345 -11.8646 5.8434 -11.8811 2.0272 -11.8513 6.2722 -11.8902
         1.8011 -10.0322 6.4317 -10.53 2.4469 -10.5284 6.9462 -10.5281
         1.9581 -9.9384 6.6852 -9.6022 2.6908 -9.6301 7.9844 -9.5815
         -20.2386 0.6416 -17.9269 0.8209 -20.9607 0.8329 -18.2366 0.8058
         -13.9561 0.724 -16.3485 1.7723 -13.9166 1.7859 -16.3954 1.7327
         -11.5502 1.1516 -14.0984 2.0156 -12.0492 2.0247 -14.2975 1.9367
         -8.6533 1.8524 -5.8471 2.9191 -8.6827 2.9598 -5.8561 2.9162
         -1.4039 -21.9178 -3.6846 -21.4145 0.0447 -22.881 -3.1428 -21.4177
         2.0642 -11.0542 5.3985 -14.098 3.0371 -10.918 6.2702 -14.1044
         2.617 -10.3073 5.7381 -10.3448 3.9949 -10.8527 6.9474 -10.331
         2.9802 -9.9268 6.5201 -9.2333 4.3003 -9.6221 7.9638 -9.225
         -19.1985 0.2781 -17.8509 0.5214 -19.5688 0.3787 -18.2218 0.517
         -15.6076 1.4042 -17.0879 2.3813 -15.8561 1.3403 -17.4014 2.3717
         -11.803 1.5828 -12.9909 3.4072 -12.6131 1.5578 -12.8108 3.3674
         -7.5332 1.7855 -6.1989 3.4941 -7.2086 3.1881 -6.2681 3.4714
         -2.9804 -21.9562 -3.4207 -19.9995 -1.9221 -22.8846 -3.0165 -19.9985
         2.3755 -11.2989 3.2618 -16.2199 5.0071 -10.9213 5.7396 -16.2282
         3.7983 -10.3596 6.3714 -10.824 5.8586 -10.9169 7.3354 -10.8143
         4.4116 -9.6323 6.9783 -7.7763 6.375 -9.5491 7.4707 -7.7566
         -18.2477 0.2277 -18.6723 -1.6687 -18.2182 0.3742 -19.1008 -1.6865
         -16.8604 1.1051 -16.7299 3.2819 -17.421 1.3358 -17.093 3.2774
         -11.767 1.3558 -11.8142 4.4155 -12.7819 1.5533 -11.4718 4.4105
         -6.6254 1.7545 -6.8628 5.4228 -6.2545 3.1824 -6.9936 5.3916
         -4.1195 -20.7607 -2.8569 -18.6607 -3.054 -21.4329 -2.5211 -18.6383
         2.2036 -14.0952 1.6102 -17.8024 5.7511 -14.1363 3.9437 -17.8314
         5.5497 -10.2167 5.4988 -10.9885 7.3069 -10.2843 5.7363 -10.9802
         5.8995 -9.059 6.2803 -6.9088 7.466 -9.2365 7.1086 -6.8759
         -18.891 -0.7939 -19.535 -2.7973 -19.5163 0.5083 -19.9981 -2.8289
         -16.1983 1.6255 -15.8801 3.2165 -15.9313 2.3809 -16.226 3.2097
         -11.4185 1.9198 -11.2203 6.7563 -12.583 3.4053 -10.8242 6.7635
         -7.0977 3.0727 -7.5723 6.8535 -7.1689 3.5576 -7.7583 6.8521
         -3.5952 -20.1465 -2.0628 -23.8565 -1.9549 -20.6433 -1.678 -23.8625
         1.9843 -14.8684 1.0405 -12.0986 5.0785 -15.1336 3.2827 -12.1006
         4.073 -10.982 4.2355 -8.2025 5.8256 -11.0878 4.4279 -8.1652
         4.889 -8.4132 5.543 -7.8428 6.4429 -8.362 5.4345 -7.8659
         -19.9916 -1.6866 -20.1567 -0.9809 -20.9184 -0.7058 -20.6244 -0.9999
         -14.7989 2.2153 -14.8391 -0.0331 -13.9914 3.4877 -15.1061 -0.0798
         -10.7532 2.7979 -11.461 0.1672 -12.0036 4.2265 -11.1389 0.0877
         -8.2385 3.9891 -8.138 1.8337 -8.6423 4.5054 -8.3506 1.7768
         -2.1221 -19.3057 -1.2354 -23.5303 0.0257 -19.5661 -0.6897 -23.5323
         1.5867 -15.4363 1.8411 -11.7055 3.093 -15.8942 3.498 -11.7013
         2.4591 -12.3263 3.359 -9.1469 4.0478 -12.5858 4.1185 -9.1875
         3.3427 -7.5315 4.6337 -8.6201 4.3553 -7.2069 4.5168 -8.5612
         -20.9102 -2.6101 -20.4985 -0.5355 -22.0427 -1.9178 -20.9431 -0.5476
         -13.3713 3.7332 -13.9616 0.4186 -12.2169 5.0349 -13.9314 0.3785
         -9.8348 4.084 -12.0207 0.6241 -10.9736 5.8897 -12.0553 0.5501
         -9.2779 4.7784 -8.4639 2.2814 -9.9337 5.9517 -8.6667 2.2308
         -0.5309 -18.309 -0.7199 -22.8782 1.5208 -18.2617 0.0598 -22.8757
         0.8061 -15.8321 2.8841 -10.9256 2.04 -16.4015 3.0701 -10.9139
         1.1761 -13.724 3.7174 -10.8271 2.4865 -14.2691 4.0254 -10.8893
         2.2727 -6.4749 4.1475 -9.6612 2.7319 -5.8641 4.3306 -9.5822
         -21.5891 -3.5847 -20.5497 0.3601 -22.8595 -3.1427 -20.9438 0.3563
         -12.0353 5.2454 -14.3037 1.3229 -10.8998 6.2846 -13.9402 1.2909
         -10.0733 5.2605 -11.5742 1.5426 -10.8432 6.9242 -12.0452 1.4751
         -8.9116 5.9263 -8.5109 3.1718 -9.6426 8.0829 -8.6667 3.1271
         -0.0076 -17.3361 -0.8576 -21.9086 0.3884 -16.9984 0.0567 -21.9007
         0.2386 -16.4426 2.82 -12.9208 1.3499 -16.9581 3.0693 -12.9695
         0.9063 -14.3771 4.0814 -10.5406 1.5692 -15.299 4.0247 -10.4837
         1.1249 -5.3355 5.1329 -9.7892 3.1991 -4.4677 4.3298 -9.7684
         -22.0119 -4.6476 -20.3109 1.5876 -23.3559 -4.4535 -20.6264 1.5659
         -10.8765 5.1113 -15.3055 1.7018 -11.4943 5.939 -15.1159 1.7072
         -10.5829 5.4034 -10.7012 2.7047 -9.713 8.5519 -11.1242 2.682
         -8.4597 7.1085 -8.2743 2.922 -8.7752 8.6053 -8.3503 2.8618
         -0.3445 -17.6999 -1.5445 -20.6384 -0.2952 -18.2136 -0.6973 -20.6239
         -0.1354 -16.3854 3.1202 -15.0706 0.6609 -16.3997 3.4964 -15.1098
         0.5168 -13.7882 3.2969 -11.1631 0.869 -14.3234 4.1347 -11.1286
         0.7267 -5.844 4.443 -8.3801 2.5212 -5.8401 4.5149 -8.3486
         -22.1721 -4.1186 -18.6796 -0.6633 -23.5254 -3.1383 -18.6407 -0.6987
         -10.7903 3.8478 -17.5892 3.4632 -11.6981 6.2716 -17.8378 3.4822
         -10.0781 5.3673 -10.5703 4.1275 -9.167 6.9876 -10.975 4.111
         -8.7593 6.1937 -6.987 4.4848 -8.5925 7.895 -6.876 4.4759
         -0.5625 -18.8417 -3.42 -19.1274 -0.5275 -19.5271 -2.8315 -19.1019
         0.293 -15.6536 2.4067 -17.0572 0.4262 -15.8679 3.2325 -17.0948
         0.4595 -12.5919 6.0377 -11.4781 0.6311 -12.6603 6.7743 -11.4622
         0.546 -6.9798 6.1222 -7.0319 2.289 -7.1809 6.8629 -6.9917
         -21.4228 -2.9517 -18.9666 -2.4887 -22.5542 -1.8935 -19.1049 -2.5302
         -11.5809 1.9639 -17.0561 3.9475 -11.9481 5.0805 -17.1004 3.9392
         -10.4192 4.551 -11.0359 5.6927 -10.5393 5.759 -11.4577 5.7189
         -10.0033 5.1862 -7.0637 7.0673 -9.5149 6.4449 -6.9924 6.6561
         0.4646 -19.8141 -3.176 -19.0906 0.8085 -20.6144 -2.5315 -19.1022
         0.7157 -14.8168 3.1102 -17.1043 1.7589 -15.0879 3.9618 -17.0945
         0.736 -11.4027 5.1723 -11.4623 2.0024 -11.1771 5.7315 -11.4624
         1.4546 -7.9823 6.1088 -7.0075 2.8847 -8.3394 7.104 -6.9921
         -20.9582 -1.8054 -19.2755 -2.5181 -21.9149 -0.6653 -19.5452 -2.5297
         -12.4071 0.5367 -16.0571 3.9403 -13.0099 3.5166 -15.8848 3.9391
         -11.1816 3.1781 -12.1335 5.7469 -10.4281 4.057 -12.6151 5.7186
         -9.5074 4.3915 -7.187 7.1222 -9.7883 4.5369 -7.1949 6.6557
         -0.0491 -20.5454 -2.6534 -20.6065 1.5311 -21.4138 -1.9151 -20.6241
         1.0599 -13.8456 4.6295 -15.134 1.6975 -14.0735 5.0501 -15.1092
         1.2809 -10.7433 4.919 -11.1157 2.703 -10.3838 5.8493 -11.1292
         2.036 -8.7663 5.4624 -8.3428 2.9165 -9.2202 5.8974 -8.349
         -20.2596 -0.7372 -19.3456 -0.712 -20.9661 0.5559 -19.6977 -0.6978
         -13.2119 0.1223 -15.1932 3.5135 -13.9507 2.4029 -14.9336 3.4817
         -12.3578 2.1621 -13.1426 4.1018 -12.0064 3.4268 -13.5902 4.1107
         -8.7461 3.3397 -7.2065 4.5315 -8.692 3.4482 -7.2852 4.4783
         -1.1319 -21.0115 -2.3034 -21.8841 0.0361 -21.902 -1.5859 -21.9009
         2.1071 -12.856 4.5128 -13.0 3.03 -12.9317 4.8364 -12.9684
         2.5706 -10.8218 5.2128 -10.4612 3.9899 -10.5346 5.5806 -10.4849
         2.9806 -9.2828 5.6486 -9.7566 4.293 -9.7684 6.1774 -9.7688
         -19.3594 0.0776 -19.143 1.5441 -19.7307 1.6028 -19.5435 1.5669
         -14.4561 0.9901 -15.6572 1.738 -14.9112 1.7336 -15.8758 1.7069
         -13.045 1.5236 -12.7552 2.7429 -13.5717 2.7341 -12.6276 2.6811
         -7.759 2.6818 -7.0724 2.9573 -7.3142 2.9506 -7.1958 2.8611
         -2.4109 -21.2028 -2.4463 -22.861 -1.5904 -22.069 -1.9086 -22.8758
         3.3042 -12.489 3.3892 -10.9195 4.7825 -12.1946 5.052 -10.9141
         4.0165 -10.8739 4.7927 -10.9013 5.5422 -10.9857 5.8306 -10.8871
         4.1117 -9.5069 5.5127 -9.5583 6.127 -9.9568 5.8982 -9.5846
         -18.3389 0.348 -18.6846 0.3861 -18.2764 1.5049 -19.1019 0.3562
         -15.9213 1.2629 -16.7253 1.3472 -16.3639 2.0675 -17.0917 1.2908
         -13.142 2.4159 -11.8067 1.5676 -14.2791 2.4712 -11.4722 1.475
         -6.6417 2.8394 -6.8527 3.1931 -5.8723 2.7144 -6.994 3.1276
         -3.7613 -21.1175 -2.8813 -23.5201 -3.1499 -21.9116 -2.5208 -23.5324
         3.6285 -13.1987 1.6235 -11.6907 6.263 -13.0043 3.9443 -11.7015
         4.9246 -10.3165 5.4729 -9.2005 6.8733 -10.4401 5.7348 -9.1834
         5.4115 -9.4304 6.2848 -8.5558 7.9799 -9.7769 7.1072 -8.5656
         -17.8652 -0.031 -18.2342 -0.52 -18.1935 1.5655 -18.6377 -0.5478
         -16.6657 1.3873 -17.4088 0.4337 -16.4618 1.7217 -17.8286 0.3786
         -12.8269 2.3109 -11.3805 0.639 -14.2633 2.7244 -10.9897 0.55
         -6.0388 2.3786 -6.7343 2.2968 -5.8241 2.9388 -6.8785 2.2305
         -4.4784 -19.6143 -3.1132 -23.8531 -3.1866 -20.018 -2.8194 -23.8625
         3.1079 -15.8899 0.9373 -12.0936 6.2764 -16.255 3.2158 -12.1007
         5.3012 -10.5848 6.5092 -8.1816 6.9984 -10.767 6.7778 -8.1614
         5.6647 -7.9365 6.5688 -7.8599 7.9684 -7.7692 6.8702 -7.8698
         -18.8922 -2.4838 -24.1805 -0.976 -19.6588 -1.6995 -24.0255 -1.0
         -15.7234 1.4353 -12.4213 -0.028 -14.9758 3.3865 -12.301 -0.0798
         -12.1114 3.5959 -7.6982 0.1722 -13.5771 4.415 -7.5262 0.0878
         -7.011 4.8308 -7.1426 1.839 -7.2509 5.4229 -7.5194 1.7775
         -3.4467 -19.0041 -1.2158 -23.3671 -1.6144 -19.1271 -1.2057 -23.3677
         2.1045 -16.4641 -0.2727 -11.5087 4.8775 -17.1136 -0.2623 -11.5033
         3.9963 -11.1958 -0.0855 -9.8479 5.6107 -11.4191 -0.0654 -9.7932
         4.7401 -7.3952 1.5902 -8.6324 6.2198 -7.005 1.6057 -8.6794
         -19.9168 -3.0773 -24.0118 -0.3115 -20.9107 -2.5407 -23.86 -0.322
         -14.6088 2.0569 -12.2153 0.645 -13.9162 4.0563 -12.1001 0.6068
         -11.0857 4.6152 -8.2013 0.8536 -12.0959 5.7142 -8.1673 0.7812
         -8.1399 5.4908 -7.6519 2.5049 -8.6362 7.0868 -7.8704 2.4562
         -2.0655 -18.3962 -0.9852 -22.8828 0.0534 -18.2468 -0.98 -22.8757
         0.7888 -16.5914 -0.0382 -10.9312 3.1025 -17.4198 -0.0331 -10.9139
         2.5504 -12.3973 0.1507 -10.8964 4.0549 -12.7658 0.1657 -10.8903
         3.4143 -6.7924 1.8221 -9.5829 4.3647 -6.276 1.8331 -9.581
         -20.7255 -3.4973 -23.5076 0.3539 -21.8759 -3.0301 -23.3653 0.3563
         -13.498 3.6836 -11.6027 1.3163 -12.9161 5.8591 -11.5027 1.2909
         -10.0548 5.7964 -9.7306 1.5364 -10.5666 7.3165 -9.795 1.4752
         -9.0722 5.9434 -8.6133 3.1623 -9.7461 7.4522 -8.6836 3.1274
         -0.7121 -18.0368 -0.2933 -22.0784 1.5911 -18.2432 -0.3027 -22.0628
         -0.2831 -16.1868 0.6632 -12.1276 1.7488 -16.4206 0.6529 -12.1994
         1.4278 -13.633 0.8597 -11.0514 2.7471 -14.2686 0.8596 -10.9788
         2.3021 -6.2528 2.5171 -9.9823 2.9672 -5.8576 2.5125 -9.9524
         -21.2876 -3.7907 -22.6732 1.4685 -22.531 -3.1449 -22.5468 1.4845
         -12.3928 5.3641 -11.8397 2.0633 -11.8155 6.2848 -11.893 2.047
         -9.7427 5.3849 -10.6026 2.436 -10.5131 6.9515 -10.5278 2.4242
         -9.5219 5.9579 -9.5645 2.6815 -9.6834 8.0796 -9.5847 2.6317
         -0.485 -17.7737 0.8578 -20.9674 0.8392 -18.2315 0.8238 -20.9419
         0.4529 -16.2611 1.815 -13.8496 1.7932 -16.3799 1.7758 -13.9369
         0.6552 -13.8435 2.041 -12.1189 2.0329 -14.3254 2.0157 -12.0455
         1.5243 -5.9256 3.2011 -8.7083 2.953 -5.8488 2.9215 -8.6647
         -21.5925 -4.0233 -21.5196 0.0964 -22.8647 -3.1331 -21.4157 0.0515
         -11.3794 3.8709 -14.0884 3.0221 -10.9072 6.27 -14.1061 3.0535
         -10.126 5.3342 -10.3398 3.9781 -10.8138 6.9609 -10.3347 3.9855
         -9.7621 6.1584 -9.2654 4.2868 -9.6754 7.8979 -9.2248 4.2494
         0.2253 -17.5591 0.7412 -19.5801 0.3805 -18.1987 0.5166 -19.5427
         0.4503 -17.0554 2.4572 -15.8106 1.3425 -17.3978 2.3874 -15.8802
         1.0363 -12.832 3.4868 -12.6552 1.5616 -12.8482 3.4107 -12.618
         1.2722 -6.0885 3.6696 -7.2489 3.1926 -6.2539 3.4902 -7.1938
         -21.6368 -3.9407 -20.0697 -1.858 -22.8725 -3.0008 -19.9969 -1.9173
         -10.881 1.8578 -16.2534 4.9851 -10.9325 5.6716 -16.2289 5.0372
         -10.2945 5.9251 -10.8369 5.8362 -10.916 7.3663 -10.8181 5.8248
         -10.2127 6.6172 -7.7569 6.3592 -9.5409 7.5021 -7.7564 5.8689
         0.1809 -18.3909 -1.5284 -18.248 0.3702 -19.0828 -1.6867 -18.2219
         0.9409 -16.6205 3.6063 -17.3709 1.3317 -17.0823 3.2759 -17.4045
         1.1832 -11.7458 4.5451 -12.8079 1.5514 -11.5141 4.4245 -12.8021
         1.2884 -6.7147 5.5535 -6.2995 3.1779 -6.9786 5.4309 -6.2659
         -20.5609 -3.385 -18.6432 -2.9925 -21.4324 -2.4995 -18.6367 -3.0255
         -13.5124 0.313 -17.9269 5.7384 -14.1541 3.8794 -17.832 5.7364
         -10.85 5.1961 -11.0073 7.287 -10.2546 5.7663 -10.984 7.0992
         -9.1279 5.9001 -6.8506 7.4803 -9.2483 7.1392 -6.8757 7.3227
         -0.734 -19.2251 -2.7059 -19.5086 0.4781 -19.9873 -2.829 -19.5431
         0.4996 -15.7781 3.5483 -15.9278 2.3614 -16.2064 3.2084 -15.8791
         1.7599 -11.1979 6.9449 -12.5969 3.3873 -10.8691 6.7741 -12.6188
         2.6484 -7.4245 7.0337 -7.175 3.5694 -7.7464 6.864 -7.1947
         -20.0312 -2.5836 -24.0141 -1.9412 -20.6511 -1.6515 -23.86 -1.9158
         -14.2134 -0.2165 -12.2166 5.0941 -15.14 3.2235 -12.1001 5.0365
         -11.5306 4.0036 -8.3258 5.8148 -11.0609 4.4493 -8.1689 5.8244
         -8.596 4.8717 -7.5213 6.4545 -8.3831 5.4551 -7.8688 5.8675
         -1.4807 -19.8261 -0.9876 -20.9108 -0.7246 -20.6197 -0.98 -20.9423
         1.0375 -14.7725 -0.0406 -13.9992 3.4587 -15.0794 -0.0331 -13.9347
         2.5786 -11.4457 0.1484 -12.0027 4.2786 -11.1826 0.1658 -12.0476
         3.4245 -8.0033 1.8198 -8.6351 4.4771 -8.3435 1.8331 -8.6657
         -19.3164 -1.7678 -23.6812 0.0181 -19.5839 -0.6624 -23.5299 0.0534
         -14.6982 0.5513 -11.809 3.1039 -15.887 3.5102 -11.7008 3.0529
         -12.7332 3.1662 -9.2967 4.06 -12.5614 4.0645 -9.1918 3.9843
         -7.8473 4.3792 -8.3005 4.3667 -7.2368 4.5305 -8.5634 4.2486
         -2.3206 -20.1658 -0.5298 -22.0368 -1.9221 -20.944 -0.5282 -22.063
         2.4944 -13.9026 0.4238 -12.2447 4.9926 -13.9056 0.4247 -12.1958
         3.7144 -12.0151 0.6172 -10.9469 6.0105 -12.0873 0.6284 -10.9826
         4.3677 -8.3612 2.2805 -9.9241 6.3591 -8.6648 2.287 -9.9531
         -18.4885 -1.2412 -23.0156 1.5294 -18.2923 0.0784 -22.8733 1.484
         -15.147 2.3243 -11.0031 2.0253 -16.377 3.078 -10.9135 2.0486
         -13.8251 2.7074 -10.8758 2.4967 -14.2482 4.0328 -10.8927 2.424
         -6.9283 3.556 -9.5007 2.7412 -5.9005 4.3386 -9.5848 2.6314
         -3.2853 -20.2375 0.3891 -22.8554 -3.133 -20.95 0.3748 -22.8759
         4.1158 -14.1473 1.3524 -10.8913 6.2702 -13.9635 1.336 -10.9144
         4.7765 -11.6588 1.5588 -10.8861 6.8656 -12.0166 1.554 -10.8872
         4.986 -8.4651 3.1967 -9.5957 8.1388 -8.6701 3.1834 -9.5845
         -17.6661 -1.2824 -22.0253 0.3947 -17.0388 0.058 -21.8986 0.356
         -16.1727 2.6206 -12.9264 1.357 -16.9172 3.0706 -12.9717 1.2905
         -13.8597 3.8143 -10.4893 1.5762 -15.2851 4.0273 -10.4867 1.4749
         -5.8937 3.9125 -9.8288 3.2081 -4.5069 4.3313 -9.7681 3.1272
         -4.4035 -20.0463 1.7644 -23.3535 -4.4315 -20.6378 1.5657 -23.3678
         3.76 -15.0518 1.801 -11.4886 5.9262 -15.1422 1.7242 -11.5035
         4.8357 -10.8466 2.7758 -9.7728 8.4621 -11.0829 2.7249 -9.7908
         6.3295 -8.314 2.9725 -8.7105 8.6203 -8.3583 2.941 -8.6821
         -17.5965 -1.8061 -20.726 -0.2909 -18.1791 -0.7072 -20.6221 -0.322
         -16.7044 2.1404 -15.0913 0.6652 -16.4073 3.4938 -15.1108 0.6068
         -13.0961 2.9236 -11.1604 0.8734 -14.3379 4.2168 -11.1318 0.7813
         -5.6222 4.1281 -8.371 2.5261 -5.8159 4.5117 -8.3484 2.457
         -4.7752 -18.7821 -0.5164 -23.5247 -3.1675 -18.6654 -0.699 -23.5324
         2.3034 -16.9441 3.5792 -11.6962 6.2563 -17.8502 3.4961 -11.7015
         5.2443 -10.7119 4.4784 -9.2306 7.0105 -10.9326 4.1126 -9.1846
         5.7606 -7.3163 4.6003 -8.5277 7.731 -6.8888 4.5145 -8.5644
         -18.6699 -3.3062 -19.1661 -0.5265 -19.4969 -2.8433 -19.1007 -0.5478
         -15.9721 1.3598 -17.1237 0.4276 -15.8546 3.3228 -17.0951 0.3783
         -11.9766 5.5198 -11.4974 0.6327 -12.6951 6.7574 -11.4658 0.5499
         -6.7862 5.7011 -6.9736 2.2907 -7.1563 6.8479 -6.9915 2.2309
         -3.6298 -18.9185 -2.3964 -22.5576 -1.9098 -19.1209 -2.5304 -22.5493
         0.667 -16.5778 4.2959 -11.9536 5.0976 -17.1211 3.9376 -11.8915
         4.2816 -11.1495 5.8521 -10.5466 5.602 -11.4171 5.7318 -10.5283
         5.0034 -7.3365 6.8298 -9.509 6.4646 -7.0029 6.6876 -9.5799
         -19.5943 -3.1452 -19.1302 0.8023 -20.5876 -2.5425 -19.1008 0.8056
         -15.1067 2.0243 -17.1637 1.7509 -15.0624 4.0528 -17.095 1.7325
         -10.8938 4.795 -11.4892 1.9972 -11.225 5.7233 -11.466 1.9364
         -7.8055 5.6548 -6.9684 2.8798 -8.3172 7.0953 -6.9919 2.9159
         -2.462 -19.0966 -2.3933 -21.921 -0.6692 -19.554 -2.53 -21.9009
         -0.5966 -15.7605 4.2671 -12.9987 3.523 -15.9099 3.9372 -12.9707
         2.9994 -12.1724 5.9023 -10.4383 3.9184 -12.5811 5.7315 -10.4821
         3.8577 -7.3445 6.8589 -9.8044 4.5494 -7.2003 6.6871 -9.7685
         -20.2889 -2.8074 -20.694 1.5462 -21.3907 -1.9204 -20.6223 1.5653
         -14.1257 3.7208 -15.1309 1.6853 -14.0391 5.0498 -15.1105 1.707
         -10.3423 4.3296 -11.1418 2.6936 -10.4402 5.9015 -11.1326 2.6814
         -8.6066 5.5999 -8.3608 2.906 -9.2022 5.94 -8.3488 2.8617
         -1.3387 -19.1023 -0.4987 -20.9752 0.5593 -19.6998 -0.698 -20.942
         -0.9289 -15.0634 3.6069 -13.9159 2.4114 -14.9332 3.4955 -13.9378
         2.0413 -13.0616 4.4013 -12.0375 3.3212 -13.5933 4.112 -12.0443
         2.9067 -7.2527 4.6301 -8.7183 3.435 -7.2844 4.5141 -8.6648
         -20.7354 -2.6382 -22.0001 0.0758 -21.8833 -1.578 -21.8986 0.0512
         -13.0794 4.0938 -12.9403 3.0122 -12.8911 4.843 -12.9706 3.0534
         -10.5889 4.2785 -10.4859 3.9709 -10.5918 5.5872 -10.4882 3.9853
         -9.1506 5.5081 -9.8182 4.2766 -9.7561 6.1843 -9.7685 4.2491
         -0.422 -18.8731 1.7904 -19.7445 1.6041 -19.5387 1.5666 -19.6961
         -0.0487 -15.5723 1.8413 -14.838 1.7377 -15.8522 1.7239 -14.9373
         1.4345 -12.6518 2.8011 -13.6357 2.7356 -12.666 2.7245 -13.5834
         2.2858 -7.0264 2.9991 -7.3537 2.9567 -7.1885 2.9407 -7.2836
         -20.9275 -2.8839 -22.9976 -1.5323 -22.0554 -1.8886 -22.8734 -1.5925
         -12.253 1.9832 -10.9942 4.7587 -12.1835 5.0657 -10.9137 4.8221
         -11.1843 4.5253 -10.7854 5.511 -11.0036 5.7678 -10.8897 5.5552
         -9.4193 5.1412 -9.5984 6.1046 -9.9516 5.9081 -9.588 6.0951
         0.0145 -18.4073 0.4077 -18.3024 1.5023 -19.0888 0.3746 -18.2351
         1.1698 -16.5947 1.3704 -16.3011 2.0502 -17.0751 1.3358 -16.4012
         1.6542 -11.7587 1.578 -14.3189 2.4689 -11.516 1.5537 -14.2886
         2.0416 -6.7306 3.2139 -5.9284 2.7138 -6.9812 3.1833 -5.8543
         -20.8663 -3.363 -23.6705 -3.0791 -21.904 -2.4971 -23.5299 -3.1509
         -12.6727 0.3181 -11.8034 6.236 -13.028 3.8821 -11.701 6.2666
         -10.9475 5.1813 -8.9456 6.8397 -10.4128 5.7576 -9.1855 6.9367
         -9.4091 5.8762 -8.6643 7.9629 -9.7797 7.1309 -8.5696 7.9609
         -0.1598 -17.9288 -0.5186 -18.1857 1.533 -18.6128 -0.5283 -18.236
         1.1444 -17.2865 0.4348 -16.4609 1.7109 -17.8237 0.4246 -16.3994
         1.2996 -11.3533 0.6286 -14.2724 2.7162 -11.0353 0.6283 -14.2896
         2.1806 -6.5746 2.2912 -5.8305 2.9298 -6.8632 2.2868 -5.8558
         -19.599 -3.6149 -24.0105 -3.1752 -20.0274 -2.7958 -23.86 -3.1491
         -15.1461 -0.332 -12.2146 6.2773 -16.2633 3.1544 -12.1002 6.2655
         -11.1377 6.1886 -8.0945 7.024 -10.7317 6.811 -8.163 6.9349
         -8.249 6.1918 -7.7619 7.9511 -7.7967 6.9014 -7.8739 7.9604
         -2.08 -23.4487 -0.9839 -19.6463 -1.7132 -24.0297 -0.98 -19.6967
         0.3304 -11.8679 -0.037 -15.0073 3.4354 -12.3042 -0.0331 -14.9341
         3.3634 -8.229 0.1522 -13.5552 4.3852 -7.5309 0.1657 -13.5862
         4.2335 -7.4606 1.8235 -7.2307 5.3934 -7.5092 1.833 -7.2853
         -19.1464 -1.1733 -23.5166 -1.6383 -19.151 -1.206 -23.3653 -1.59
         -15.5731 -0.2428 -11.6075 4.8961 -17.1059 -0.2626 -11.5028 4.8211
         -11.6885 -0.0229 -9.935 5.6294 -11.3856 -0.066 -9.7965 5.5539
         -7.8399 1.6389 -8.3839 6.239 -7.0436 1.6053 -8.682 6.0939
         -2.5439 -23.2958 -0.3031 -20.9011 -2.5414 -23.8641 -0.3028 -20.9425
         0.9111 -11.6882 0.654 -13.95 4.1122 -12.1032 0.6528 -13.9354
         4.2889 -8.6104 0.8498 -12.0647 5.6694 -8.1645 0.8595 -12.0467
         5.0733 -8.0182 2.5083 -8.6145 7.0435 -7.8682 2.5124 -8.6661
         -18.6815 -0.9665 -23.0224 0.0252 -18.2877 -0.9801 -22.8734 0.0534
         -15.5964 -0.0339 -11.0174 3.1172 -17.3934 -0.0332 -10.9135 3.0526
         -12.7747 0.1884 -11.0063 4.0693 -12.7348 0.1654 -10.8938 3.9841
         -7.3236 1.8457 -9.3357 4.3803 -6.3177 1.8329 -9.5835 4.2487
         -2.9572 -22.8386 0.3816 -21.8692 -3.0211 -23.3691 0.3748 -21.9012
         2.4536 -11.1548 1.3456 -12.9461 5.9161 -11.5054 1.3359 -12.969
         5.1194 -9.9881 1.5514 -10.5325 7.261 -9.7925 1.5539 -10.4841
         5.3378 -8.9244 3.1898 -9.7279 7.417 -8.6825 3.1833 -9.7695
         -18.2843 -0.3472 -22.2008 1.5712 -18.247 -0.3024 -22.0606 1.5666
         -15.4113 0.5903 -12.2541 1.7601 -16.4246 0.6532 -12.2022 1.7064
         -13.7631 0.8232 -10.8576 2.7534 -14.2462 0.8598 -10.9811 2.6809
         -6.7495 2.4599 -10.0218 2.9798 -5.8833 2.5127 -9.9518 2.8607
         -3.437 -22.0837 1.5256 -22.5269 -3.1394 -22.5502 1.5018 -22.5494
         4.1264 -11.9805 2.2609 -11.8495 6.2714 -11.8914 2.0463 -11.8901
         4.6947 -10.288 2.4967 -10.5001 6.935 -10.5297 2.4675 -10.5288
         4.9161 -9.7109 2.726 -9.6428 8.1316 -9.5843 2.7113 -9.5816
         -17.8946 0.6836 -21.0648 0.8468 -18.2162 0.8247 -20.9401 0.8053
         -16.0863 1.546 -13.9044 1.8046 -16.3565 1.7769 -13.9383 1.732
         -13.4584 1.8799 -12.0605 2.0415 -14.3528 2.0165 -12.0477 1.936
         -6.1709 1.9847 -8.7008 2.9333 -5.8467 2.9295 -8.6644 2.9164
         -4.0805 -21.0458 0.2319 -22.8632 -3.1374 -21.4184 0.0512 -22.876
         2.5185 -14.0678 3.1196 -10.9002 6.244 -14.1055 3.069 -10.9144
         5.0191 -10.352 4.0862 -10.8702 6.9538 -10.3349 4.0243 -10.8874
         5.9335 -9.1409 4.3722 -9.6144 7.7368 -9.2257 4.3296 -9.5843
         -17.4008 -0.1888 -19.6409 0.3841 -18.1591 0.5227 -19.5413 0.3558
         -17.2015 2.0942 -15.8572 1.347 -17.387 2.3894 -15.8808 1.2904
         -12.3875 2.428 -12.6685 1.5658 -12.8917 3.4128 -12.6205 1.4747
         -5.8965 3.1121 -7.1874 3.199 -6.2323 3.5005 -7.1935 3.1278
         -4.5236 -19.7711 -1.7721 -22.8734 -3.0078 -19.9988 -1.9177 -22.8759
         0.6639 -16.0212 5.1411 -10.9658 5.5078 -16.2296 5.0498 -10.9143
         5.6771 -10.766 6.011 -10.9158 7.3832 -10.8187 5.8248 -10.8883
         6.1733 -7.881 6.2242 -9.5058 7.5187 -7.7561 5.8958 -9.5832
         -18.2192 -2.0822 -18.2272 0.3693 -19.0483 -1.6826 -18.2204 0.356
         -16.735 2.1754 -17.481 1.3305 -17.0614 3.2849 -17.4051 1.2905
         -11.3454 3.9687 -12.8508 1.5511 -11.5671 4.4276 -12.8047 1.475
         -6.4788 5.2624 -6.2268 3.1768 -6.9529 5.4343 -6.2659 3.1269
         -4.0312 -18.7439 -2.8976 -21.4394 -2.5 -18.6365 -3.0256 -21.418
         -0.7655 -17.2817 6.0962 -14.1613 3.7261 -17.8348 5.7339 -14.1054
         4.9139 -10.9085 7.2795 -10.2409 5.7868 -10.9847 7.1001 -10.3294
         5.6441 -7.167 7.499 -9.2632 7.1627 -6.8746 7.3306 -9.2252
         -19.035 -3.0223 -19.5721 0.4766 -19.9593 -2.826 -19.5416 0.5162
         -15.8861 2.1 -15.9301 2.3505 -16.1762 3.2175 -15.8799 2.3714
         -10.8451 6.0792 -12.6595 3.3786 -10.9274 6.7785 -12.6215 3.3702
         -7.2631 6.2329 -7.1741 3.5678 -7.7246 6.8685 -7.1944 3.472
         -3.1406 -23.2977 -1.741 -20.6624 -1.6447 -23.8642 -1.916 -20.6242
         -1.2435 -11.6875 5.2143 -15.1371 3.0778 -12.1033 5.0494 -15.1108
         3.801 -8.8344 6.0561 -11.054 4.4626 -8.1731 5.8248 -11.127
         4.6329 -7.7755 6.1267 -8.4099 5.4671 -7.8595 5.8949 -8.3488
         -19.6185 -0.9684 -21.0077 -0.7039 -20.5974 -0.9803 -20.9402 -0.6995
         -14.8829 -0.0354 -13.9167 3.4393 -15.0427 -0.0333 -13.9363 3.482
         -11.1608 0.1872 -12.0992 4.3241 -11.24 0.1653 -12.0506 4.1135
         -7.9113 1.8446 -8.6823 4.4597 -8.3287 1.8328 -8.6654 4.474
         -2.2128 -23.0028 0.3067 -19.6012 -0.6543 -23.534 0.053 -19.5429
         -0.4794 -11.3403 3.1816 -15.8714 3.5078 -11.7038 3.0684 -15.8812
         3.0105 -9.6314 4.1375 -12.5612 3.927 -9.1942 4.0236 -12.6166
         3.8196 -8.6023 4.4342 -7.2766 4.5352 -8.5568 4.329 -7.1939
         -19.9588 -0.5686 -22.1578 -1.8806 -20.928 -0.5283 -22.0607 -1.9179
         -13.8622 0.3693 -12.049 4.9646 -13.8695 0.4246 -12.1975 5.0372
         -11.9408 0.5969 -11.1059 6.0535 -12.1327 0.6281 -10.9865 5.828
         -8.3428 2.2459 -10.0041 6.3403 -8.6587 2.2868 -9.9526 5.8683
         -1.5157 -22.4068 1.5716 -18.3217 0.0768 -22.8773 1.5014 -18.2353
         1.2159 -11.0318 2.3821 -16.3383 3.0776 -10.916 2.048 -16.4019
         2.5731 -10.6473 2.5384 -14.2582 4.03 -10.8914 2.4671 -14.2877
         3.3536 -9.7445 2.7725 -5.9532 4.3403 -9.5834 2.7107 -5.8544
         -20.0566 0.2445 -22.9936 -3.0764 -20.9411 0.3751 -22.8735 -3.1512
         -13.6405 1.1843 -10.9927 6.2506 -13.9816 1.3363 -10.9139 6.2666
         -12.0826 1.4309 -10.6649 6.8285 -11.9949 1.5541 -10.8888 6.9364
         -8.536 2.624 -9.7337 8.19 -8.6741 3.1836 -9.5888 7.9636
         -1.3184 -21.5195 0.412 -17.1008 0.0403 -21.9021 0.3745 -16.9814
         2.4607 -12.9411 1.3742 -16.8548 3.061 -12.9699 1.3356 -16.9771
         2.5731 -10.585 1.5822 -15.2845 4.018 -10.4873 1.5537 -15.2874
         3.4619 -9.6399 3.2174 -4.5706 4.3236 -9.7698 3.1829 -4.4692
         -19.9255 0.7234 -23.5023 -4.3691 -20.6371 1.5723 -23.3654 -4.4649
         -14.4215 1.4656 -11.6003 8.3901 -15.1615 1.7251 -11.503 5.9293
         -11.3823 2.6848 -9.5658 8.6791 -11.0479 2.7259 -9.7924 8.5597
         -8.4959 3.3131 -8.7939 8.9974 -8.3737 2.9418 -8.6859 8.5638
         -1.6008 -20.3649 -0.2878 -18.1512 -0.729 -20.625 -0.3029 -18.2368
         0.9796 -14.9333 0.6685 -16.4541 3.4752 -15.1099 0.6526 -16.3991
         2.7437 -11.1322 0.8652 -14.3155 4.2549 -11.1329 0.8595 -14.2891
         3.6018 -8.4041 2.522 -5.7603 4.4933 -8.349 2.5122 -5.8575
         -19.0399 -1.2272 -23.6773 -3.2304 -18.707 -0.6939 -23.53 -3.1475
         -15.9333 2.954 -11.8071 6.2682 -17.8264 3.498 -11.7011 6.2656
         -11.2381 3.1609 -9.2114 7.0355 -10.8955 4.1223 -9.1867 6.9349
         -7.8218 4.3712 -8.391 7.6624 -6.9314 4.5167 -8.5684 7.9664
         -2.6768 -19.0547 -0.5258 -19.4803 -2.8416 -19.1025 -0.5284 -19.5437
         0.2659 -16.7044 0.4278 -15.8838 3.3756 -17.0957 0.4245 -15.8796
         5.0803 -11.3931 0.6213 -12.6755 6.7073 -11.4669 0.6281 -12.6178
         5.4508 -7.2567 2.2843 -7.1099 6.8067 -6.9908 2.2868 -7.1963
         -19.079 -2.7761 -22.6935 -1.9645 -19.1355 -2.5269 -22.547 -1.9148
         -15.6745 2.826 -12.0262 5.1161 -17.1243 3.9472 -11.8943 5.0361
         -11.6377 5.0904 -10.6118 5.5392 -11.3807 5.7343 -10.5279 5.8299
         -7.7893 6.0166 -9.3243 6.4881 -7.0379 6.6911 -9.583 5.8677
         -2.6039 -18.971 0.8356 -20.5769 -2.546 -19.1012 0.8236 -20.6247
         0.8847 -16.8058 1.7948 -15.0837 4.1036 -17.0972 1.7756 -15.1098
         4.5052 -11.3608 2.0187 -11.2037 5.6925 -11.4666 2.0155 -11.1281
         5.2687 -7.2179 3.1041 -8.281 7.0649 -6.9907 2.921 -8.3505
         -19.1183 -2.8104 -22.0432 -0.7113 -19.5535 -2.5268 -21.8987 -0.6973
         -15.0369 2.7816 -13.0941 3.5328 -15.9264 3.9461 -12.9729 3.4811
         -12.5615 5.2444 -10.2681 3.8657 -12.549 5.7361 -10.4851 4.1156
         -7.6722 6.1366 -9.8352 4.5648 -7.2201 6.6921 -9.7682 4.4734
         -2.4992 -20.2769 1.7116 -21.3841 -1.9295 -20.6239 1.565 -21.4183
         2.4874 -15.069 1.7443 -14.0588 5.0313 -15.1114 1.724 -14.1046
         4.0602 -11.0437 2.7532 -10.4155 5.9831 -11.1327 2.7246 -10.3304
         4.9681 -8.3384 2.9526 -9.176 6.3975 -8.3487 2.9408 -9.2264
         -19.0366 -1.2725 -21.0759 0.5285 -19.689 -0.693 -20.9402 0.5174
         -14.7006 2.8622 -14.0336 2.4234 -14.9252 3.4988 -13.9393 2.371
         -13.1226 3.2379 -11.9016 3.276 -13.6005 4.1195 -12.0466 3.3682
         -7.4321 4.4873 -8.7026 3.4444 -7.2888 4.5176 -8.6645 3.4693
         -2.6382 -21.4492 0.1767 -21.8803 -1.5823 -21.9011 0.0508 -21.9013
         2.8648 -13.0622 3.1067 -12.9193 4.8357 -12.9702 3.0689 -12.9688
         4.0695 -10.4755 4.0758 -10.5585 5.5745 -10.4878 4.0242 -10.4842
         4.6004 -9.579 4.3595 -9.7398 6.1787 -9.7695 4.3295 -9.7696
         -18.7636 0.7116 -19.8132 1.5791 -19.519 1.5738 -19.6946 1.5665
         -15.5666 1.557 -14.9561 1.7451 -15.8181 1.7259 -14.9381 1.7063
         -12.3582 2.7718 -13.5623 2.7388 -12.7136 2.7267 -13.5852 2.6807
         -7.0496 3.0893 -7.2806 2.9659 -7.1791 2.9426 -7.2833 2.8605
         -3.1708 -22.3563 -1.4778 -22.0556 -1.8843 -22.8766 -1.5927 -22.0632
         0.8487 -11.0306 4.9002 -12.236 5.0639 -10.9156 4.8353 -12.1963
         4.3219 -10.6003 5.6864 -10.9468 5.6106 -10.887 5.5792 -10.982
         4.8329 -9.7384 6.235 -9.9452 6.4349 -9.5881 6.1764 -9.9532
         -18.2724 0.3099 -18.3175 1.5042 -19.0598 0.3757 -18.2341 1.4837
         -16.6238 1.2438 -16.3752 2.0194 -17.0467 1.3369 -16.4011 2.0483
         -11.4247 1.4954 -14.3628 2.4731 -11.572 1.5547 -14.2903 2.4243
         -6.5799 2.6643 -5.806 2.7163 -6.9605 3.1844 -5.8541 2.6313
         -3.8802 -22.9725 -3.0422 -21.9075 -2.4909 -23.5336 -3.151 -21.9011
         -0.7295 -11.3111 6.4274 -13.0488 3.7285 -11.7036 6.2678 -12.9695
         4.9378 -9.2295 7.0502 -10.3871 5.7702 -9.1784 6.9454 -10.4834
         5.561 -9.0034 8.3371 -9.7836 7.1462 -8.5727 7.9586 -9.769
         -17.6436 -0.5291 -18.2183 1.5109 -18.5702 -0.5279 -18.2347 1.5662
         -17.4564 0.4065 -16.4536 1.7069 -17.808 0.4249 -16.3997 1.7065
         -11.0137 0.6358 -14.3867 2.7161 -11.0942 0.6285 -14.2914 2.6808
         -6.3063 2.2822 -5.8079 2.9264 -6.8375 2.2871 -5.8557 2.8614
         -4.2616 -23.2876 -2.9985 -20.0409 -2.7905 -23.864 -3.1492 -19.999
         -1.3539 -11.6776 6.4007 -16.2666 3.0041 -12.1032 6.2669 -16.229
         5.8106 -8.622 7.187 -10.7136 6.8339 -8.1545 6.9441 -10.8126
         5.9427 -7.9942 8.2869 -7.8244 6.9203 -7.8776 7.9583 -7.7569
         -22.9775 -0.9552 -19.7197 -1.7026 -24.0363 -0.9801 -19.6949 -1.6875
         -11.5295 -0.0229 -14.879 3.4939 -12.3093 -0.0332 -14.9344 3.28
         -8.182 0.1997 -13.744 4.3703 -7.5646 0.1654 -13.589 4.4099
         -7.9557 1.8569 -7.2674 5.3787 -7.4666 1.8329 -7.285 5.3917
         -1.1601 -22.8464 -1.3689 -19.1744 -1.2061 -23.3695 -1.5903 -19.1024
         -0.2663 -11.151 5.0089 -17.0971 -0.2623 -11.5057 4.8343 -17.0956
         -0.0626 -10.3534 5.7292 -11.3698 -0.066 -9.8007 5.5779 -11.4606
         1.485 -8.5067 6.3404 -7.0844 1.6053 -8.6734 6.1754 -6.9919
         -22.8346 -0.3539 -21.0011 -2.513 -23.8707 -0.303 -20.9403 -2.5311
         -11.3635 0.5853 -13.8276 4.1709 -12.1081 0.6527 -13.9359 3.9419
         -8.7922 0.8162 -12.2065 5.6463 -8.1703 0.8592 -12.0504 5.7184
         -8.2291 2.458 -8.682 7.0183 -7.8543 2.5122 -8.6658 6.6577
         -0.9674 -22.4125 0.341 -18.3273 -0.9801 -22.8776 0.0529 -18.2225
         -0.0718 -11.2648 3.1889 -17.3674 -0.0328 -10.9162 3.0682 -17.405
         0.1343 -10.6434 4.1432 -12.7198 0.1655 -10.8974 4.0234 -12.8008
         1.8025 -9.4755 4.4414 -6.3578 1.8331 -9.5764 4.3287 -6.2661
         -22.4079 0.239 -21.9896 -2.9911 -23.3754 0.3747 -21.8989 -3.0262
         -10.8707 1.1825 -12.8284 5.9811 -11.5097 1.3359 -12.9699 5.7393
         -10.312 1.4253 -10.6311 7.2302 -9.7939 1.5537 -10.4883 7.321
         -8.8532 2.565 -9.8156 7.4033 -8.6746 3.1833 -9.7691 7.4517
         -0.3891 -21.6835 1.8018 -18.257 -0.3019 -22.0646 1.5664 -18.2367
         0.5092 -12.3307 1.8955 -16.4336 0.6541 -12.203 1.7234 -16.3996
         0.7257 -10.9299 2.815 -14.2251 0.8603 -10.9783 2.7245 -14.2886
         2.0114 -9.8043 3.0104 -5.888 2.5134 -9.954 2.9402 -5.856
         -21.704 1.0794 -22.6613 -3.1459 -22.5561 1.5021 -22.547 -3.1494
         -12.2102 1.2429 -11.7137 6.2834 -11.8929 2.0522 -11.8913 6.2659
         -10.0697 2.1879 -10.5982 6.9358 -10.5331 2.4679 -10.5284 6.9343
         -9.555 2.458 -9.724 8.1945 -9.5782 2.7116 -9.5854 7.9638
         0.5722 -20.6756 0.8703 -18.204 0.826 -20.9438 0.8233 -18.2373
         0.9979 -13.8188 1.8262 -16.3904 1.7789 -13.9365 1.7751 -16.3986
         1.5268 -12.0184 2.0535 -14.3272 2.0179 -12.0488 2.0151 -14.2891
         1.7004 -8.6803 3.2571 -5.7937 2.9392 -8.6657 2.921 -5.8577
         -20.7388 -0.6176 -23.005 -3.2032 -21.4237 0.0564 -22.8736 -3.1475
         -14.1635 2.6479 -10.9985 6.2591 -14.108 3.0699 -10.914 6.2655
         -10.147 3.9115 -10.7435 6.9597 -10.3298 4.0254 -10.8889 6.9339
         -9.0087 4.9252 -9.6511 7.6715 -9.2277 4.3304 -9.5886 7.9662
         -0.8559 -19.4337 0.4001 -18.0943 0.5296 -19.5443 0.3744 -18.2253
         1.5292 -15.5809 1.3627 -17.4597 2.3918 -15.8795 1.3355 -17.4019
         1.9334 -12.5262 1.5701 -12.8773 3.4159 -12.6222 1.5535 -12.8013
         3.0948 -7.4111 3.2063 -6.1497 3.513 -7.1935 3.1829 -6.2697
         -19.5639 -2.2868 -23.0157 -3.101 -20.0031 -1.9134 -22.8735 -3.0223
         -15.974 4.4079 -11.0033 5.4391 -16.2332 5.0516 -10.9138 5.7423
         -10.5315 4.7013 -10.9955 7.4155 -10.8141 5.8356 -10.8909 7.3198
         -7.839 5.2548 -9.3609 7.5307 -7.7562 5.8983 -9.5866 7.4464
         -2.596 -18.3488 0.3883 -19.0227 -1.6791 -18.2202 0.3745 -19.1036
         1.3188 -16.8901 1.3521 -17.091 3.307 -17.4075 1.3356 -17.0942
         3.6286 -12.648 1.5584 -11.553 4.4318 -12.806 1.5537 -11.4612
         4.866 -6.6043 3.1965 -6.8807 5.4389 -6.2644 3.1831 -6.9952
         -18.6599 -3.2281 -21.5554 -2.5789 -18.6381 -3.0225 -21.416 -2.5276
         -17.0962 4.5284 -14.2413 3.658 -17.8412 5.7437 -14.1071 3.9447
         -10.6639 6.4292 -10.1009 5.8109 -10.9802 7.336 -10.333 5.7173
         -7.1699 6.4448 -9.2762 7.1898 -6.8731 7.4666 -9.225 6.6545
         -3.4682 -19.2484 0.6077 -19.9468 -2.8246 -19.5418 0.5158 -19.9996
         1.2392 -15.8233 2.4195 -16.193 3.2401 -15.882 2.3871 -16.2283
         5.5901 -12.4416 3.6471 -10.9116 6.7842 -12.6222 3.4104 -10.8132
         5.7188 -7.2928 4.1612 -7.6723 6.8734 -7.1932 3.4897 -7.7592
         -22.8356 -2.3607 -20.7591 -1.7027 -23.8713 -1.912 -20.6224 -1.6851
         -11.3659 4.4984 -15.2379 3.0152 -12.1087 5.054 -15.1119 3.2827
         -8.6149 4.7401 -10.9312 4.4797 -8.202 5.8319 -11.1301 4.4086
         -8.4179 5.4605 -8.3748 5.4838 -7.822 5.9 -8.3486 5.3902
         -0.9691 -20.5165 -0.6467 -20.5912 -0.9808 -20.9416 -0.6998 -20.6248
         -0.0734 -14.0237 3.544 -15.0584 -0.0335 -13.9366 3.4958 -15.1095
         0.1328 -11.9254 4.5678 -11.2216 0.1648 -12.0513 4.115 -11.1282
         1.7838 -8.5518 4.7439 -8.2944 1.8323 -8.6652 4.5142 -8.3506
         -22.5635 -0.6364 -19.6665 -0.6981 -23.542 0.0593 -19.5415 -0.6973
         -11.0528 2.8581 -15.9935 3.5159 -11.7099 3.0721 -15.8817 3.4809
         -9.5452 4.1158 -12.4759 3.8812 -9.2143 4.0271 -12.6191 4.1158
         -9.0232 4.5481 -7.1803 4.5465 -8.5298 4.3327 -7.1936 4.4731
         -0.601 -21.5628 -1.8798 -20.9265 -0.53 -22.063 -1.9183 -20.9428
         0.298 -12.2255 5.11 -13.9019 0.4234 -12.1941 5.0501 -13.9351
         0.5101 -11.107 5.998 -12.0962 0.6265 -10.9893 5.8286 -12.0468
         2.1673 -9.6966 6.4402 -8.6406 2.2854 -9.9536 5.8955 -8.6663
         -22.0108 1.1639 -18.3297 0.0398 -22.8859 1.5037 -18.2343 0.0533
         -11.1462 1.3994 -16.4874 3.0833 -10.9223 2.0568 -16.4019 3.0522
         -10.4202 2.6126 -14.2271 4.0369 -10.8914 2.4693 -14.2894 3.9838
         -9.8015 3.7314 -5.7864 4.3471 -9.5787 2.713 -5.8542 4.2485
         0.1517 -22.3318 -3.1115 -20.944 0.3726 -22.8765 -3.1513 -20.9427
         1.0512 -10.7719 6.4644 -14.0129 1.3345 -10.9157 6.268 -13.9354
         1.2809 -10.5674 7.0317 -11.9574 1.5518 -10.8825 6.9453 -12.0465
         1.7797 -10.0279 8.5462 -8.6712 3.1822 -9.5929 7.9623 -8.6659
         -21.188 0.343 -16.9885 0.0123 -21.9113 0.3757 -16.9806 0.0529
         -13.0382 1.2826 -16.9739 3.06 -12.9645 1.3369 -16.9761 3.0524
         -10.4704 1.5281 -15.4184 4.0204 -10.4889 1.5547 -15.2897 3.984
         -9.4963 2.5785 -4.4259 4.3231 -9.7757 3.1843 -4.4691 4.2485
         0.0478 -22.8061 -4.3326 -20.6448 1.5819 -23.369 -4.4648 -20.6245
         1.2819 -11.112 6.0646 -15.1746 1.7223 -11.5055 5.9294 -15.11
         2.2055 -9.8426 8.7757 -11.0251 2.7238 -9.7849 8.5575 -11.1277
         2.3662 -9.0719 8.8403 -8.386 2.9395 -8.6902 8.5633 -8.3493
         -20.1211 -0.3016 -18.2116 -0.7392 -20.6346 -0.3024 -18.2353 -0.6986
         -14.9269 0.6356 -16.3806 3.4692 -15.105 0.6532 -16.3983 3.4813
         -10.9357 0.8681 -14.4739 4.3188 -11.1338 0.8598 -14.2921 4.1134
         -8.3565 2.5088 -5.822 4.4878 -8.3536 2.5126 -5.8572 4.4772
         -1.752 -22.9783 -2.9538 -18.7465 -0.6857 -23.534 -3.1474 -18.639
         2.108 -11.3087 6.3769 -17.8044 3.4962 -11.7039 6.267 -17.832
         2.8169 -9.7331 7.1981 -10.8759 4.1518 -9.185 6.944 -10.9785
         3.7082 -8.4681 8.1057 -6.9731 4.5146 -8.566 7.9623 -6.8761
         -18.9293 -0.5338 -19.557 -2.8173 -19.1113 -0.5283 -19.5418 -2.8299
         -16.5703 0.4024 -15.8441 3.4389 -17.0936 0.4246 -15.8793 3.2123
         -11.1503 0.6308 -12.7806 6.6785 -11.4648 0.6282 -12.6217 6.7628
         -7.277 2.2786 -7.1815 6.7874 -6.992 2.2868 -7.1961 6.8511
         -3.2102 -22.1006 -1.6728 -19.1543 -2.5224 -22.5511 -1.9151 -19.1028
         1.9342 -12.2928 5.2307 -17.1236 3.9729 -11.8986 5.0484 -17.0951
         4.6483 -10.2754 5.9534 -11.3603 5.7311 -10.5305 5.8275 -11.4608
         5.5153 -9.3152 6.061 -7.0688 6.6909 -9.5751 6.4119 -6.9923
         -18.8248 0.6677 -20.6744 -2.5332 -19.099 0.8235 -20.6226 -2.5306
         -16.709 1.5278 -15.0388 4.1705 -17.1082 1.7757 -15.1099 3.9419
         -11.0943 1.8282 -11.2858 5.6762 -11.4608 2.0152 -11.1328 5.7178
         -7.2018 1.8646 -8.3616 7.0496 -6.9872 2.9256 -8.3505 6.6554
         -3.2999 -21.5336 -0.4174 -19.5617 -2.5265 -21.9029 -0.6976 -19.5437
         1.8966 -13.2007 3.6234 -15.9403 3.9667 -12.9767 3.4954 -15.8796
         4.842 -10.2558 4.2489 -12.5241 5.749 -10.4783 4.1172 -12.6176
         5.7021 -9.6238 4.6491 -7.2286 6.7053 -9.7704 4.5137 -7.195
         -20.0085 0.5899 -21.4985 -1.9409 -20.6238 1.5687 -21.4161 -1.9163
         -15.1428 1.4501 -13.9911 5.0286 -15.122 1.724 -14.1052 5.0359
         -10.7696 2.4171 -10.4913 6.0582 -11.1234 2.725 -10.3353 5.8275
         -8.2357 2.6712 -9.2615 6.3914 -8.3459 2.9407 -9.2262 5.8678
         -1.9075 -20.683 0.8239 -19.6899 -0.6901 -20.9443 0.5171 -19.6972
         2.0175 -14.008 2.4781 -14.9642 3.5086 -13.9419 2.3866 -14.934
         2.9544 -11.7954 3.5103 -13.5566 4.1346 -12.0422 3.41 -13.5857
         4.1614 -8.6645 3.6257 -7.2724 4.5274 -8.6661 3.4911 -7.2855
         -21.099 -0.7035 -22.0054 -1.6185 -21.9029 0.0544 -21.899 -1.5902
         -13.2993 2.6372 -12.8217 4.8384 -12.9791 3.0692 -12.9697 4.8205
         -10.2033 3.9019 -10.6318 5.5865 -10.4764 4.0249 -10.4885 5.5536
         -9.3886 5.2093 -9.8227 6.1826 -9.7687 4.3298 -9.7692 6.0935
         -0.0946 -19.577 1.7853 -19.5123 1.5822 -19.6983 1.5663 -19.5441
         1.4088 -14.8268 1.8704 -15.8367 1.7328 -14.9377 1.7233 -15.8791
         2.2737 -13.3009 2.7985 -12.6955 2.7342 -13.5849 2.7243 -12.6179
         2.5705 -7.4612 2.9942 -7.136 2.9498 -7.2838 2.9401 -7.1965
         -21.9474 -2.0928 -22.1849 -1.9416 -22.8799 -1.5891 -22.061 -1.9148
         -11.4363 4.2012 -12.1946 5.0728 -10.9168 4.8361 -12.198 5.0353
         -10.3393 4.9582 -10.9345 5.5575 -10.8887 5.5813 -10.9859 5.8303
         -9.4687 5.5489 -10.0162 6.4468 -9.581 6.1772 -9.9528 5.8674
         0.2395 -18.3176 1.5423 -19.0397 0.3802 -18.2369 1.5011 -19.1038
         1.1229 -16.0336 2.3064 -17.0703 1.3417 -16.3996 2.0477 -17.094
         1.3675 -14.0158 2.5117 -11.5591 1.5593 -14.2925 2.4668 -11.4612
         1.7086 -6.2313 2.7428 -6.8897 3.1886 -5.8531 2.7106 -6.9953
         -22.5256 -3.3868 -22.0339 -2.5696 -23.5384 -3.1478 -21.8989 -2.5276
         -11.0032 5.7012 -13.1009 3.6639 -11.7065 6.2719 -12.9717 3.9448
         -9.7847 6.016 -10.2658 5.7872 -9.1564 6.9468 -10.4867 5.7171
         -8.7403 6.7396 -9.8318 7.1669 -8.588 7.9683 -9.7687 6.6541
         -0.5485 -18.0695 1.7267 -18.5061 -0.5254 -18.2333 1.5659 -18.6418
         0.3485 -16.294 1.7537 -17.8737 0.4278 -16.4022 1.7235 -17.8292
         0.5616 -14.0013 2.762 -11.0813 0.631 -14.2938 2.7242 -10.979
         2.1366 -6.1048 2.9623 -6.7519 2.2897 -5.8535 2.9404 -6.8798
         -22.823 -3.4575 -20.1241 -2.8818 -23.87 -3.1459 -19.9974 -2.8261
         -11.3492 5.7353 -16.3748 2.9359 -12.1075 6.2708 -16.2297 3.2153
         -8.6716 6.4453 -10.5995 6.869 -8.1161 6.9517 -10.8164 6.7611
         -8.3461 6.721 -7.7564 6.9334 -7.9083 7.9675 -7.7567 6.8505
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 1.39 1.68 2
         "" "electrons" 337.65 402.35 10
         "" "update_pot" 1.39 1.48 8
         "" "forces" 1.48 2.29 10
         "" "stress" 3.9 6.39 10
         "init_run" "wfcinit" 1.11 1.35 2
         "init_run" "potinit" 0.08 0.11 2
         "electrons" "c_bands" 301.98 357.53 136
         "electrons" "sum_band" 34.17 42.89 136
         "electrons" "v_of_rho" 0.43 0.72 145
         "electrons" "newd" 0.79 0.91 145
         "electrons" "mix_rho" 0.22 0.25 136
         "c_bands" "init_us_2" 3.21 5.8 21756
         "c_bands" "cegterg" 296.73 350.66 10064
         "sum_band" "sum_band:bec" 0.13 0.44 10064
         "sum_band" "addusdens" 0.87 0.99 136
         "*egterg" "h_psi" 243.45 300.37 72271
         "*egterg" "s_psi" 5.51 4.3 72271
         "*egterg" "g_psi" 1.03 1.17 62059
         "*egterg" "cdiaghg" 12.01 12.56 71457
         "h_psi" "h_psi:pot" 242.38 299.04 72271
         "h_psi" "h_psi:calbec" 7.21 7.52 72271
         "h_psi" "vloc_psi" 229.07 285.54 72271
         "h_psi" "add_vuspsi" 5.92 5.66 72271
         "General routines" "calbec" 9.91 10.24 86035
         "General routines" "fft" 0.68 0.78 931
         "General routines" "fftw" 242.89 296.07 716666
         "General routines" "davcio" 0.0 0.01 148
         "Parallel routines" "fft_scatter" 69.24 72.55 717597
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse relax CO output" begin
    url = "https://raw.githubusercontent.com/maxhutch/deprecated-quantum-espresso/master/PW/tests/relax-damped.ref"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 2,
        alat = 12.0,
        omega = 1728.0,
        nat = 2,
        ntyp = 2,
        nelec = 10.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 5,
        ecutwfc = 24.0,
        ecutrho = 144.0,
        ecutfock = nothing,
        conv_thr = 1.0e-6,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)",
        nstep = 50,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Sum" 1649 1101 277
         "gvecs" "Sum" 50541 27609 3407
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [0.0 0.0 0.0 2.0], cryst = nothing)

    @test all(isempty, parse_stress(str))

    @test parse_cell_parameters(str) == [[
        12.0 0.0 0.0
        0.0 12.0 0.0
        0.0 0.0 12.0
    ]]

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.161309101, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.05503841, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.111613831, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.178918345, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.166035881, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.140753228, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.115110591, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.127180324, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.144570629, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.142564627, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.139519983, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr",QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.139767533, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("C", [2.139767533, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 24.0 DavidsonDiagonalization() 0.01 2.0 0.7 1.1 -43.00560028 -43.13946473 0.20142084
             1 2 24.0 DavidsonDiagonalization() 0.00201 4.0 0.7 1.2 -42.97192905 -43.22189611 0.69794621
             1 3 24.0 DavidsonDiagonalization() 0.00201 3.0 0.7 1.4 -43.09499395 -43.09749186 0.00768862
             1 4 24.0 DavidsonDiagonalization() 7.69e-5 2.0 0.7 1.5 -43.09571104 -43.09617585 0.00118904
             1 5 24.0 DavidsonDiagonalization() 1.19e-5 3.0 0.7 1.7 -43.09622618 -43.09637952 0.00054718
             1 6 24.0 DavidsonDiagonalization() 5.47e-6 1.0 0.7 1.8 -43.09619459 -43.09625737 0.000193
             1 7 24.0 DavidsonDiagonalization() 1.93e-6 3.0 0.7 1.9 -43.0962549 -43.09626006 1.788e-5
             1 8 24.0 DavidsonDiagonalization() 1.79e-7 2.0 0.7 2.1 -43.09625733 -43.09625777 2.56e-6
             1 9 24.0 DavidsonDiagonalization() 2.56e-8 3.0 0.7 2.2 nothing nothing nothing
             2 1 24.0 DavidsonDiagonalization() 1.0e-6 5.0 0.7 2.6 -43.10825672 -43.11074971 0.00435174
             2 2 24.0 DavidsonDiagonalization() 4.35e-5 2.0 0.7 2.7 -43.10912901 -43.10942463 0.00053892
             2 3 24.0 DavidsonDiagonalization() 5.39e-6 2.0 0.7 2.9 -43.10924328 -43.10925158 2.323e-5
             2 4 24.0 DavidsonDiagonalization() 2.32e-7 4.0 0.7 3.0 -43.10925024 -43.10928148 0.00012258
             2 5 24.0 DavidsonDiagonalization() 2.32e-7 3.0 0.7 3.2 -43.10925169 -43.10925836 1.614e-5
             2 6 24.0 DavidsonDiagonalization() 1.61e-7 3.0 0.7 3.3 nothing nothing nothing
             3 1 24.0 DavidsonDiagonalization() 1.0e-6 5.0 0.7 3.7 -43.09901792 -43.10284311 0.00652404
             3 2 24.0 DavidsonDiagonalization() 6.52e-5 2.0 0.7 3.8 -43.10034879 -43.10058877 0.00048248
             3 3 24.0 DavidsonDiagonalization() 4.82e-6 2.0 0.7 4.0 -43.10043294 -43.10046987 6.432e-5
             3 4 24.0 DavidsonDiagonalization() 6.43e-7 3.0 0.7 4.1 -43.10044299 -43.10046877 6.082e-5
             3 5 24.0 DavidsonDiagonalization() 6.08e-7 2.0 0.7 4.2 nothing nothing nothing
             4 1 24.0 DavidsonDiagonalization() 1.0e-6 5.0 0.7 4.6 -43.10834553 -43.10952579 0.00199952
             4 2 24.0 DavidsonDiagonalization() 2.0e-5 2.0 0.7 4.8 -43.10876348 -43.10883933 0.00015055
             4 3 24.0 DavidsonDiagonalization() 1.51e-6 2.0 0.7 4.9 -43.10879034 -43.10880265 2.306e-5
             4 4 24.0 DavidsonDiagonalization() 2.31e-7 3.0 0.7 5.0 -43.10879483 -43.10880208 1.729e-5
             4 5 24.0 DavidsonDiagonalization() 1.73e-7 2.0 0.7 5.2 nothing nothing nothing
             5 1 24.0 DavidsonDiagonalization() 1.0e-6 5.0 0.7 5.6 -43.10753695 -43.10895232 0.00243803
             5 2 24.0 DavidsonDiagonalization() 2.44e-5 2.0 0.7 5.7 -43.10804643 -43.10816059 0.0002229
             5 3 24.0 DavidsonDiagonalization() 2.23e-6 2.0 0.7 5.8 -43.10808784 -43.10809655 1.679e-5
             5 4 24.0 DavidsonDiagonalization() 1.68e-7 4.0 0.7 6.0 -43.10808706 -43.10810564 5.311e-5
             5 5 24.0 DavidsonDiagonalization() 1.68e-7 3.0 0.7 6.1 nothing nothing nothing
             6 1 24.0 DavidsonDiagonalization() 1.0e-6 3.0 0.7 6.5 -43.10898535 -43.10903872 9.066e-5
             6 2 24.0 DavidsonDiagonalization() 9.07e-7 2.0 0.7 6.6 -43.1090036 -43.10901149 1.401e-5
             6 3 24.0 DavidsonDiagonalization() 1.4e-7 2.0 0.7 6.8 nothing nothing nothing
             7 1 24.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 7.2 -43.10969039 -43.10988487 0.00033653
             7 2 24.0 DavidsonDiagonalization() 3.37e-6 2.0 0.7 7.3 -43.10975976 -43.1097782 3.482e-5
             7 3 24.0 DavidsonDiagonalization() 3.48e-7 2.0 0.7 7.4 -43.10976646 -43.10976749 2.17e-6
             7 4 24.0 DavidsonDiagonalization() 2.17e-8 3.0 0.7 7.6 -43.10976648 -43.10976919 8.37e-6
             7 5 24.0 DavidsonDiagonalization() 2.17e-8 4.0 0.7 7.7 nothing nothing nothing
             8 1 24.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 8.1 -43.10894342 -43.10915529 0.00036499
             8 2 24.0 DavidsonDiagonalization() 3.65e-6 2.0 0.7 8.2 -43.10901923 -43.109033 2.773e-5
             8 3 24.0 DavidsonDiagonalization() 2.77e-7 2.0 0.7 8.4 -43.10902396 -43.10902655 4.42e-6
             8 4 24.0 DavidsonDiagonalization() 4.42e-8 4.0 0.7 8.5 -43.1090247 -43.10902628 3.82e-6
             8 5 24.0 DavidsonDiagonalization() 3.82e-8 2.0 0.7 8.6 nothing nothing nothing
             9 1 24.0 DavidsonDiagonalization() 1.0e-6 3.0 0.7 9.0 -43.10955328 -43.10960346 8.548e-5
             9 2 24.0 DavidsonDiagonalization() 8.55e-7 2.0 0.7 9.2 -43.10957201 -43.10957598 7.68e-6
             9 3 24.0 DavidsonDiagonalization() 7.68e-8 2.0 0.7 9.3 nothing nothing nothing
             10 1 24.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 9.7 -43.10970375 -43.10980381 0.00017021
             10 2 24.0 DavidsonDiagonalization() 1.7e-6 2.0 0.7 9.8 -43.10973953 -43.10975061 2.032e-5
             10 3 24.0 DavidsonDiagonalization() 2.03e-7 2.0 0.7 10.0 nothing nothing nothing
             11 1 24.0 DavidsonDiagonalization() 1.0e-6 2.0 0.7 10.4 -43.10975884 -43.10976234 5.12e-6
             11 2 24.0 DavidsonDiagonalization() 5.12e-8 2.0 0.7 10.5 -43.10975993 -43.10976118 2.18e-6
             11 3 24.0 DavidsonDiagonalization() 2.18e-8 2.0 0.7 10.7 nothing nothing nothing
             12 1 24.0 DavidsonDiagonalization() 1.0e-6 3.0 0.7 11.1 -43.10976664 -43.10976925 4.59e-6
             12 2 24.0 DavidsonDiagonalization() 4.59e-8 2.0 0.7 11.3 nothing nothing nothing
             13 1 24.0 DavidsonDiagonalization() 1.0e-6 2.0 0.7 11.7 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-43.09625738,
        -43.0962577,
        3.9e-7,
        nothing,
        nothing,
        nothing,
        nothing,)

    @test parse_version(str) == "5.2.1"

    @test parse_parallel_info(str) == ("Serial version", 1)

    @test parse_fft_dimensions(str) == (25271, (nr1 = 45, nr2 = 45, nr3 = 45))

    @test parse_bands(str) == ([
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
        ],
        [
         -27.899 -11.819 -13.383 -8.3482 -11.3791
         -13.4027 -8.2731 -11.265 -28.922 -11.3791
         -10.8557 -29.051 -11.265 -13.3805 -8.3809
         -10.8557 -13.3798 -8.4098 -11.4526 -28.8213
         -8.5036 -11.5296 -28.8126 -11.4526 -13.3817
         -28.647 -11.5296 -13.3815 -8.3634 -11.3917
         -13.3852 -8.3451 -11.3867 -28.7792 -11.3917
         -11.289 -28.5063 -11.3867 -13.3817 -8.3738
         -11.289 -13.3873 -8.3802 -11.3676 -28.8231
         -8.4016 -11.2059 -29.0227 -11.3676 -13.383
         -29.5199 -11.2059 -13.3812 -8.383 -11.3935
         -13.3829 -8.422 -11.5131 -28.7994 -11.3935
         -11.819 -28.6088 -11.5131 -13.3835 -8.381
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 0.86 0.88 1
         "" "electrons" 7.88 7.95 13
         "" "update_pot" 0.9 0.93 12
         "" "forces" 1.03 1.02 13
         "init_run" "wfcinit" 0.0 0.01 1
         "init_run" "potinit" 0.04 0.05 1
         "electrons" "c_bands" 1.2 1.21 58
         "electrons" "sum_band" 3.48 3.5 58
         "electrons" "v_of_rho" 0.89 0.89 68
         "electrons" "newd" 2.08 2.1 68
         "electrons" "mix_rho" 0.36 0.37 58
         "c_bands" "init_us_2" 0.09 0.1 117
         "c_bands" "regterg" 1.1 1.09 58
         "sum_band" "sum_band:bec" 0.0 0.0 58
         "sum_band" "addusdens" 2.8 2.81 58
         "*egterg" "h_psi" 0.92 0.88 213
         "*egterg" "s_psi" 0.02 0.02 213
         "*egterg" "g_psi" 0.05 0.04 154
         "*egterg" "rdiaghg" 0.02 0.02 197
         "h_psi" "add_vuspsi" 0.02 0.02 213
         "General routines" "calbec" 0.05 0.05 323
         "General routines" "fft" 0.88 0.94 610
         "General routines" "ffts" 0.17 0.13 126
         "General routines" "fftw" 0.76 0.73 1276
         "General routines" "interpolate" 0.49 0.48 126
         "General routines" "davcio" 0.0 0.0 13
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "/home/giannozz/trunk/espresso/PW/tests/relax-damped.in"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse relax H2O output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/cluster_example/reference/h2o.out-12"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 3,
        alat = 12.0,
        omega = 1728.0,
        nat = 3,
        ntyp = 2,
        nelec = 8.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 8,
        ecutwfc = 30.0,
        ecutrho = 120.0,
        ecutfock = nothing,
        conv_thr = 1.0e-7,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA PW PBX PBC ( 1  4  3  4 0 0)",
        nstep = 50,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Sum" 1369 1369 349
         "gvecs" "Sum" 38401 38401 4801
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [0.0 0.0 0.0 2.0], cryst = nothing)

    @test all(isempty, parse_stress(str))

    @test isempty(parse_cell_parameters(str))

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.071762144, 1.071762144, 1.079345568], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.071762144, 1.071762144, 1.079345568], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.056119834, 1.056119834, 1.077852061], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.056119834, 1.056119834, 1.077852061], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.043633568, 1.043633568, 1.087578917], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.043633568, 1.043633568, 1.087578917], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.033898772, 1.033898772, 1.105355925], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.033898772, 1.033898772, 1.105355925], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.028292878, 1.028292878, 1.124227412], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.028292878, 1.028292878, 1.124227412], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.029112625, 1.029112625, 1.126614785], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.029112625, 1.029112625, 1.126614785], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("O", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.029112625, 1.029112625, 1.126614785], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.029112625, 1.029112625, 1.126614785], [1, 1, 1])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 30.0 DavidsonDiagonalization() 0.01 8.0 0.7 2.4 -43.77716297 -44.16056246 0.48854594
             1 2 30.0 DavidsonDiagonalization() 0.00611 2.0 0.7 3.1 -43.8883681 -44.11207947 0.4291763
             1 3 30.0 DavidsonDiagonalization() 0.00536 2.0 0.7 3.7 -43.98513387 -43.9870297 0.0046026
             1 4 30.0 DavidsonDiagonalization() 5.75e-5 14.0 0.7 4.4 -43.98712771 -43.98722722 0.00034056
             1 5 30.0 DavidsonDiagonalization() 4.26e-6 3.0 0.7 5.0 -43.98713333 -43.98714858 6.155e-5
             1 6 30.0 DavidsonDiagonalization() 7.69e-7 2.0 0.7 5.6 -43.98713711 -43.98713768 3.74e-6
             1 7 30.0 DavidsonDiagonalization() 4.68e-8 2.0 0.7 6.3 nothing nothing nothing
             2 1 30.0 DavidsonDiagonalization() 1.0e-6 9.0 0.7 8.3 -43.99392398 -44.00235061 0.0114905
             2 2 30.0 DavidsonDiagonalization() 0.000144 2.0 0.7 8.9 -43.99647022 -44.00286894 0.01269539
             2 3 30.0 DavidsonDiagonalization() 0.000144 2.0 0.7 9.6 -43.99912253 -43.99912355 4.86e-5
             2 4 30.0 DavidsonDiagonalization() 6.07e-7 6.0 0.7 10.3 -43.999157 -43.99915768 3.64e-6
             2 5 30.0 DavidsonDiagonalization() 4.55e-8 1.0 0.7 10.9 -43.99915732 -43.99915735 2.9e-7
             2 6 30.0 DavidsonDiagonalization() 3.64e-9 3.0 0.7 11.5 nothing nothing nothing
             3 1 30.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 13.5 -43.99991449 -44.00002502 0.00017069
             3 2 30.0 DavidsonDiagonalization() 2.13e-6 3.0 0.7 14.1 -43.99994728 -44.00002898 0.00015697
             3 3 30.0 DavidsonDiagonalization() 1.96e-6 2.0 0.7 14.7 -43.99998314 -43.99998331 1.67e-6
             3 4 30.0 DavidsonDiagonalization() 2.09e-8 2.0 0.7 15.4 -43.99998378 -43.99998378 4.0e-8
             3 5 30.0 DavidsonDiagonalization() 5.42e-10 2.0 0.7 16.0 nothing nothing nothing
             4 1 30.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 18.0 -44.00029805 -44.0003096 3.838e-5
             4 2 30.0 DavidsonDiagonalization() 4.8e-7 2.0 0.7 18.6 -44.00030292 -44.00031303 1.848e-5
             4 3 30.0 DavidsonDiagonalization() 2.31e-7 2.0 0.7 19.3 -44.00030821 -44.00030867 1.8e-6
             4 4 30.0 DavidsonDiagonalization() 2.25e-8 2.0 0.7 19.9 -44.00030857 -44.00030856 1.0e-8
             4 5 30.0 DavidsonDiagonalization() 1.29e-10 3.0 0.7 20.5 nothing nothing nothing
             5 1 30.0 DavidsonDiagonalization() 1.0e-6 5.0 0.7 22.5 -44.00053011 -44.00052409 1.746e-5
             5 2 30.0 DavidsonDiagonalization() 2.18e-7 2.0 0.7 23.1 -44.00053186 -44.00053217 3.04e-6
             5 3 30.0 DavidsonDiagonalization() 3.8e-8 2.0 0.7 23.7 -44.0005324 -44.00053298 1.37e-6
             5 4 30.0 DavidsonDiagonalization() 1.71e-8 2.0 0.7 24.4 -44.00053271 -44.00053271 1.0e-8
             5 5 30.0 DavidsonDiagonalization() 1.75e-10 2.0 0.7 25.0 nothing nothing nothing
             6 1 30.0 DavidsonDiagonalization() 1.0e-6 3.0 0.7 27.0 -44.00063016 -44.00064419 3.586e-5
             6 2 30.0 DavidsonDiagonalization() 4.48e-7 2.0 0.7 27.6 -44.00063601 -44.00065001 2.798e-5
             6 3 30.0 DavidsonDiagonalization() 3.5e-7 2.0 0.7 28.2 -44.00064284 -44.0006427 8.4e-7
             6 4 30.0 DavidsonDiagonalization() 1.04e-8 2.0 0.7 28.9 -44.00064301 -44.00064301 7.9e-9
             6 5 30.0 DavidsonDiagonalization() 9.89e-11 2.0 0.7 29.5 nothing nothing nothing
             7 1 30.0 DavidsonDiagonalization() 1.0e-6 2.0 0.7 31.5 -44.00064977 -44.00065223 3.19e-6
             7 2 30.0 DavidsonDiagonalization() 3.98e-8 2.0 0.7 32.1 -44.00065043 -44.00065198 3.01e-6
             7 3 30.0 DavidsonDiagonalization() 3.77e-8 2.0 0.7 32.7 -44.0006511 -44.0006511 3.0e-8
             7 4 30.0 DavidsonDiagonalization() 3.16e-10 3.0 0.7 33.4 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-43.98713779,
        -43.98713782,
        8.0e-8,
        -152.747884,
        [-83.31808498, 43.20085074, -8.51939717, 14.56351319],
        nothing,
        nothing,)

    @test parse_version(str) == "6.2"

    @test parse_parallel_info(str) == ("Serial version", 1)

    @test parse_fft_dimensions(str) == (19201, (nr1 = 45, nr2 = 45, nr3 = 45))

    @test parse_bands(str) == ([
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
        ],
        [
         -25.7682 2.6719 1.9889 1.7196 -1.4183 -7.1592 -9.2013 -13.0172
         -13.8217 -24.9814 2.6085 2.0181 1.7413 -1.4185 -7.1593 -9.203
         -9.0546 -13.125 -25.1028 2.6215 2.0263 1.7436 -1.4249 -7.1569
         -7.2702 -8.9632 -13.1874 -25.1742 2.6347 2.0192 1.7316 -1.4277
         -1.3182 -7.1093 -9.0127 -13.1787 -25.2037 2.6457 2.004 1.7255
         1.9473 -1.4483 -7.135 -9.0749 -13.116 -25.1968 2.6539 2.0
         2.1707 1.6724 -1.4273 -7.151 -9.1444 -13.0324 -25.184 2.6531
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 1.15 1.15 1
         "" "electrons" 23.97 23.98 7
         "" "update_pot" 3.66 3.66 6
         "" "forces" 2.71 2.71 7
         "init_run" "wfcinit" 0.02 0.01 1
         "init_run" "potinit" 0.57 0.57 1
         "electrons" "c_bands" 1.56 1.57 38
         "electrons" "sum_band" 1.11 1.11 38
         "electrons" "v_of_rho" 13.25 13.26 44
         "electrons" "newd" 0.63 0.62 44
         "electrons" "PAW_pot" 12.29 12.29 50
         "electrons" "mix_rho" 0.13 0.13 38
         "c_bands" "init_us_2" 0.05 0.07 83
         "c_bands" "regterg" 1.52 1.52 38
         "sum_band" "sum_band:bec" 0.0 0.0 44
         "sum_band" "addusdens" 0.7 0.72 38
         "*egterg" "h_psi" 1.32 1.35 155
         "*egterg" "s_psi" 0.01 0.01 155
         "*egterg" "g_psi" 0.01 0.01 116
         "*egterg" "rdiaghg" 0.02 0.03 147
         "h_psi" "h_psi:pot" 1.31 1.34 155
         "h_psi" "h_psi:calbec" 0.03 0.02 155
         "h_psi" "vloc_psi" 1.26 1.3 155
         "h_psi" "add_vuspsi" 0.02 0.01 155
         "General routines" "calbec" 0.06 0.04 227
         "General routines" "fft" 0.9 0.88 658
         "General routines" "fftw" 1.33 1.37 1072
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse relax NH4 output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/cluster_example/reference/nh4%2B.out-12"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 5,
        alat = 12.0,
        omega = 1728.0,
        nat = 5,
        ntyp = 2,
        nelec = 8.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 8,
        ecutwfc = 30.0,
        ecutrho = 120.0,
        ecutfock = nothing,
        conv_thr = 1.0e-7,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA PW PBX PBC ( 1  4  3  4 0 0)",
        nstep = 50,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Sum" 1369 1369 349
         "gvecs" "Sum" 38401 38401 4801
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [0.0 0.0 0.0 2.0], cryst = nothing)

    @test all(isempty, parse_stress(str))

    @test isempty(parse_cell_parameters(str))

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("N", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.154573639, 1.154573639, 1.154573639], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.154573639, 1.154573639, 1.154573639], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.154573639, 1.154573639, 1.154573639], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.154573639, 1.154573639, 1.154573639], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("N", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.137353203, 1.137353203, 1.137353203], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.137353203, 1.137353203, 1.137353203], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.137353203, 1.137353203, 1.137353203], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.137353203, 1.137353203, 1.137353203], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("N", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129023295, 1.129023295, 1.129023295], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129023295, 1.129023295, 1.129023295], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129023295, 1.129023295, 1.129023295], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129023295, 1.129023295, 1.129023295], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("N", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1])]),
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("bohr", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("N", [0.0, 0.0, 0.0], [0, 0, 0]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("H", [1.129577564, 1.129577564, 1.129577564], [1, 1, 1])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 30.0 DavidsonDiagonalization() 0.01 3.0 0.7 2.8 -31.57695851 -33.30284118 2.34289398
             1 2 30.0 DavidsonDiagonalization() 0.01 2.0 0.7 3.4 -32.20709529 -32.59490336 0.73212995
             1 3 30.0 DavidsonDiagonalization() 0.00915 2.0 0.7 4.1 -32.34001777 -32.3467286 0.01300493
             1 4 30.0 DavidsonDiagonalization() 0.000163 5.0 0.7 4.7 -32.34428091 -32.34504157 0.0016397
             1 5 30.0 DavidsonDiagonalization() 2.05e-5 3.0 0.7 5.4 -32.34435726 -32.34436948 2.886e-5
             1 6 30.0 DavidsonDiagonalization() 3.61e-7 4.0 0.7 6.1 -32.3443684 -32.34438431 3.706e-5
             1 7 30.0 DavidsonDiagonalization() 3.61e-7 1.0 0.7 6.7 -32.34437177 -32.34437186 4.5e-7
             1 8 30.0 DavidsonDiagonalization() 5.66e-9 3.0 0.7 7.4 nothing nothing nothing
             2 1 30.0 DavidsonDiagonalization() 1.0e-6 13.0 0.7 9.5 -32.40858339 -32.47321777 0.10053878
             2 2 30.0 DavidsonDiagonalization() 0.00126 2.0 0.7 10.2 -32.43449921 -32.46769906 0.06689062
             2 3 30.0 DavidsonDiagonalization() 0.000836 1.0 0.7 10.8 -32.44749827 -32.44741407 0.0003074
             2 4 30.0 DavidsonDiagonalization() 3.84e-6 3.0 0.7 11.5 -32.44777934 -32.44778967 5.521e-5
             2 5 30.0 DavidsonDiagonalization() 6.9e-7 1.0 0.7 12.2 -32.44777248 -32.44778205 2.146e-5
             2 6 30.0 DavidsonDiagonalization() 2.68e-7 2.0 0.7 12.8 -32.44777767 -32.44777965 4.63e-6
             2 7 30.0 DavidsonDiagonalization() 5.79e-8 1.0 0.7 13.5 -32.44777827 -32.44777831 1.1e-7
             2 8 30.0 DavidsonDiagonalization() 1.33e-9 3.0 0.7 14.1 nothing nothing nothing
             3 1 30.0 DavidsonDiagonalization() 1.0e-6 4.0 0.7 16.2 -32.45016819 -32.45062893 0.00077953
             3 2 30.0 DavidsonDiagonalization() 9.74e-6 2.0 0.7 16.8 -32.45036021 -32.45058894 0.00045495
             3 3 30.0 DavidsonDiagonalization() 5.69e-6 2.0 0.7 17.5 -32.45045353 -32.45045249 3.68e-6
             3 4 30.0 DavidsonDiagonalization() 4.6e-8 2.0 0.7 18.1 -32.45045437 -32.45045439 8.0e-8
             3 5 30.0 DavidsonDiagonalization() 9.49e-10 2.0 0.7 18.8 -32.45045441 -32.45045441 3.0e-8
             3 6 30.0 DavidsonDiagonalization() 3.31e-10 1.0 0.7 19.4 nothing nothing nothing
             4 1 30.0 DavidsonDiagonalization() 1.0e-6 3.0 0.7 21.4 -32.45068208 -32.45079172 0.00018723
             4 2 30.0 DavidsonDiagonalization() 2.34e-6 2.0 0.7 22.1 -32.45072843 -32.45078037 0.0001025
             4 3 30.0 DavidsonDiagonalization() 1.28e-6 2.0 0.7 22.7 -32.45074971 -32.45074947 8.8e-7
             4 4 30.0 DavidsonDiagonalization() 1.1e-8 2.0 0.7 23.4 -32.45074993 -32.45074993 1.0e-8
             4 5 30.0 DavidsonDiagonalization() 1.65e-10 3.0 0.7 24.0 nothing nothing nothing
             5 1 30.0 DavidsonDiagonalization() 1.0e-6 2.0 0.7 26.1 -32.45075107 -32.45075183 9.2e-7
             5 2 30.0 DavidsonDiagonalization() 1.15e-8 2.0 0.7 26.7 -32.4507513 -32.45075156 4.9e-7
             5 3 30.0 DavidsonDiagonalization() 6.18e-9 2.0 0.7 27.4 -32.45075141 -32.45075141 4.2e-9
             5 4 30.0 DavidsonDiagonalization() 5.26e-11 2.0 0.7 28.0 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-32.34437204,
        -32.34437208,
        5.0e-8,
        -113.643173,
        [-82.06686425, 38.91703901, -8.21266386, 27.33665144],
        nothing,
        nothing,)

    @test parse_version(str) == "6.2"

    @test parse_parallel_info(str) == ("Serial version", 1)

    @test parse_fft_dimensions(str) == (19201, (nr1 = 45, nr2 = 45, nr3 = 45))

    @test parse_bands(str) == ([
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
         0.0 0.0 0.0
        ],
        [
         -33.6496 -3.7177 -20.5654 -4.7342 -7.443 -20.8479 -4.5444 -20.842
         -22.3919 -3.7177 -20.5654 -31.4759 -4.6059 -20.8479 -4.5444 -7.4142
         -22.3919 -3.7177 -7.5142 -20.7553 -4.6059 -20.8479 -31.5898 -4.5487
         -22.3919 -31.2273 -4.7342 -20.7553 -4.6059 -7.412 -20.842 -4.5487
         -7.041 -20.5654 -4.7342 -20.7553 -31.5977 -4.5444 -20.842 -4.5487
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 1.22 1.31 1
         "" "electrons" 20.66 20.75 5
         "" "update_pot" 2.52 2.52 4
         "" "forces" 2.04 2.04 5
         "init_run" "wfcinit" 0.01 0.02 1
         "init_run" "potinit" 0.65 0.73 1
         "electrons" "c_bands" 1.39 1.48 32
         "electrons" "sum_band" 0.97 0.98 32
         "electrons" "v_of_rho" 11.13 11.2 36
         "electrons" "newd" 0.53 0.54 36
         "electrons" "PAW_pot" 10.2 10.2 40
         "electrons" "mix_rho" 0.14 0.14 32
         "c_bands" "init_us_2" 0.08 0.06 69
         "c_bands" "regterg" 1.34 1.43 32
         "sum_band" "sum_band:bec" 0.0 0.0 36
         "sum_band" "addusdens" 0.62 0.62 32
         "*egterg" "h_psi" 1.18 1.2 116
         "*egterg" "s_psi" 0.02 0.02 116
         "*egterg" "g_psi" 0.02 0.01 83
         "*egterg" "rdiaghg" 0.02 0.1 110
         "h_psi" "h_psi:pot" 1.17 1.2 116
         "h_psi" "h_psi:calbec" 0.02 0.02 116
         "h_psi" "vloc_psi" 1.13 1.15 116
         "h_psi" "add_vuspsi" 0.02 0.02 116
         "General routines" "calbec" 0.03 0.04 172
         "General routines" "fft" 0.7 0.74 530
         "General routines" "fftw" 1.19 1.21 916
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse relax MoS output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/gatefield/reference/single_%2B0.10.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 3,
        alat = 5.9716,
        omega = 2213.0132,
        nat = 3,
        ntyp = 2,
        nelec = 25.9,
        nelup = nothing,
        neldw = nothing,
        nbnd = 17,
        ecutwfc = 50.0,
        ecutrho = 410.0,
        ecutfock = nothing,
        conv_thr = 1.0e-9,
        mixing_beta = 0.7,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "LDA ( 1  1  0  0 0 0)",
        nstep = 300,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 253 124 37
         "gvecs" "Min" 77613 26474 4339
         "sticks" "Max" 255 125 38
         "gvecs" "Max" 77637 26515 4364
         "sticks" "Sum" 1015 499 151
         "gvecs" "Sum" 310487 105989 17427
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
            0.0 0.0 0.0 0.0078125
            0.0 0.0721688 0.0 0.046875
            0.0 0.1443376 0.0 0.046875
            0.0 0.2165064 0.0 0.046875
            0.0 0.2886751 0.0 0.046875
            0.0 0.3608439 0.0 0.046875
            0.0 0.4330127 0.0 0.046875
            0.0 0.5051815 0.0 0.046875
            0.0 -0.5773503 0.0 0.0234375
            0.0625 0.1082532 0.0 0.046875
            0.0625 0.180422 0.0 0.09375
            0.0625 0.2525907 0.0 0.09375
            0.0625 0.3247595 0.0 0.09375
            0.0625 0.3969283 0.0 0.09375
            0.0625 0.4690971 0.0 0.09375
            0.0625 0.5412659 0.0 0.09375
            0.125 0.2165064 0.0 0.046875
            0.125 0.2886751 0.0 0.09375
            0.125 0.3608439 0.0 0.09375
            0.125 0.4330127 0.0 0.09375
            0.125 0.5051815 0.0 0.09375
            0.125 0.5773503 0.0 0.046875
            0.1875 0.3247595 0.0 0.046875
            0.1875 0.3969283 0.0 0.09375
            0.1875 0.4690971 0.0 0.09375
            0.1875 0.5412659 0.0 0.09375
            0.25 0.4330127 0.0 0.046875
            0.25 0.5051815 0.0 0.09375
            0.25 0.5773503 0.0 0.046875
            0.3125 0.5412659 0.0 0.046875
        ],
        cryst = [
            0.0 0.0 0.0 0.0078125
            0.0 0.0625 0.0 0.046875
            0.0 0.125 0.0 0.046875
            0.0 0.1875 0.0 0.046875
            0.0 0.25 0.0 0.046875
            0.0 0.3125 0.0 0.046875
            0.0 0.375 0.0 0.046875
            0.0 0.4375 0.0 0.046875
            0.0 -0.5 0.0 0.0234375
            0.0625 0.0625 0.0 0.046875
            0.0625 0.125 0.0 0.09375
            0.0625 0.1875 0.0 0.09375
            0.0625 0.25 0.0 0.09375
            0.0625 0.3125 0.0 0.09375
            0.0625 0.375 0.0 0.09375
            0.0625 0.4375 0.0 0.09375
            0.125 0.125 0.0 0.046875
            0.125 0.1875 0.0 0.09375
            0.125 0.25 0.0 0.09375
            0.125 0.3125 0.0 0.09375
            0.125 0.375 0.0 0.09375
            0.125 0.4375 0.0 0.046875
            0.1875 0.1875 0.0 0.046875
            0.1875 0.25 0.0 0.09375
            0.1875 0.3125 0.0 0.09375
            0.1875 0.375 0.0 0.09375
            0.25 0.25 0.0 0.046875
            0.25 0.3125 0.0 0.09375
            0.25 0.375 0.0 0.046875
            0.3125 0.3125 0.0 0.046875
        ],)

    @test all(isempty, parse_stress(str))

    @test isempty(parse_cell_parameters(str))

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("alat", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("S", [0.5, 0.28867513, 1.86331695], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("Mo", [0.0, 0.57735027, 2.350404949], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("S", [0.0, 0.57735027, 2.838203782], [1, 1, 1])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 50.0 DavidsonDiagonalization() 0.01 3.0 0.7 10.1 -180.35021564 -181.15850969 0.95068425
             1 2 50.0 DavidsonDiagonalization() 0.00367 5.1 0.7 19.4 -177.05589334 -184.83602299 52.96038171
             1 3 50.0 DavidsonDiagonalization() 0.00367 4.7 0.7 28.7 -180.94555842 -181.08909805 0.67171122
             1 4 50.0 DavidsonDiagonalization() 0.00259 2.0 0.7 34.3 -181.00685439 -181.05351213 0.19909785
             1 5 50.0 DavidsonDiagonalization() 0.000769 1.2 0.7 39.8 -181.02313707 -181.02761805 0.02297994
             1 6 50.0 DavidsonDiagonalization() 8.87e-5 3.0 0.7 46.0 -181.0250963 -181.02591062 0.00231326
             1 7 50.0 DavidsonDiagonalization() 8.93e-6 3.0 0.7 53.1 -181.02558567 -181.0257633 0.00049502
             1 8 50.0 DavidsonDiagonalization() 1.91e-6 2.1 0.7 58.9 -181.0256008 -181.02566663 0.00010457
             1 9 50.0 DavidsonDiagonalization() 4.04e-7 2.4 0.7 64.9 -181.02561878 -181.0256278 1.47e-6
             1 10 50.0 DavidsonDiagonalization() 5.69e-9 4.0 0.7 73.5 -181.02562083 -181.02561807 1.11e-6
             1 11 50.0 DavidsonDiagonalization() 4.3e-9 1.6 0.7 79.0 -181.02562089 -181.02562103 1.0e-7
             1 12 50.0 DavidsonDiagonalization() 3.92e-10 2.6 0.7 85.3 -181.02562091 -181.02562113 5.0e-9
             1 13 50.0 DavidsonDiagonalization() 1.95e-11 3.0 0.7 92.7 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-181.02562091,
        -181.02562123,
        8.8e-10,
        nothing,
        [-2649.72274889, 1320.02437487, -42.50453186, 1179.75742402],
        nothing,
        nothing,)

    @test parse_version(str) == "6.1"

    @test parse_parallel_info(str) == ("Parallel version (MPI)", 4)

    @test parse_fft_dimensions(str) == (310487, (nr1 = 40, nr2 = 40, nr3 = 480))

    @test parse_bands(str) == ([
         0.0 0.0625 0.125
         0.0 0.1804 0.5052
         0.0 0.0 0.0
         0.0 0.0625 0.125
         0.0722 0.2526 0.5774
         0.0 0.0 0.0
         0.0 0.0625 0.1875
         0.1443 0.3248 0.3248
         0.0 0.0 0.0
         0.0 0.0625 0.1875
         0.2165 0.3969 0.3969
         0.0 0.0 0.0
         0.0 0.0625 0.1875
         0.2887 0.4691 0.4691
         0.0 0.0 0.0
         0.0 0.0625 0.1875
         0.3608 0.5413 0.5413
         0.0 0.0 0.0
         0.0 0.125 0.25
         0.433 0.2165 0.433
         0.0 0.0 0.0
         0.0 0.125 0.25
         0.5052 0.2887 0.5052
         0.0 0.0 0.0
         0.0 0.125 0.25
         -0.5774 0.3608 0.5774
         0.0 0.0 0.0
         0.0625 0.125 0.3125
         0.1083 0.433 0.5413
         0.0 0.0 0.0
        ],
        [
         -64.7819 -1.2874 -6.3067 -16.2998 -39.1479 -1.2315 -6.2811 -9.1997 -38.9482 -0.0012 -4.9766 -8.7019 -38.8842 0.4732 -4.5836 -7.6257 -16.0426
         -39.0266 -1.0905 -6.2323 -9.3376 -38.9452 -0.0262 -4.9045 -7.9744 -38.8475 0.3004 -4.7871 -7.8813 -16.1887 -64.7443 -1.7815 -7.298 -15.8795
         -38.7512 -1.0196 -5.0293 -8.1597 -38.8425 0.4439 -4.2625 -7.8771 -16.4436 -64.7651 -1.9696 -7.3061 -15.9812 -39.0709 -1.1695 -7.0076 -9.5251
         -38.7512 -0.9175 -4.3849 -7.937 -16.3557 -64.7772 -1.6326 -6.9413 -16.0988 -38.9932 -1.4907 -6.0391 -9.7235 -38.9518 -0.6909 -5.5503 -9.1006
         -18.2121 -64.7757 -1.6819 -6.9177 -16.0691 -39.0173 -1.4252 -5.8011 -9.7156 -38.9109 -1.2815 -5.4725 -9.1176 -38.8834 0.647 -4.6059 -8.0056
         -16.775 -39.0143 -1.3629 -5.536 -9.8643 -38.7979 -1.2898 -5.2374 -8.7418 -38.796 0.2069 -4.7507 -8.094 -16.5515 -64.7368 -1.8015 -7.2664
         -9.7541 -38.8132 -1.3052 -5.4681 -8.8083 -38.7624 -0.5098 -4.892 -8.0257 -17.5406 -64.7501 -1.7705 -7.2728 -16.1244 -39.1193 -1.155 -7.036
         -6.7889 -38.7652 -0.4403 -5.0104 -8.1089 -18.0323 -64.764 -1.8549 -7.2555 -16.5164 -39.0414 -1.4143 -5.9482 -9.3214 -38.9368 -0.8369 -5.9486
         -6.7889 -17.9742 -64.7607 -1.7257 -7.2799 -16.7059 -38.9909 -1.5575 -5.6822 -9.1533 -38.9633 -0.6173 -5.7545 -8.6769 -38.9236 0.7506 -4.2161
         -5.6493 -16.6836 -38.9845 -1.6552 -5.7602 -9.5572 -38.9254 -0.9568 -5.4296 -7.797 -38.8445 0.4815 -4.4386 -7.7134 -16.1338 -64.7362 -2.2026
         -5.6162 -9.4985 -38.9604 -0.6327 -5.0433 -7.154 -38.795 0.141 -4.7436 -7.2362 -16.8468 -64.7394 -1.7866 -7.3724 -15.9399 -39.0928 -0.7898
         -5.6162 -7.2584 -38.7989 0.1453 -4.7713 -6.5995 -17.4912 -64.7483 -1.7666 -6.793 -16.2473 -39.1314 -1.0658 -6.5024 -9.461 -38.9581 -0.5207
         -3.945 -6.6363 -17.3478 -64.7457 -1.7785 -6.4799 -16.4979 -39.0726 -1.3824 -6.1555 -9.2694 -38.942 -0.1714 -5.3078 -9.2876 -38.9355 0.8938
         -1.1886 -6.536 -16.4438 -39.1044 -1.3124 -5.9266 -9.1527 -38.9597 -0.2368 -4.8545 -8.236 -38.8793 0.3327 -4.7836 -7.9988 -16.1056 -64.7344
         -1.1886 -5.9654 -9.1605 -38.9545 -0.0893 -5.0141 -7.8399 -38.8334 0.3636 -4.5216 -7.8518 -16.2746 -64.7507 -1.7973 -7.1899 -15.9156 -39.0718
         -0.8565 -5.0047 -7.9227 -38.8323 0.4139 -4.058 -7.4072 -16.7537 -64.7394 -1.833 -7.2775 -16.0191 -39.0273 -1.4109 -6.6327 -9.4221 -39.0004
         -0.8565 -4.1066 -7.6606 -16.6124 -64.7394 -1.437 -6.6554 -16.217 -39.1546 -1.4821 -6.0289 -9.6622 -38.9645 -0.9188 -5.799 -9.1799 -38.932
         -64.7803 -1.4709 -6.6269 -16.1668 -39.1631 -1.25 -6.0749 -9.4539 -38.942 -1.4158 -5.2 -9.0101 -38.8515 0.6118 -4.382 -7.8186 -16.0113
         -39.0234 -1.3135 -5.8928 -9.6302 -38.942 -1.0883 -5.0075 -8.3875 -38.8551 -0.129 -4.8808 -8.0209 -16.8791 -64.7394 -1.931 -7.2622 -15.8566
         -38.7674 -1.1331 -5.2319 -8.5457 -38.8461 -0.8448 -4.5926 -7.9288 -16.2659 -64.7576 -1.9325 -7.2878 -16.2573 -39.1026 -0.8599 -7.1142 -9.6018
         -38.7549 -0.7828 -4.7556 -7.9945 -16.2623 -64.7715 -1.8009 -7.1368 -16.0279 -38.9783 -1.469 -5.9748 -9.1875 -38.942 -0.5541 -5.8188 -8.8632
         -18.1513 -64.7689 -1.7462 -7.1407 -16.0314 -39.0059 -1.4784 -5.6279 -9.873 -38.9767 -0.9686 -5.6821 -8.1556 -38.9092 0.6343 -4.356 -8.0724
         -16.7516 -39.0006 -1.5785 -5.6512 -9.9527 -38.8536 -1.2876 -5.4632 -8.9345 -38.82 0.421 -4.5459 -7.8145 -16.2834 -64.7394 -2.0511 -7.5466
         -9.6857 -38.8812 -1.0185 -5.2362 -8.9011 -38.7765 -0.1519 -4.9262 -8.1307 -17.2077 -64.7437 -1.7468 -7.3704 -16.0102 -39.0877 -0.8219 -6.9547
         -6.9211 -38.7807 -0.1127 -4.98 -8.1659 -17.8051 -64.7558 -1.8232 -7.3111 -16.3875 -39.0955 -1.207 -6.2045 -9.4178 -38.9421 -0.7019 -6.0271
         -6.7238 -17.6972 -64.7526 -1.7842 -7.3272 -16.6185 -39.0023 -1.5169 -5.7775 -9.1041 -38.9505 -0.3157 -5.0654 -9.0785 -38.9246 0.8539 -4.1327
         -6.1206 -16.5772 -39.0387 -1.4859 -5.7957 -9.3399 -38.9747 -0.5728 -5.3158 -7.9157 -38.8654 0.4095 -4.8682 -7.8221 -16.287 -64.735 -2.3367
         -5.7135 -9.2642 -38.9682 -0.2942 -4.9781 -7.5168 -38.815 0.3201 -4.5656 -7.6843 -16.5132 -64.7379 -1.9824 -7.2846 -16.0064 -39.0928 -0.8745
         -5.2287 -7.6572 -38.817 0.3184 -4.657 -6.9239 -17.1246 -64.7425 -1.748 -7.1454 -16.1165 -39.144 -1.4473 -6.6299 -9.3053 -38.9717 -0.3918
         -3.9775 -7.156 -16.9671 -64.741 -1.7658 -6.3757 -16.3578 -39.1259 -1.2303 -6.0117 -9.4942 -38.939 -1.1503 -5.5964 -9.1035 -38.9333 1.0046
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 2.9 3.0 1
         "" "electrons" 87.76 89.4 1
         "" "forces" 0.83 0.88 1
         "init_run" "wfcinit" 2.07 2.1 1
         "init_run" "wfcinit:atom" 0.03 0.03 30
         "init_run" "wfcinit:wfcr" 1.92 1.95 30
         "init_run" "potinit" 0.16 0.16 1
         "electrons" "c_bands" 64.67 65.15 13
         "electrons" "sum_band" 19.29 20.02 13
         "electrons" "v_of_rho" 0.38 0.4 14
         "electrons" "v_h" 0.16 0.17 14
         "electrons" "v_xc" 0.12 0.13 15
         "electrons" "newd" 2.93 3.4 14
         "electrons" "mix_rho" 0.47 0.47 13
         "c_bands" "init_us_2" 1.25 1.26 840
         "c_bands" "cegterg" 62.42 62.87 390
         "sum_band" "sum_band:bec" 0.01 0.01 390
         "sum_band" "addusdens" 4.79 5.42 13
         "*egterg" "h_psi" 56.04 56.46 1548
         "*egterg" "s_psi" 1.08 1.09 1548
         "*egterg" "g_psi" 0.22 0.22 1128
         "*egterg" "cdiaghg" 0.81 0.82 1518
         "*egterg" "cegterg:over" 1.95 1.97 1128
         "*egterg" "cegterg:upda" 1.45 1.46 1128
         "*egterg" "cegterg:last" 0.66 0.67 420
         "h_psi" "h_psi:pot" 55.69 56.1 1548
         "h_psi" "h_psi:calbec" 1.52 1.53 1548
         "h_psi" "vloc_psi" 53.1 53.5 1548
         "h_psi" "add_vuspsi" 1.06 1.07 1548
         "General routines" "calbec" 2.11 2.13 2058
         "General routines" "fft" 1.34 1.35 126
         "General routines" "ffts" 0.04 0.04 27
         "General routines" "fftw" 48.23 48.51 44506
         "General routines" "interpolate" 0.35 0.35 27
         "Parallel routines" "fft_scatt_xy" 3.85 3.88 44659
         "Parallel routines" "fft_scatt_yz" 9.8 9.87 44659
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end

@testset "Parse relax Al output" begin
    url = "https://raw.githubusercontent.com/QEF/q-e/master/PW/examples/ESM_example/reference/Al001_bc3_v00.out"
    str = open(download(url), "r") do io
        read(io, String)
    end

    @test isnothing(tryparse(SubroutineError, str))

    @test_throws Meta.ParseError parse(SubroutineError, str)

    @test tryparse(Preamble, str) == parse(Preamble, str) == Preamble(ibrav = 4,
        alat = 10.8223,
        omega = 2655.9321,
        nat = 4,
        ntyp = 1,
        nelec = 12.0,
        nelup = nothing,
        neldw = nothing,
        nbnd = 10,
        ecutwfc = 20.0,
        ecutrho = 80.0,
        ecutfock = nothing,
        conv_thr = 1.0e-6,
        mixing_beta = 0.3,
        mixing_ndim = 8,
        mixing_mode = "plain",
        xc = "SLA  PW   PBE  PBE ( 1  4  3  4 0 0)",
        nstep = 50,
    )

    @test parse_fft_base_info(str) == DataFrame([
         "sticks" "Min" 187 187 54
         "gvecs" "Min" 8037 8037 1261
         "sticks" "Max" 188 188 56
         "gvecs" "Max" 8044 8044 1262
         "sticks" "Sum" 749 749 221
         "gvecs" "Sum" 32157 32157 5047
        ],
        [:kind, :stats, :dense, :smooth, :PW],
    )

    @test parse_ibz(str) == (cart = [
            0.0833333 0.0833333 0.0 0.1111111
            0.0833333 0.25 0.0 0.1111111
            0.0833333 0.4166667 0.0 0.1111111
            0.0833333 -0.4166667 0.0 0.1111111
            0.0833333 -0.25 0.0 0.1111111
            0.0833333 -0.0833333 0.0 0.1111111
            0.25 0.0833333 0.0 0.1111111
            0.25 0.25 0.0 0.1111111
            0.25 0.4166667 0.0 0.1111111
            0.25 -0.4166667 0.0 0.1111111
            0.25 -0.25 0.0 0.1111111
            0.25 -0.0833333 0.0 0.1111111
            0.4166667 0.0833333 0.0 0.1111111
            0.4166667 0.25 0.0 0.1111111
            0.4166667 0.4166667 0.0 0.1111111
            0.4166667 -0.4166667 0.0 0.1111111
            0.4166667 -0.25 0.0 0.1111111
            0.4166667 -0.0833333 0.0 0.1111111
        ],
        cryst = nothing,)

    @test all(isempty, parse_stress(str))

    @test isempty(parse_cell_parameters(str))

    # @test parse_atomic_positions(str) == QuantumESPRESSOBase.Cards.AtomicPositionsCard[
    #     QuantumESPRESSOBase.Cards.AtomicPositionsCard{String,Array{QuantumESPRESSOBase.Cards.AtomicPosition,1}}("angstrom", QuantumESPRESSOBase.Cards.AtomicPosition[QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("Al", [0.0, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("Al", [2.863450038, 0.0, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("Al", [0.0, 2.863450038, 0.0], [1, 1, 1]), QuantumESPRESSOBase.Cards.AtomicPosition{String,Array{Float64,1},Array{Int64,1}}("Al", [2.863450038, 2.863450038, 0.0], [1, 1, 1])])
    # ]

    @test parse_iteration_head(str) == groupby(DataFrame([
             1 1 20.0 DavidsonDiagonalization() 0.01 2.9 0.3 2.6 -49.52188732 -49.53876762 0.0697686
             1 2 20.0 DavidsonDiagonalization() 0.000581 1.0 0.3 3.1 -49.51177719 -49.52342336 0.02663719
             1 3 20.0 DavidsonDiagonalization() 0.000222 1.0 0.3 3.7 -49.51482345 -49.51499369 0.0005093
             1 4 20.0 DavidsonDiagonalization() 4.24e-6 13.3 0.3 5.5 -49.51663021 -49.51666358 0.00011418
             1 5 20.0 DavidsonDiagonalization() 9.51e-7 6.4 0.3 6.5 -49.51662162 -49.51673099 0.00018936
             1 6 20.0 DavidsonDiagonalization() 9.51e-7 3.6 0.3 7.2 -49.51663875 -49.51672804 0.00017723
             1 7 20.0 DavidsonDiagonalization() 9.51e-7 1.0 0.3 7.7 nothing nothing nothing
            ],
            [:n, :i, :ecut, :diag, :ethr, :avg, :β, :t, :ε, :hf, :δ],
        ),
        :n,
    )

    @test parse_electrons_energies(str, :converged) == (-49.51665947,
        -49.51665955,
        7.1e-7,
        nothing,
        [-320.16789804, 160.63877064, -39.2194214, 149.23211567],
        nothing,
        -0.00022634,)

    @test parse_version(str) == "6.1"

    @test parse_parallel_info(str) == ("Parallel version (MPI & OpenMP)", 1)

    @test parse_fft_dimensions(str) == (32157, (nr1 = 32, nr2 = 32, nr3 = 72))

    @test parse_bands(str) == ([
         0.0833 0.25 0.4167
         0.0833 0.0833 0.0833
         0.0 0.0 0.0
         0.0833 0.25 0.4167
         0.25 0.25 0.25
         0.0 0.0 0.0
         0.0833 0.25 0.4167
         0.4167 0.4167 0.4167
         0.0 0.0 0.0
         0.0833 0.25 0.4167
         -0.4167 -0.4167 -0.4167
         0.0 0.0 0.0
         0.0833 0.25 0.4167
         -0.25 -0.25 -0.25
         0.0 0.0 0.0
         0.0833 0.25 0.4167
         -0.0833 -0.0833 -0.0833
         0.0 0.0 0.0
        ],
        [
         -11.3527 -3.0064 -4.6019 -6.0688 -8.6654 -10.3632 -2.7725 -4.6019 -3.8942 -8.179
         -7.7771 -1.8457 -3.9029 -5.3121 -6.5108 -9.6289 -2.7724 -3.9029 -3.1982 -7.4569
         -7.7771 -10.6103 -2.6623 -4.2914 -4.837 -8.179 -11.1046 -2.6623 -3.1982 -4.3652
         -6.0688 -9.8749 -1.3875 -3.4185 -4.2753 -7.4569 -8.9098 -1.3875 -2.5046 -3.8449
         -6.0688 -7.0473 -11.1045 -2.1246 -4.2752 -4.3652 -7.533 -10.3632 -2.1147 -3.6672
         -5.3121 -6.3276 -8.9098 -2.0766 -2.9617 -3.8448 -5.8411 -9.6289 -2.0504 -3.2626
         -4.2914 -5.3896 -7.533 -11.1046 -2.7727 -3.6672 -5.3904 -8.179 -9.8711 -2.4627
         -3.4185 -4.7336 -5.8411 -8.9098 -2.7726 -3.2626 -5.0744 -7.4569 -9.1394 -2.3068
         -2.1245 -4.6019 -5.3904 -7.533 -10.3632 -2.4627 -4.4973 -4.3652 -9.1394 -10.6103
         -2.0766 -3.9029 -5.0744 -5.8411 -9.6289 -2.3067 -3.9372 -3.8448 -8.4114 -9.8749
         -11.1046 -2.6623 -4.4973 -5.3904 -8.179 -10.8568 -3.0064 -3.6672 -3.8942 -7.0473
         -8.9098 -1.3875 -3.9372 -5.0744 -7.4569 -8.6654 -1.8457 -3.2623 -3.1982 -6.3276
         -7.533 -10.6103 -3.0064 -4.4973 -4.3652 -8.6654 -10.6103 -2.4624 -3.1982 -5.3896
         -5.8411 -9.8749 -1.8457 -3.9372 -3.8449 -6.5108 -9.8749 -2.3069 -2.5042 -4.7336
         -5.3904 -7.0473 -11.3527 -3.0064 -3.6672 -4.837 -7.0473 -9.8711 -2.1147 -4.6019
         -5.0744 -6.3276 -7.7771 -1.8457 -3.2623 -4.2752 -6.3276 -9.1394 -2.0499 -3.9029
         -4.4973 -5.3896 -7.7771 -10.8568 -2.4628 -4.2752 -5.3896 -9.1394 -10.3632 -2.6623
         -3.9372 -4.7337 -6.0688 -8.6654 -2.3069 -2.9617 -4.7336 -8.4114 -9.6289 -1.3875
        ],)

    @test parse_clock(str) == DataFrame([
         "" "init_run" 0.47 0.54 1
         "" "electrons" 6.26 7.18 1
         "" "forces" 0.12 0.22 1
         "init_run" "wfcinit" 0.32 0.37 1
         "init_run" "potinit" 0.05 0.06 1
         "electrons" "c_bands" 5.48 6.2 8
         "electrons" "sum_band" 0.62 0.77 8
         "electrons" "v_of_rho" 0.13 0.18 8
         "electrons" "newd" 0.03 0.04 8
         "electrons" "mix_rho" 0.01 0.02 8
         "c_bands" "init_us_2" 0.07 0.1 324
         "c_bands" "cegterg" 5.38 6.08 144
         "sum_band" "sum_band:bec" 0.0 0.01 144
         "sum_band" "addusdens" 0.04 0.04 8
         "*egterg" "h_psi" 4.27 5.19 835
         "*egterg" "s_psi" 0.14 0.08 835
         "*egterg" "g_psi" 0.03 0.03 673
         "*egterg" "cdiaghg" 0.22 0.23 799
         "h_psi" "h_psi:pot" 4.24 5.15 835
         "h_psi" "h_psi:calbec" 0.16 0.16 835
         "h_psi" "vloc_psi" 3.93 4.87 835
         "h_psi" "add_vuspsi" 0.14 0.11 835
         "General routines" "calbec" 0.22 0.22 1051
         "General routines" "fft" 0.09 0.1 119
         "General routines" "fftw" 4.03 4.93 12632
         "General routines" "davcio" 0.0 0.0 18
         "Parallel routines" "fft_scatter" 1.32 1.38 12751
        ],
        [:subroutine, :item, :CPU, :wall, :calls],
    )

    @test whatinput(str) == "standard input"

    @test isrelaxed(str) == true

    @test isjobdone(str) == true
end
