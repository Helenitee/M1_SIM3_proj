using Plots 

# 1st line of each new "mesh" is the non-optimized compilation
# 2sh line of each new "mesh" is the 1. parallelized before each for loop
#                                and 2. definition of variables computes multiple times
# BUT we didn't presice the nb of threads we were using!
data = [
#    mesh Nsteps DOF  Ttot    Tmoy      GFLP
        1 580 244288 23.654141 0.039810 0.001310;
        1 580 244288 8.697748 0.013957 0.003563;
        1 580 244288 7.875630 0.012540 0.003935;
        2 877 242640 31.708632 0.035399 0.001169;
        2 877 242640 12.056385 0.012935 0.003074;
        2 877 242640 10.403350 0.011051 0.003563;
        3 1238 256880 56.144169 0.044611 0.000870;
        3 1238 256880 20.456454 0.015730 0.002387;
        3 1238 256880 18.093384 0.013822 0.002699;
        4 1833 258300 112.240059 0.060395 0.000548;
        4 1833 258300 38.750534 0.020338 0.001587;
        4 1833 258300 35.598210 0.018619 0.001728;
        5 1719 251328 133.176513 0.076642 0.000560;
        5 1719 251328 45.144566 0.025443 0.001653;
        5 1719 251328 42.343810 0.023818 0.001762;
        6 2123 253344 227.536071 0.106247 0.000407;
        6 2123 253344 75.665758 0.034706 0.001225;
        6 2123 253344 71.591680 0.032790 0.001295;
        7 2211 228480 309.429414 0.138886 0.000329;
        7 2211 228480 100.699086 0.044464 0.001010;
        7 2211 228480 97.519817 0.043039 0.001043;
        ]

mesh           = @view data[1:3:end  ,1]
Tmoy_no_paral  = @view data[1:3:end-2,5]
Tmoy_paral_1K  = @view data[2:3:end-1,5]
Tmoy_paral_2K  = @view data[3:3:end  ,5]
Ttot_no_paral  = @view data[1:3:end-2,5]
Ttot_paral_1K  = @view data[2:3:end-1,4]
Ttot_paral_2K  = @view data[3:3:end  ,4]
gflop_no_paral = @view data[1:3:end-2,5]
gflop_paral_1K = @view data[1:3:end-1,6]
gflop_paral_2K = @view data[2:3:end  ,6]

# Tmoy
begin
    plot(mesh, Tmoy_no_paral,
        label = "Tmoy no paral", 
        marker = :point,
        lw = 2,
        yaxis = ("time", font(10)),
        xaxis = ("polyomial degree", font(10)),
        xticks=(1:1:7),
        # ytricks = log2,
        minorgrid = true)
    plot!(mesh, Tmoy_paral_1K, label = "Tmoy paral in the 1st K-loop", lw = 2, marker = :point)
    plot!(mesh, Tmoy_paral_2K, label = "Tmoy paral in both K-loops", lw = 2, marker = :point)
end
savefig("figures/perf1_comp_paral_Tmoy.pdf")

#Ttot
begin
    plot(mesh, Ttot_no_paral,
        label = "Ttot no paral", 
        marker = :point,
        lw = 2,
        yaxis = ("time", font(10)),
        xaxis = ("polyomial degree", font(10)),
        xticks=(1:1:7),
        # ytricks = log2,
        minorgrid = true)
    plot!(mesh, Ttot_paral_1K, label = "Ttot paral in the 1st K-loop", lw = 2, marker = :point)
    plot!(mesh, Ttot_paral_2K, label = "Ttot paral in both K-loops", lw = 2, marker = :point)
end
savefig("figures/perf1_comp_paral_Ttot.pdf")

#GFLOP/s
begin
    plot(mesh, gflop_no_paral,
        label = "GFLOP/s no paral", 
        marker = :point,
        lw = 2,
        yaxis = ("time", font(10)),
        xaxis = ("polyomial degree", font(10)),
        xticks=(1:1:7),
        # ytricks = log2,
        minorgrid = true)
    plot!(mesh, gflop_paral_1K, label = "GFLOP/s paral in the 1st K-loop", lw = 2, marker = :point)
    plot!(mesh, gflop_paral_2K, label = "GFLOP/s paral in both K-loops", lw = 2, marker = :point)
end
savefig("figures/perf1_comp_paral_gflops.pdf")