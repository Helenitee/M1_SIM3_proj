using Plots
#       mesh Nsteps DOF Ttot Tmoy GFLP
data = [7 2211 228480 291.24  0.128446 0.176346;
        6 2123 253344 206.039 0.0954487 0.284039;
        5 1719 251328 121.175 0.0704912 0.430391;
        4 1833 258300 103.59  0.0565022 0.649805;
        3 1238 256880 53.4068 0.0431394 1.1016;
        2 877  242640 30.5197 0.0347999 1.98713;
        1 580  244288 23.5373 0.0405813 3.52152]


data2 = [7 2211 228480 285.653  0.129196 0.179795;
        6 2123 253344 204.943 0.0965345 0.285558;
        5 1719 251328 121.693 0.0707926 0.428559;
        4 1833 258300 105.937  0.0577944 0.635276;
        3 1238 256880 54.3709 0.0439182  1.08207;
        2 877  242640 31.1017 0.0354636  1.94995;
        1 580  244288 23.616 0.0407171 3.50978]

mesh   = @view data[:,1]
Nsteps = @view data[:,2]
DOF    = @view data[:,3]
Ttot   = @view data[:,4]
Tmoy   = @view data[:,5]
GFLP   = @view data[:,6]

Ttot2   = @view data2[:,4]
Tmoy2   = @view data2[:,5]
GFLP2   = @view data2[:,6]


begin
    plot(mesh,GFLP,
        title = "Comparison of the arithmetic performance as functions of the polynomial degree",
        titlefont=font(9),
        label = "Unoptimized", 
        marker = :point,
        lw = 2,
        yaxis = ("Arithmetic Performance (GFLOP/s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, GFLP2,
        label = "Pre-calculated",
        marker = :point,
        lw = 2)
end
savefig("figures/perf2GFLOP_comparaison.pdf")

begin
    plot(mesh,Ttot,
        title = "Comparison of the total runtime as a function of the polynomial degree",
        titlefont=font(9),
        label = "Unoptimized",
        marker = :point,
        lw = 2,
        yaxis = ("Time (s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Ttot2,
        label = "Pre-calculated",
        marker = :point
    )
end
savefig("figures/perf2tot_comparaison.pdf")
begin
    plot(mesh,Tmoy, 
        title = "Average runtime per iteration as a function of the Polynomial Degree",
        titlefont=font(9),
        label = "Unoptimized",
        marker = :point,
        lw = 2,
        yaxis = ("Time (s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Tmoy2,
        label = "Pre-calculated",
        marker = :point)
end
savefig("figures/perf2Tmoy_comparaison.pdf")
