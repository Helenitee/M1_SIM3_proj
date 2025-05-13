using Plots
#       mesh Nsteps DOF Ttot Tmoy GFLP
data = [7 2211 228480 291.24  0.128446 0.176346;
        6 2123 253344 206.039 0.0954487 0.284039;
        5 1719 251328 121.175 0.0704912 0.430391;
        4 1833 258300 103.59  0.0565022 0.649805;
        3 1238 256880 53.4068 0.0431394 1.1016;
        2 877  242640 30.5197 0.0347999 1.98713;
        1 580  244288 23.5373 0.0405813 3.52152]

mesh   = @view data[:,1]
Nsteps = @view data[:,2]
DOF    = @view data[:,3]
Ttot   = @view data[:,4]
Tmoy   = @view data[:,5]
GFLP   = @view data[:,6]

begin
    plot(mesh,GFLP,
        title = "Arithmetic Performance as a function of the Polynomial Degree",
        titlefont=font(10),
        label = "GFLOP", 
        marker = :point,
        lw = 2,
        yaxis = ("GFLOP/s", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("Desktop/SIM3/Projet/perf1GFLOP.pdf")

begin
    plot(mesh,Ttot,
        title = "Total runtime as a function of the Polynomial Degree",
        titlefont=font(10),
        marker = :point,
        lw = 2,
        yaxis = ("Time", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("Desktop/SIM3/Projet/perf1Ttot.pdf")
begin
    plot(mesh,Tmoy, 
        title = "Average runtime per iteration as a function of the Polynomial Degree",
        titlefont=font(10),
        label = "",
        marker = :point,
        lw = 2,
        yaxis = ("Time (s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("Desktop/SIM3/Projet/perf1Tmoy.pdf")