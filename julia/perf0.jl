using Plots
#       mesh Nsteps DOF Ttot Tmoy GFLP
data = [1 580 244288 24.339927 0.040992 0.001273;
        2 877 242640 31.513164 0.035178 0.001176;
        3 1238 256880 54.110842 0.042969 0.000903;
        4 1833 258300 107.849607 0.057982 0.000570;
        5 1719 251328 124.515850 0.071614 0.000599;
        6 2123 253344 210.613216 0.098274 0.000440;
        7 2211 228480 291.822604 0.130924 0.000348
]

mesh   = @view data[:,1]
Nsteps = @view data[:,2]
DOF    = @view data[:,3]
Ttot   = @view data[:,4]
Tmoy   = @view data[:,5]
GFLP   = @view data[:,6]


begin
    plot(mesh,GFLP,
        # title = "Arithmetic Performances as functions of the Polynomial Degree",
        titlefont=font(10),
        label = "", 
        marker = :point,
        lw = 2,
        yaxis = ("GFLOP/s", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("figures/perf0_glop.pdf")

begin
    plot(mesh,Ttot,
        # title = "Total runtime as a function of the Polynomial Degree",
        titlefont=font(10),
        label = "",
        marker = :point,
        lw = 2,
        yaxis = ("Time (s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("figures/perf0_Ttot.pdf")

begin
    plot(mesh,Tmoy, 
        # title = "Average runtime per iteration as a function of the Polynomial Degree",
        titlefont=font(10),
        label = "",
        marker = :point,
        lw = 2,
        yaxis = ("Time (s)", font(10)),
        xaxis = ("Polynomial Degree", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
end
savefig("figures/perf0_Tmoy.pdf")
