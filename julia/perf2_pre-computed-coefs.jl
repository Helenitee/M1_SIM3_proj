using Plots
#     mesh Nsteps DOF     Ttot    Tmoy    GFLP
data1 = [1 580 244288 24.339927 0.040992 0.001273;
         2 877 242640 31.513164 0.035178 0.001176;
         3 1238 256880 54.110842 0.042969 0.000903;
         4 1833 258300 107.849607 0.057982 0.000570;
         5 1719 251328 124.515850 0.071614 0.000599;
         6 2123 253344 210.613216 0.098274 0.000440;
         7 2211 228480 291.822604 0.130924 0.000348
] # non pre computed coefs

data2 = [1 580 244288 24.099809 0.040577 0.001286;
         2 877 242640 33.429988 0.037335 0.001109;
         3 1238 256880 56.977587 0.045253 0.000857;
         4 1833 258300 108.190225 0.058222 0.000569;
         5 1719 251328 124.823212 0.071797 0.000598;
         6 2123 253344 208.069190 0.097052 0.000445;
         7 2211 228480 290.004796 0.130101 0.000351
] # pre-computed coefs

mesh  = @view data1[:,1]
Ttot  = @view data1[:,4]
Tmoy  = @view data1[:,5]
gflop  = @view data1[:,6]
Ttot2 = @view data2[:,4]
Tmoy2 = @view data2[:,5]
gflop2 = @view data2[:,6]

begin
    plot(mesh,gflop,
        # title = "Comparison of the arithmetic performance as functions of the polynomial degree",
        # titlefont=font(9),
        label = "Non optimisé", 
        marker = :point,
        lw = 2,
        yaxis = ("GFLOP/s", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, gflop2,
        label = "coefficients pre-calculés",
        marker = :point,
        lw = 2)
end
savefig("figures/perf2_comp_pre-computed_gflop.pdf")

begin
    plot(mesh,Ttot,
        # title = "Comparison of the total runtime as a function of the polynomial degree",
        # titlefont=font(9),
        label = "Non optimisé", 
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Ttot2,
        label = "coefficients pre-calculés",
        marker = :point, lw=2
    )
end
savefig("figures/perf2_comp_pre-computed_Ttot_.pdf")

begin
    plot(mesh,Tmoy, 
        # title = "Average runtime per iteration as a function of the Polynomial Degree",
        # titlefont=font(9),
        label = "Non optimisé",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Tmoy2,
        label = "coefficients pre-calculés",
        marker = :point, lw=2)
end
savefig("figures/perf2_comp_pre-computed_Tmoy.pdf")
