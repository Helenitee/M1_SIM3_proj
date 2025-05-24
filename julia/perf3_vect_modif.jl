using Plots
#     mesh Nsteps DOF     Ttot    Tmoy    GFLP
data1 = [1 580 244288 24.099809 0.040577 0.001286;
         2 877 242640 33.429988 0.037335 0.001109;
         3 1238 256880 56.977587 0.045253 0.000857;
         4 1833 258300 108.190225 0.058222 0.000569;
         5 1719 251328 124.823212 0.071797 0.000598;
         6 2123 253344 208.069190 0.097052 0.000445;
         7 2211 228480 290.004796 0.130101 0.000351
] # pre-computed coefs 

data2 = [1 580 244288 23.635822 0.039781 0.001311;
         2 877 242640 31.820934 0.035529 0.001165;
         3 1238 256880 55.900908 0.044417 0.000874;
         4 1833 258300 112.727663 0.060695 0.000546;
         5 1719 251328 133.581571 0.076892 0.000559;
         6 2123 253344 224.076248 0.104618 0.000414;
         7 2211 228480 308.919124 0.138635 0.000329
] # pre-computed coefs + vectors s_p, ... not used!

mesh   = @view data1[:, 1]
Ttot1  = @view data1[:, 4]
Tmoy1  = @view data1[:, 5]
gflop1 = @view data1[:, 6]
Ttot2  = @view data2[:, 4]
Tmoy2  = @view data2[:, 5]
gflop2 = @view data2[:, 6]


begin
    plot(mesh,gflop2,
        # title = "Comparison of the arithmetic performance as functions of the polynomial degree",
        # titlefont=font(9),
        label = "vecteurs s_p, s_u, s_v, s_w non utilisés",
        marker = :point,
        lw = 2,
        yaxis = ("GFLOP/s", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, gflop1, label = "coefficients pre-calculés", marker = :point, lw = 2)
end
savefig("figures/perf3_comp_vect_gflop.pdf")

begin
    plot(mesh,Ttot2,
        # title = "Comparison of the arithmetic performance as functions of the polynomial degree",
        # titlefont=font(9),
        label = "vecteurs s_p, s_u, s_v, s_w non utilisés",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Ttot1, label = "coefficients pre-calculés", marker = :point, lw = 2)
end
savefig("figures/perf3_comp_vect_Ttot_.pdf")

begin
    plot(mesh,Tmoy2,
        # title = "Comparison of the arithmetic performance as functions of the polynomial degree",
        # titlefont=font(9),
        label = "vecteurs s_p, s_u, s_v, s_w non utilisés",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré Polynomial", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Tmoy1, label = "coefficients pre-calculés", marker = :point, lw = 2)
end
savefig("figures/perf3_comp_vect_Tmoy.pdf")
