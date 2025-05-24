using Plots

datanoPARA = [7 2211 228480 306.111 0.138449 6.69808;
            6 2123 253344 220.112 0.103679 7.08522;
            5 1719 251328 130.253 0.0757726 6.60264;
            4 1833 258300 109.709 0.0598523 5.61878;
            3 1238 256880 55.5406 0.044863 4.60215;
            2 877 242640 31.2653 0.0356501 3.17335;
            1 580 244288 23.4482 0.0404278 1.57107]

dataPARA = [7 2211 228480 95.7693 0.0433148 21.4093;
            6 2123 253344 69.4934 0.0327334 22.4416;
            5 1719 251328 40.9982 0.0238499 20.9769;
            4 1833 258300 34.4877 0.0188147 17.874;
            3 1238 256880 17.1585 0.0138597 14.8968;
            2 877 242640 9.65181 0.0110053 10.2795;
            1 580 244288 7.34264 0.0126595 5.01708]

mesh   = @view datanoPARA[:,1]
Nsteps = @view datanoPARA[:,2]
DOF    = @view datanoPARA[:,3]
Ttot   = @view datanoPARA[:,4]
Tmoy   = @view datanoPARA[:,5]
GFLP   = @view datanoPARA[:,6]

Ttotpara   = @view dataPARA[:,4]
Tmoypara   = @view dataPARA[:,5]
GFLPpara   = @view dataPARA[:,6]


begin
    plot(mesh,GFLP,
        title = " ",
        titlefont=font(9),
        label = "Non parallélisé", 
        marker = :point,
        lw = 2,
        yaxis = ("Performance arithmétique (GFLOP/s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, GFLPpara,
        label = "Parallèle sur k",
        marker = :point,
        lw = 2)
end
savefig("figures/perf4GFLOP_comp.pdf")

begin
    plot(mesh,Ttot,
        title = " ",
        titlefont=font(9),
        label = "Non parallélisé",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Ttotpara,
        label = "Parallèle sur k",
        marker = :point
    )
end
savefig("figures/perf4tot_comp.pdf")
begin
    plot(mesh,Tmoy, 
        title = " ",
        titlefont=font(9),
        label = "Non parallélisé",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Tmoypara,
        label = "Parallèle sur k",
        marker = :point)
end
savefig("figures/perf4Tmoy_comp.pdf")
