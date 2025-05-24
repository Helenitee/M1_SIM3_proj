using Plots
#       mesh Nsteps DOF Ttot Tmoy GFLP
dataLOOPnm = [7 2211 228480 314.625 0.141237 6.51682;
         6 2123 253344 229.891 0.107354 6.78382;
         5 1719 251328 135.742 0.0781168 6.33564;
         4 1833 258300 114.73 0.0617603 5.3729;
         3 1238 256880 57.4848 0.0456947 4.4465;
         2 877 242640 32.7719 0.0366122 3.02746;
         1 580 244288 24.3567 0.0410248 1.51246
] # pre-computed coefs + vectors s_p, ... not used!



dataLOOPSmn = [7 2211 228480 306.111 0.138449 6.69808;
            6 2123 253344 220.112 0.103679 7.08522;
            5 1719 251328 130.253 0.0757726 6.60264;
            4 1833 258300 109.709 0.0598523 5.61878;
            3 1238 256880 55.5406 0.044863 4.60215;
            2 877 242640 31.2653 0.0356501 3.17335;
            1 580 244288 23.4482 0.0404278 1.57107] #LOOP SWAPPED

mesh   = @view dataLOOPnm[:,1]
Nsteps = @view dataLOOPnm[:,2]
DOF    = @view dataLOOPnm[:,3]
Ttot   = @view dataLOOPnm[:,4]
Tmoy   = @view dataLOOPnm[:,5]
GFLP   = @view dataLOOPnm[:,6]

Ttotswap   = @view dataLOOPSmn[:,4]
Tmoyswap   = @view dataLOOPSmn[:,5]
GFLPswap   = @view dataLOOPSmn[:,6]


begin
    plot(mesh,GFLP,
        title = " ",
        titlefont=font(9),
        label = "Non échangées", 
        marker = :point,
        lw = 2,
        yaxis = ("Performance arithmétique (GFLOP/s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, GFLPswap,
        label = "Échangées",
        marker = :point,
        lw = 2)
end
savefig("figures/perf3GFLOP_comp.pdf")

begin
    plot(mesh,Ttot,
        title = " ",
        titlefont=font(9),
        label = "Non échangées",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Ttotswap,
        label = "Échangées",
        marker = :point
    )
end
savefig("figures/perf3tot_comparaison.pdf")
begin
    plot(mesh,Tmoy, 
        title = " ",
        titlefont=font(9),
        label = "Non échangées",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, Tmoyswap,
        label = "Échangées",
        marker = :point)
end
savefig("figures/perf3Tmoy_comparaison.pdf")
