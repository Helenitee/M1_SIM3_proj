using Plots

dataNOOPTI = [7 2211 228480 294.276 0.13203 6.96745;
        6 2123 253344 209.753 0.0978694 7.43513;
        5 1719 251328 129.139 0.0743028 6.65963;
        4 1833 258300 107.567 0.0578833 5.7307;
        3 1238 256880 54.1156 0.0429734 4.72334;
        2 877 242640 31.2516 0.0348786 3.17474;
        1 580 244288 24.0365 0.040471 1.53261]


dataOPTI = [7 2211 228480 95.7693 0.0433148 21.4093;
            6 2123 253344 69.4934 0.0327334 22.4416;
            5 1719 251328 40.9982 0.0238499 20.9769;
            4 1833 258300 34.4877 0.0188147 17.874;
            3 1238 256880 17.1585 0.0138597 14.8968;
            2 877 242640 9.65181 0.0110053 10.2795;
            1 580 244288 7.34264 0.0126595 5.01708]

dataBLASv1 = [7 2211 228480 256.064 0.115813 8.00721;
            6 2123 253344 275.654 0.129842 5.65759;
            5 1719 251328 220.966 0.128543 3.89208;
            4 1833 258300 231.668 0.126387 2.66085;
            3 1238 256880 168.386 0.136014 1.51798;
            2 877 242640 117.391 0.133855 0.845175;
            1 580 244288 88.557 0.152684 0.415988]


dataBLASv2 = [7 2211 228480 61.7357 0.0279219 33.2118;
            6 2123 253344 53.5882 0.0252416 29.1023;
            5 1719 251328 34.5013 0.0200704 24.927;
            4 1833 258300 36.7017 0.0200226 16.7957;
            3 1238 256880 25.9523 0.0209629 9.8491;
            2 877 242640 22.5621 0.0257263 4.39746;
            1 580 244288 29.4599 0.0507926 1.25047]

mesh   = @view dataNOOPTI[:,1]
Nsteps = @view dataNOOPTI[:,2]
DOF    = @view dataNOOPTI[:,3]
Ttot   = @view dataNOOPTI[:,4]
Tmoy   = @view dataNOOPTI[:,5]
GFLP   = @view dataNOOPTI[:,6]

TtotOPTI   = @view dataOPTI[:,4]
TmoyOPTI   = @view dataOPTI[:,5]
GFLPOPTI   = @view dataOPTI[:,6]

TtotBLAS = @view dataBLASv1[:,4]
TmoyBLAS = @view dataBLASv1[:,5]
GFLPBLAS = @view dataBLASv1[:,6]


TtotBLAS2 = @view dataBLASv2[:,4]
TmoyBLAS2 = @view dataBLASv2[:,5]
GFLPBLAS2 = @view dataBLASv2[:,6]

begin
    plot(mesh,GFLP,
        title = " ",
        titlefont=font(9),
        label = "Original", 
        marker = :point,
        lw = 2,
        yaxis = ("Performance arithmétique (GFLOP/s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, GFLPOPTI,
        label = "Optimisé",
        marker = :point,
        lw = 2)
    plot!(mesh, GFLPBLAS,
        label = "BLAS 1",
        marker = :point,
        lw = 2)
    plot!(mesh, GFLPBLAS2,
        label = "BLAS 2",
        marker = :point,
        lw = 2)
end
savefig("figures/perf5GFLOP_comp.pdf")

begin
    plot(mesh,Ttot,
        title = " ",
        titlefont=font(9),
        label = "Original",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, TtotOPTI,
        label = "Optimisé",
        marker = :point,
        lw = 2)
    plot!(mesh, TtotBLAS,
        label = "BLAS 1",
        marker =:point,
        lw = 2)
    plot!(mesh, TtotBLAS2,
        label = "BLAS 2",
        marker =:point,
        lw = 2)
end
savefig("figures/perf5Ttot_comp.pdf")
begin
    plot(mesh,Tmoy, 
        title = " ",
        titlefont=font(9),
        label = "Original",
        marker = :point,
        lw = 2,
        yaxis = ("Temps (s)", font(10)),
        xaxis = ("Degré des fonctions de base", font(10)),
        xticks=(1:1:7),
        minorgrid = true)
    plot!(mesh, TmoyOPTI,
        label = "Optimisé",
        marker = :point,
        lw = 2)
    plot!(mesh, TmoyBLAS,
        label = "BLAS 1",
        marker = :point,
        lw = 2)
    plot!(mesh, TmoyBLAS2,
        label = "BLAS 2",
        marker = :point,
        lw = 2)
end
savefig("figures/perf5Tmoy_comp.pdf")
