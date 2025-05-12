import numpy as np
import matplotlib.pyplot as plt

data = np.array([
    [7, 2211, 228480, 291.24,   0.128446, 0.176346],
    [6, 2123, 253344, 206.039,  0.0954487, 0.284039],
    [5, 1719, 251328, 121.175,  0.0704912, 0.430391],
    [4, 1833, 258300, 103.59,   0.0565022, 0.649805],
    [3, 1238, 256880, 53.4068,  0.0431394, 1.1016],
    [2, 877,  242640, 30.5197,  0.0347999, 1.98713],
    [1, 580,  244288, 23.5373,  0.0405813, 3.52152]
])

plot_GFLOPs = 0
plot_Tmoy = 0
plot_Ttot = 1



mesh = data[:,0]
Nsteps = data[:,1]
DOF    = data[:,2]
Ttot   = data[:,3]
Tmoy   = data[:,4]
GFLP   = data[:,5]

## Plot des GFLOP/s en fonction du degré des polynômes
if plot_GFLOPs:
    # Version anglaise
    fig = plt.figure()
    plt.plot(mesh, GFLP, "--o" )
    plt.title("Computational Performance as a function of the Polynomial Degree")
    plt.xlabel("Polynomial Degree")
    plt.ylabel("Computational Performance (GFLOP/s)")
    plt.grid()
    plt.show()


    # Version française
    fig = plt.figure()
    plt.plot(mesh, GFLP, "--o" )
    plt.title("Performance Arithmétique en fonction du Degré des Polynômes")
    plt.xlabel("Degré des polynômes")
    plt.ylabel("Performance Arithmétique (GFLOP/s)")
    plt.grid()
    plt.savefig("D:/valen/Documents/Projects/python/Projet_SIM3_Graphes/perf1_GFLOP_fr.pdf")
    plt.show()


## Plot du temps d'exécution moyen par itération en fonction du degré des polynômes
if plot_Tmoy:
    # Version française
    fig = plt.figure()
    plt.plot(mesh, Tmoy, "--o" )
    plt.title("Temps d'exécution moyen par itération en fonction du degré des polynômes")
    plt.xlabel("Degré des polynômes")
    plt.ylabel("Temps d'exécution moyen par itération (s)")
    plt.grid()
    plt.savefig("D:/valen/Documents/Projects/python/Projet_SIM3_Graphes/perf1_tempsmoy_fr.pdf")
    plt.show()


if plot_Ttot:
    # Version française
    fig = plt.figure()
    plt.plot(mesh, Ttot, "--o")
    plt.title("Temps d'exécution total du code en fonction du degré des polynômes")
    plt.xlabel("Degré des polynômes")
    plt.ylabel("Temps d'exécution total (s)")
    plt.grid()
    plt.savefig("D:/valen/Documents/Projects/python/Projet_SIM3_Graphes/perf1_tempstot_fr.pdf")
    plt.show()