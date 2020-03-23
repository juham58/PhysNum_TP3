from scipy.integrate import solve_ivp
import scipy.constants as cst
from No1_d import graph
from No1_e import vdp, rangetemps

if __name__ == "__main__":
    trange = [0, 8 * cst.pi]
    t = rangetemps(1000000, trange)

    solution1 = solve_ivp(vdp, [0, 8 * cst.pi], [1.0, 0], t_eval=t, method='RK45')
    sol1 = solution1.y

    solution2 = solve_ivp(vdp, [0, 8 * cst.pi], [2.0, 0], t_eval=t, method='RK45')
    sol2 = solution2.y

    solution3 = solve_ivp(vdp, [0, 8 * cst.pi], [3.0, 0], t_eval=t, method='RK45')
    sol3 = solution3.y

#Par curiosité, nous avons même testé avec une vitesse initiale v(0) = 3:
    solution4 = solve_ivp(vdp, [0, 8 * cst.pi], [3.0, 3.0], t_eval=t, method='RK45')
    sol4 = solution4.y

#graphique avec x(0) = 1, 2, 3
    graph(sol1[0], sol1[1], "Graphique de l'espace des phases du système de van der Pol avec x(0) = 1.0 et v(0) = 0.0")
    graph(sol2[0], sol2[1], "Graphique de l'espace des phases du système de van der Pol avec x(0) = 2.0 et v(0) = 0.0")
    graph(sol3[0], sol3[1], "Graphique de l'espace des phases du système de van der Pol avec x(0) = 3.0 et v(0) = 0.0")
#graphique d'une vitesse initiale non-nulle
    graph(sol4[0], sol4[1], "Graphique de l'espace des phases du système de van der Pol avec x(0) = 3.0 et v(0) = 3.0")