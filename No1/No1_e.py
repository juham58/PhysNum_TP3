import numpy as np
from scipy.integrate import solve_ivp
import scipy.constants as cst
from No1_d import graph

#vdp avec epsilon=1 déjà dans l'équation de x3
def vdp(t,x):
    x1 = x[0]
    x2 = x[1]
    x3 = -x1 - (x1 ** 2 - 1) * x2
    return np.array([x2, x3], float)

#fonction qui fait la liste de points d'intervalle h entre t=0 et t=8*pi
def rangetemps(N, trange):
   return np.arange(trange[0], trange[1], (trange[1] - trange[0]) / N)

#on pose le range, on trouve notre liste de temps et on solutionne pour obtenir notre graphique
if __name__ == "__main__":
    trange = [0, 8 * cst.pi]
    t = rangetemps(1000000, trange)

    solution = solve_ivp(vdp, [0, 8 * cst.pi], [0.5, 0], t_eval=t, method='RK45')
    sol = solution.y

    graph(sol[0], sol[1], "Graphique de l'espace des phases du système de van der Pol avec x(0) = 0.5 et v(0) = 0.0")