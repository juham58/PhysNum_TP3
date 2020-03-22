from math import sin
from numpy import arange
from pylab import plot,xlabel,ylabel,show

def RK4 (fonct_u, fonct_v, val_init, range):
    t0, t1 = range[0], range[1]
    N = 10000
    h = (t1-t0)/N
    tpoints = np.arange(t0, t1, h)
    upoints, vpoints = [], []
    u, v = val_init[0], val_init[1]
    for t in tpoints:
        upoints.append(u)
        vpoints.append(v)
        k1u = h*fonct_u(u, v, t)
        k2u = h*fonct_u(u, v+0.5*k1u, t+0.5*h)
        k3u = h*fonct_u(u, v+0.5*k2u, t+0.5*h)
        k4u = h*fonct_u(u, v+k3u, t+h)

        k1v = h*fonct_v(u, v, t)
        k2v = h*fonct_v(u+0.5*k1v, v, t+0.5*h)
        k3v = h*fonct_v(x+0.5*k2v, v, t+0.5*h)
        k4v = h*fonct_v(x +k3v, v, t+h)

        u += (k1u +2*k2u +2*k3u +k4u)/6
        v += (k1v +2*k2v +2*k3v +k4v)/6

    return upoints, vpoints, tpoints

def graph(fct, axe_x):
    plt.figure(figsize=(16,8))
    plt.xlabel("Valeur de x")
    plt.ylabel("Resultat de l'equation")
    plt.title("Zeros d'un systeme d'equations lineaires")
    plt.grid()
    plt.plot(axe_x, fct(axe_x), "k-")
    plt.show()

RK4()