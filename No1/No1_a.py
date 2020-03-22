import numpy as np
import matplotlib.pyplot as plt

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
        k3v = h*fonct_v(u+0.5*k2v, v, t+0.5*h)
        k4v = h*fonct_v(u +k3v, v, t+h)

        u += (k1u +2*k2u +2*k3u +k4u)/6
        v += (k1v +2*k2v +2*k3v +k4v)/6

    return upoints, vpoints, tpoints

def graph(axe_y, axe_x, titre):
    plt.figure(figsize=(16,8))
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.title(titre)
    plt.grid()
    plt.plot(axe_x, axe_y, "k-")
    plt.show()

def u_pt(u, v, t):
    return 1.0*v

def v_pt(u, v, t):
    return - 1.0*u

if __name__ == "__main__":
    RK4_a = RK4(u_pt, v_pt, [1.0, 0.0], [0.0, 50.0])
    graph(RK4_a[0], RK4_a[2], "Graphique de la position en fonction du temps d'un oscillateur harmonique en 1D")