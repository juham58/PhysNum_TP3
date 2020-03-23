import matplotlib.pyplot as plt
from No1_c import RK4, v_pt
from No1_a import u_pt

def graph(axe_y, axe_x, titre):
    plt.figure(figsize=(16,8))
    plt.xlabel("x(t)")
    plt.ylabel("v(t)")
    plt.title(titre)
    plt.grid()
    plt.plot(axe_x, axe_y, "k-")
    plt.show()

if __name__ == "__main__":
    RK4_c1 = RK4(u_pt, v_pt, [1.0, 0.0], [0.0, 50.0], 10000)
    RK4_c5 = RK4(u_pt, v_pt, [5.0, 0.0], [0.0, 50.0], 100000)

    #Nous prenons RK4 pour x(0) = 1 et RK4 pour x(0) = 5
    graph(RK4_c1[1], RK4_c1[0], "Graphique de l'espace de phase d'un oscillateur anharmonique avec x(0) = 1")
    graph(RK4_c5[1], RK4_c5[0], "Graphique de l'espace de phase d'un oscillateur anharmonique avec x(0) = 5")