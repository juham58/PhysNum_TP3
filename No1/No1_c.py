from No1_a import RK4, graph, u_pt

def v_pt(u, v, t):
    return -(u**3.0)

def periode(fonction, val_init):
    liste1 = fonction[0]
    liste2 = fonction[2]
    for i in range(len(liste1)):
        if liste1[i] > val_init:
            return liste2[i]


if __name__ == "__main__":
    RK4_c1 = RK4(u_pt, v_pt, [1.0, 0.0], [0.0, 50.0], 10000)
    RK4_c3 = RK4(u_pt, v_pt, [3.0, 0.0], [0.0, 50.0], 100000)
    RK4_c5 = RK4(u_pt, v_pt, [5.0, 0.0], [0.0, 50.0], 100000)

    print("La periode pour x(0) = 1 est de " + str(periode(RK4_c1, 1.0)) + " secondes.")
    print("La periode pour x(0) = 3 est de " + str(periode(RK4_c3, 3.0)) + " secondes.")
    print("La periode pour x(0) = 5 est de " + str(periode(RK4_c5, 5.0)) + " secondes.")

    graph(RK4_c1[0], RK4_c1[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=1")
    graph(RK4_c3[0], RK4_c3[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=3")
    graph(RK4_c5[0], RK4_c5[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=5")