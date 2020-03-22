from No1_a import RK4, graph, u_pt

def v_pt(u, v, t):
    return -(u**3.0)

RK4_c1 = RK4(u_pt, v_pt, [1.0, 0.0], [0.0, 50.0])
RK4_c10 = RK4(u_pt, v_pt, [10.0, 0.0], [0.0, 50.0])
RK4_c50 = RK4(u_pt, v_pt, [50.0, 0.0], [0.0, 50.0])

graph(RK4_c1[0], RK4_c1[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=1")
graph(RK4_c1[0], RK4_c1[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=10")
graph(RK4_c1[0], RK4_c1[2],"Graphique de la position en fonction du temps d'un oscillateur anharmonique en 1D avec x(0)=50")