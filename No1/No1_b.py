from No1_a import RK4, graph, u_pt, v_pt

#On change les conditions initiales dans les parametres de RK4
RK4_b = RK4(u_pt, v_pt, [5.0, 10.0], [0.0, 50.0])
graph(RK4_b[0], RK4_b[2], "Graphique de la position en fonction du temps d'un oscillateur harmonique en 1D avec x(0)=5 et v(0)=10")