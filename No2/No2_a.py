import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# définitions des masses des corps
m_A = 3
m_B = 4
m_C = 5

# définition de la constante gravitationnelle
G = 4*constants.pi**2

# définition des conditions initiales
r_Ai = np.array([1.0, 3.0])
r_Bi = np.array([-2.0, -1.0])
r_Ci = np.array([1.0, -1.0])

v_Ai = np.array([0.0, 0.0])
v_Bi = np.array([0.0, 0.0])
v_Ci = np.array([0.0, 0.0])


def F(corps, r_A, r_B, r_C):
    if corps == "A":
        return -G*(m_B*((r_A-r_B)/(np.linalg.norm(r_A-r_B)**3))
                   + m_C*((r_A-r_C)/(np.linalg.norm(r_A-r_C)**3)))

    if corps == "B":
        return -G*(m_A*((r_B-r_A)/(np.linalg.norm(r_B-r_A)**3))
                   + m_C*((r_B-r_C)/(np.linalg.norm(r_B-r_C)**3)))

    if corps == "C":
        return -G*(m_A*((r_C-r_A)/(np.linalg.norm(r_C-r_A)**3))
                   + m_B*((r_C-r_B)/(np.linalg.norm(r_C-r_B)**3)))


def mouton_3_corps(t_i, t_f, N):
    t_points = np.linspace(t_i, t_f, N)
    rA_arr = np.zeros((len(t_points), 2))
    rB_arr = np.zeros((len(t_points), 2))
    rC_arr = np.zeros((len(t_points), 2))
    h = (t_f-t_i)/N

    # calcul du point r(t+h/2) avec Runge-Kutta d'ordre 4
    k1_A_v = 0.5*h*F("A", r_Ai, r_Bi, r_Ci)
    k1_B_v = 0.5*h*F("B", r_Ai, r_Bi, r_Ci)
    k1_C_v = 0.5*h*F("C", r_Ai, r_Bi, r_Ci)

    k2_A_v = 0.5*h*F("A", r_Ai + 0.5*k1_A_v, r_Bi + 0.5*k1_B_v, r_Ci + 0.5*k1_C_v)
    k2_B_v = 0.5*h*F("B", r_Ai + 0.5*k1_A_v, r_Bi + 0.5*k1_B_v, r_Ci + 0.5*k1_C_v)
    k2_C_v = 0.5*h*F("C", r_Ai + 0.5*k1_A_v, r_Bi + 0.5*k1_B_v, r_Ci + 0.5*k1_C_v)

    k3_A_v = 0.5*h*F("A", r_Ai + 0.5*k2_A_v, r_Bi + 0.5*k2_B_v, r_Ci + 0.5*k2_C_v)
    k3_B_v = 0.5*h*F("B", r_Ai + 0.5*k2_A_v, r_Bi + 0.5*k2_B_v, r_Ci + 0.5*k2_C_v)
    k3_C_v = 0.5*h*F("C", r_Ai + 0.5*k2_A_v, r_Bi + 0.5*k2_B_v, r_Ci + 0.5*k2_C_v)

    k4_A_v = 0.5*h*F("A", r_Ai + 0.5*k3_A_v, r_Bi + 0.5*k3_B_v, r_Ci + 0.5*k3_C_v)
    k4_B_v = 0.5*h*F("B", r_Ai + 0.5*k3_A_v, r_Bi + 0.5*k3_B_v, r_Ci + 0.5*k3_C_v)
    k4_C_v = 0.5*h*F("C", r_Ai + 0.5*k3_A_v, r_Bi + 0.5*k3_B_v, r_Ci + 0.5*k3_C_v)

    v_A_demie = v_Ai + 1/6*(k1_A_v+2*k2_A_v+2*k3_A_v+k4_A_v)
    v_B_demie = v_Bi + 1/6*(k1_B_v+2*k2_B_v+2*k3_B_v+k4_B_v)
    v_C_demie = v_Ci + 1/6*(k1_C_v+2*k2_C_v+2*k3_C_v+k4_C_v)

    r_A_demie = r_Ai + 0.5*h*v_A_demie
    r_B_demie = r_Bi + 0.5*h*v_B_demie
    r_C_demie = r_Ci + 0.5*h*v_C_demie

    v_A = v_Ai
    v_B = v_Bi
    v_C = v_Ci

    r_A = r_Ai
    r_B = r_Bi
    r_C = r_Ci

    rA_arr[0][0], rA_arr[0][1] = r_A[0], r_A[1]
    rB_arr[0][0], rB_arr[0][1] = r_B[0], r_B[1]
    rC_arr[0][0], rC_arr[0][1] = r_C[0], r_C[1]
    # début des calculs par sauts
    for i, t in enumerate(t_points[1:]):


        #print(100*compteur_t/(t_f-t_i), "\n")

        v_A = v_A + h*F("A", r_A_demie, r_B_demie, r_C_demie)
        v_B = v_B + h*F("B", r_A_demie, r_B_demie, r_C_demie)
        v_C = v_C + h*F("C", r_A_demie, r_B_demie, r_C_demie)

        r_A = r_A + h*v_A_demie
        r_B = r_B + h*v_B_demie
        r_C = r_C + h*v_C_demie

        # premier calcul pour un demi-saut
        v_A_demie = v_A_demie + h*F("A", r_A, r_B, r_C)
        v_B_demie = v_B_demie + h*F("B", r_A, r_B, r_C)
        v_C_demie = v_C_demie + h*F("C", r_A, r_B, r_C)

        r_A_demie = r_A_demie + h*v_A
        r_B_demie = r_B_demie + h*v_B
        r_C_demie = r_C_demie + h*v_C

        # deuxième calcul pour compléter le saut
        rA_arr[i+1][0], rA_arr[i+1][1] = r_A[0], r_A[1]
        rB_arr[i+1][0], rB_arr[i+1][1] = r_B[0], r_B[1]
        rC_arr[i+1][0], rC_arr[i+1][1] = r_C[0], r_C[1]

        print((m_A*r_A+m_B*r_B+m_C*r_C)/(m_A + m_B + m_C))  # centre de masse
    return {"A": rA_arr, "B": rB_arr, "C": rC_arr, "t": t_points}


np.set_printoptions(threshold=np.inf)
def graph_3_corps(t_i, t_f, N):
    mouton = mouton_3_corps(t_i, t_f, N)
    #print(mouton["C"])
    fig, ax = plt.subplots()
    ax.set(xlim=(-5, 5), ylim=(-5, 5))

    # point_A, = ax.plot(r_Ai, 'b.')
    # point_B, = ax.plot(r_Bi, 'g.')
    # point_C, = ax.plot(r_Ci, 'r.')
    ligne_A, = ax.plot(r_Ai[0], r_Ai[1], 'b-', label="Corps A")
    ligne_B, = ax.plot(r_Bi[0], r_Bi[1], 'g-', label="Corps B")
    ligne_C, = ax.plot(r_Ci[0], r_Ci[1], 'r-', label="Corps C")
    titre = ax.set_title("Mouvement des trois corps à t= {}".format(0))

    # anim_point_A = lambda i: point_A.set_data(mouton["A"][i])
    # anim_point_B = lambda i: point_B.set_data(mouton["B"][i])
    # anim_point_C = lambda i: point_C.set_data(mouton["C"][i])
    anim_ligne_A = lambda i: ligne_A.set_data(mouton["A"][:i, 0], mouton["A"][:i, 1])
    anim_ligne_B = lambda i: ligne_B.set_data(mouton["B"][:i, 0], mouton["B"][:i, 1])
    anim_ligne_C = lambda i: ligne_C.set_data(mouton["C"][:i, 0], mouton["C"][:i, 1])
    anim_titre = lambda i: ax.set_title("Mouvement des trois corps\nà t= {}".format(round(mouton["t"][i], 3)))

    frames_anim = len(mouton["t"])
    # graph_anim_A = FuncAnimation(fig, anim_point_A, frames=frames_anim, interval=1)
    # graph_anim_B = FuncAnimation(fig, anim_point_B, frames=frames_anim, interval=1)
    # graph_anim_C = FuncAnimation(fig, anim_point_C, frames=frames_anim, interval=1)
    graph_anim_A = FuncAnimation(fig, anim_ligne_A, frames=frames_anim, interval=1)
    graph_anim_B = FuncAnimation(fig, anim_ligne_B, frames=frames_anim, interval=1)
    graph_anim_C = FuncAnimation(fig, anim_ligne_C, frames=frames_anim, interval=1)
    graph_anim_titre = FuncAnimation(fig, anim_titre, frames=frames_anim, interval=1)
    plt.legend()
    plt.show()


graph_3_corps(0, 1, 750)
