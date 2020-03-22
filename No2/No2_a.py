import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# définitions des masses des corps
m_A = 3
m_B = 4
m_C = 5


# définition des conditions initiales
r_Ai = np.array([1.0, 3.0])
r_Bi = np.array([-2.0, -1.0])
r_Ci = np.array([1.0, -1.0])


def F(corps, r_A, r_B, r_C):
    if corps == "A":
        return -4*constants.pi**2*(m_B*((r_A-r_B)/(np.linalg.norm(r_A-r_B)**3))+m_C*((r_A-r_C)/(np.linalg.norm(r_A-r_C)**3)))

    if corps == "B":
        return -4*constants.pi**2*(m_A*((r_B-r_A)/(np.linalg.norm(r_B-r_A)**3))+m_C*((r_B-r_C)/(np.linalg.norm(r_B-r_C)**3)))

    if corps == "C":
        return -4*constants.pi**2*(m_A*((r_C-r_A)/(np.linalg.norm(r_C-r_A)**3))+m_B*((r_C-r_B)/(np.linalg.norm(r_C-r_B)**3)))


def mouton_3_corps(t_i, t_f, N):
    r_A_1, r_B_1, r_C_1 = r_Ai, r_Bi, r_Ci
    t_points = np.linspace(t_i, t_f, N)
    rA_arr = np.zeros((len(t_points), 2))
    rB_arr = np.zeros((len(t_points), 2))
    rC_arr = np.zeros((len(t_points), 2))
    h = (t_f-t_i)/N

    k1_A, k1_B, k1_C = h*F("A", r_Ai, r_Bi, r_Ci), h*F("B", r_Ai, r_Bi, r_Ci), h*F("C", r_Ai, r_Bi, r_Ci)
    k2_A, k2_B, k2_C = h*F("A", r_Ai + 0.5*k1_A, r_Bi + 0.5*k1_B, r_Ci + 0.5*k1_C), h*F("B", r_Ai + 0.5*k1_A, r_Bi + 0.5*k1_B, r_Ci + 0.5*k1_C), h*F("C", r_Ai + 0.5*k1_A, r_Bi + 0.5*k1_B, r_Ci + 0.5*k1_C)
    r_A_2, r_B_2, r_C_2 = r_Ai + k2_A, r_Bi + k2_B, r_Ci + k2_C

    r_A, r_B, r_C = r_A_1 + h*F("A", r_A_2, r_B_2, r_C_2), r_B_1 + h*F("B", r_A_2, r_B_2, r_C_2), r_C_1 + h*F("C", r_A_2, r_B_2, r_C_2)
    for i, t in enumerate(t_points):
        rA_arr[i][0], rA_arr[i][1] = r_A[0], r_A[1]
        rB_arr[i][0], rB_arr[i][1] = r_B[0], r_B[1]
        rC_arr[i][0], rC_arr[i][1] = r_C[0], r_C[1]

        # premier calcul pour un demi-saut
        r_A_1, r_B_1, r_C_1 = r_A_2, r_B_2, r_C_2
        r_A_2, r_B_2, r_C_2 = r_A, r_B, r_C
        r_A, r_B, r_C = r_A_1 + h*F("A", r_A_2, r_B_2, r_C_2), r_B_1 + h*F("B", r_A_2, r_B_2, r_C_2), r_C_1 + h*F("C", r_A_2, r_B_2, r_C_2)

        # deuxième calcul pour compléter le saut
        r_A_1, r_B_1, r_C_1 = r_A_2, r_B_2, r_C_2
        r_A_2, r_B_2, r_C_2 = r_A, r_B, r_C
        r_A, r_B, r_C = r_A_1 + h*F("A", r_A_2, r_B_2, r_C_2), r_B_1 + h*F("B", r_A_2, r_B_2, r_C_2), r_C_1 + h*F("C", r_A_2, r_B_2, r_C_2)

    return {"A": rA_arr, "B": rB_arr, "C": rC_arr, "t": t_points}


def graph_3_corps(t_i, t_f, N):
    data_ligne_A = []
    data_ligne_B = []
    data_ligne_C = [r_Ci]

    mouton = mouton_3_corps(t_i, t_f, N)
    print(mouton["A"][:,0])
    fig, ax = plt.subplots()
    ax.set(xlim=(-5, 5), ylim=(-5, 5))

    #point_A, = ax.plot(r_Ai, 'b.')
    #point_B, = ax.plot(r_Bi, 'g.')
    #point_C, = ax.plot(r_Ci, 'r.')
    ligne_A, = ax.plot(r_Ai[0], r_Ai[1], 'b-')
    ligne_B, = ax.plot(r_Bi[0], r_Bi[1],'g-')
    ligne_C, = ax.plot(r_Ci[0], r_Ci[1], 'r-')

    #anim_point_A = lambda i: point_A.set_data(mouton["A"][i])
    #anim_point_B = lambda i: point_B.set_data(mouton["B"][i])
    #anim_point_C = lambda i: point_C.set_data(mouton["C"][i])
    anim_ligne_A = lambda i: ligne_A.set_data(mouton["A"][:i, 0], mouton["A"][:i, 1])
    anim_ligne_B = lambda i: ligne_B.set_data(mouton["B"][:i, 0], mouton["B"][:i, 1])
    anim_ligne_C = lambda i: ligne_C.set_data(mouton["C"][:i, 0], mouton["C"][:i, 1])

    frames_anim = len(mouton["t"])
    #graph_anim_A = FuncAnimation(fig, anim_point_A, frames=frames_anim, interval=1)
    #graph_anim_B = FuncAnimation(fig, anim_point_B, frames=frames_anim, interval=1)
    #graph_anim_C = FuncAnimation(fig, anim_point_C, frames=frames_anim, interval=1)
    graph_anim_A = FuncAnimation(fig, anim_ligne_A, frames=frames_anim, interval=1)
    graph_anim_B = FuncAnimation(fig, anim_ligne_B, frames=frames_anim, interval=1)
    graph_anim_C = FuncAnimation(fig, anim_ligne_C, frames=frames_anim, interval=1)
    plt.show()


graph_3_corps(0, 1, 8000)
