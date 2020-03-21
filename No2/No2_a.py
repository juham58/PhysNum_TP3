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


def RK4_3_corps(t_i, t_f, N):
    r_A, r_B, r_C = r_Ai, r_Bi, r_Ci
    t_points = np.linspace(t_i, t_f, N)
    rA_arr = np.zeros((len(t_points), 2))
    rB_arr = np.zeros((len(t_points), 2))
    rC_arr = np.zeros((len(t_points), 2))
    h = (t_f-t_i)/N
    for i, t in enumerate(t_points):
        rA_arr[i][0], rA_arr[i][1] = r_A[0], r_A[1]
        rB_arr[i][0], rB_arr[i][1] = r_B[0], r_B[1]
        rC_arr[i][0], rC_arr[i][1] = r_C[0], r_C[1]
        k1_A, k1_B, k1_C = h*F("A", r_A, r_B, r_C), h*F("B", r_A, r_B, r_C), h*F("C", r_A, r_B, r_C)
        k2_A, k2_B, k2_C = h*F("A", r_A + 0.5*k1_A, r_B + 0.5*k1_B, r_C + 0.5*k1_C), h*F("B", r_A + 0.5*k1_A, r_B + 0.5*k1_B, r_C + 0.5*k1_C), h*F("B", r_A + 0.5*k1_A, r_B + 0.5*k1_B, r_C + 0.5*k1_C)
        k3_A, k3_B, k3_C = h*F("A", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C), h*F("B", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C), h*F("B", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C)
        k4_A, k4_B, k4_C = h*F("A", r_A+k3_A, r_B+k3_B, r_C+k3_C), h*F("B", r_A+k3_A, r_B+k3_B, r_C+k3_C), h*F("C", r_A+k3_A, r_B+k3_B, r_C+k3_C)
        r_A += (k1_A + 2*k2_A + 2*k3_A + k4_A)/6
        r_B += (k1_B + 2*k2_B + 2*k3_B + k4_B)/6
        r_C += (k1_C + 2*k2_C + 2*k3_C + k4_C)/6
    return {"A": rA_arr, "B": rB_arr, "C": rC_arr, "t": t_points}


#data_A = [r_Ai]
#data_B = [r_Bi]
#data_C = [r_Ci]



def graph_3_corps(t_i, t_f, N):
    RK4 = RK4_3_corps(t_i, t_f, N)
    fig, ax = plt.subplots()
    ax.set(xlim=(-25, 25), ylim=(-25, 25))
    point_A, = ax.plot(r_Ai, 'o')
    point_B, = ax.plot(r_Bi, 'o')
    point_C, = ax.plot(r_Ci, 'o')
    animation_A = lambda i: point_A.set_data(RK4["A"][i])
    animation_B = lambda i: point_B.set_data(RK4["B"][i])
    animation_C = lambda i: point_C.set_data(RK4["C"][i])
    frames_anim = len(RK4["t"])
    graph_anim_A = FuncAnimation(fig, animation_A, frames=frames_anim, interval=100)
    graph_anim_B = FuncAnimation(fig, animation_B, frames=frames_anim, interval=100)
    graph_anim_C = FuncAnimation(fig, animation_C, frames=frames_anim, interval=100)
    plt.show()

graph_3_corps(0, 1, 200)
