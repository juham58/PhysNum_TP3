import numpy as np
from scipy import constants
import matplotlib.pyplot as plt


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
    rA_liste = []
    rB_liste = []
    rC_liste = []
    h = (t_f-t_i)/N
    for t in t_points:
        rA_liste.append(r_A), rB_liste.append(r_B), rC_liste.append(r_C)
        k1_A, k1_B, k1_C = h*F("A", r_A, r_B, r_C), h*F("B", r_A, r_B, r_C), h*F("C", r_A, r_B, r_C)
        k2_A, k2_B, k2_C = h*F("A", r_A+0.5*k1_A, r_B+0.5*k1_B, r_C+0.5*k1_C), h*F("B", r_A+0.5*k1_A, r_B+0.5*k1_B, r_C+0.5*k1_C), h*F("B", r_A+0.5*k1_A, r_B+0.5*k1_B, r_C+0.5*k1_C)
        k3_A, k3_B, k3_C = h*F("A", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C), h*F("B", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C), h*F("B", r_A+0.5*k2_A, r_B+0.5*k2_B, r_C+0.5*k2_C)
        k4_A, k4_B, k4_C = h*F("A", r_A+k3_A, r_B+k3_B, r_C+k3_C), h*F("B", r_A+k3_A, r_B+k3_B, r_C+k3_C), h*F("C", r_A+k3_A, r_B+k3_B, r_C+k3_C)
        r_A += (k1_A + 2*k2_A + 2*k3_A + k4_A)/6
        r_B += (k1_B + 2*k2_B + 2*k3_B + k4_B)/6
        r_C += (k1_C + 2*k2_C + 2*k3_C + k4_C)/6
        print(r_A)
    return {"A": rA_liste, "B": rB_liste, "C": rC_liste, "t": t_points}

def graph_3_corps(t_i, t_f, N):
    RK4 = RK4_3_corps(t_i, t_f, N)
    plt.figure(figsize=(16,8))
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    #plt.plot(RK4["A"][0], RK4["A"][1])
    #plt.plot(RK4["B"][0], RK4["B"][1])
    plt.plot(RK4["C"][0], RK4["C"][1])
    plt.show()

graph_3_corps(0, 1, 100)
