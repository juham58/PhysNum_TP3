import matplotlib.pyplot as plt
from No1_e import vdp, rangetemps
from scipy.integrate import solve_ivp
import scipy.constants as cst


if __name__ == "__main__":
    trange = [0, 8 * cst.pi]
    t = rangetemps(1000000, trange)

    solution1 = solve_ivp(vdp, [0, 8 * cst.pi], [0.5, 0], t_eval=t, method='RK45')
    sol1 = solution1.y

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol1[0], sol1[1], t)
    ax.set_zlabel("t")
    plt.title("Présentation 3D du système van der Pol qui évolue dans le temps avec x(0) = 0.5")
    plt.xlabel("x(t)")
    plt.ylabel("v(t)")

    plt.grid()
    plt.show()