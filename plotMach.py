import numpy as np
import matplotlib.pyplot as plt


def plotMach(node, elem, M, blade):
    X, Y = node.T

    f = plt.figure(figsize=(12, 6))
    plt.tripcolor(X, Y, elem - 1, facecolors=M, shading='flat', edgecolor='black')
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.tick_params()
    plt.set_cmap('jet')
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('Converged Mach number for ' + blade, fontsize=16)
    f.tight_layout()
    plt.axis('equal')
    plt.show()
    plt.close(f)


def plotCp(Uinf, U, B2E, Bn, elem, nodes, alpha):
    U = np.loadtxt(U)
    B2E = np.loadtxt(B2E)
    Bn = np.loadtxt(Bn)
    elem = np.loadtxt(elem)
    nodes = np.loadtxt(nodes)
    Fx = 0
    Fy = 0
    for i in range(len(B2E)):
        if B2E[i][2] != 4:
            ielem = B2E[i][0]
            r = U[ielem][0]
            q = np.sqrt(U[ielem][1] ** 2 + U[ielem][2] ** 2) / r
            p = (gamma - 1) * (U[3] - 0.5 * r * q ** 2)

            face = B2E[i][1]
            if face == 1:
                node1 = elem[i][1]
                node2 = elem[i][2]
            elif face == 2:
                node1 = elem[i][0]
                node2 = elem[i][2]
            else:
                node1 = elem[i][0]
                node2 = elem[i][1]

            node1Coord = nodes[node1]
            node2Coord = nodes[node2]
            l = np.sqrt((node2Coord[0]-node1Coord[0])**2 + (node2Coord[1]-node1Coord[1])**2)

            nx = Bn[i][0]
            ny = Bn[i][1]
            Fx += l*p * -nx
            Fy += l * p * -ny

    D = np.cos(alpha)*Fx - np.sin(alpha)*Fy
    L = np.sin(alpha)*Fx + np.cos(alpha)*Fy
    rinf = Uinf[0]
    uinf = Uinf[1]/rinf
    vinf = Uinf[2]/rinf
    Umag = np.sqrt(uinf**2 + vinf**2)
    c = 1
    cl = L/(1/2 * rinf * Umag**2 * c)
    cd = D / (1 / 2 * rinf * Umag ** 2 * c)

