import numpy as np
import matplotlib.pyplot as plt


def plotMach(U, node, elem):
    X, Y = node.T
    # read the state matrix
    gamma = 1.4
    r = U[:, 0]
    u = U[:, 1] / r
    v = U[:, 2] / r
    q = np.sqrt(u ** 2 + v ** 2)
    p = (gamma - 1) * (U[:, 3] - 0.5 * r * q ** 2)
    c = np.sqrt(gamma * p / r)
    M = q / c

    f = plt.figure(figsize=(12, 6))
    plt.tripcolor(X, Y, triangles=elem, facecolors=M, shading='flat', edgecolor='black')
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.tick_params()
    plt.set_cmap('jet')
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('Converged Mach number', fontsize=16)
    f.tight_layout()
    plt.axis('equal')
    plt.show()
    plt.close(f)


# plot the cp and calculates the lift and drag coefficient
# INPUTS: U = state matrix (4xNe)
#         Uinf = free-stream state (1x4)
#         B2E = boundary face to element mapping matrix (Nbx3)
#         Bn = boundary face normal vector (Nbx2)
#         elem = nodes to elements mapping matrix (Nex3)
#         nodes = nodes coordinates matrix (number of nodes x 2)
#         alpha = AoA
# OUTPUTS: cd = drag coefficient
#          cl = lift coefficient
def outputs(U, Uinf, B2E, Bn, elem, nodes, alpha):
    U = np.loadtxt(U)
    Uinf = np.loadtxt(Uinf)
    B2E = np.loadtxt(B2E)
    Bn = np.loadtxt(Bn)
    elem = np.loadtxt(elem)
    nodes = np.loadtxt(nodes)

    # read the state matrix
    gamma = 1.4
    r = U[:, 0]
    u = U[:, 1] / r
    v = U[:, 2] / r
    q = np.sqrt(u ** 2 + v ** 2)
    p = (gamma - 1) * (U[:, 3] - 0.5 * r * q ** 2)

    # free-stream state
    c = 1
    rinf = Uinf[0]
    uinf = Uinf[1] / rinf
    vinf = Uinf[2] / rinf
    qinf = np.sqrt(uinf ** 2 + vinf ** 2)
    pinf = (gamma - 1) * (Uinf[:, 3] - 0.5 * rinf * qinf ** 2)

    Fx = 0
    Fy = 0
    pB = np.array([])
    x = np.array([])
    for i in range(len(B2E)):
        if B2E[i][2] != 4:
            ielem = B2E[i][0]
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
            l = np.sqrt((node2Coord[0] - node1Coord[0]) ** 2 + (node2Coord[1] - node1Coord[1]) ** 2)
            x = np.append(x, (node1Coord[0] + node2Coord[0]) / 2)
            pB = np.append(pB, p)
            
            nx = Bn[i][0]
            ny = Bn[i][1]
            Fx += l * p[ielem] * -nx
            Fy += l * p[ielem] * -ny
    D = np.cos(alpha) * Fx - np.sin(alpha) * Fy
    L = np.sin(alpha) * Fx + np.cos(alpha) * Fy

    cl = L / (1 / 2 * rinf * qinf ** 2 * c)
    cd = D / (1 / 2 * rinf * qinf ** 2 * c)
    cp = (pB - pinf) / (1/2*qinf**2*c)
    # cp plot
    f = plt.figure(figsize=(8, 8))
    plt.xlabel('x', fontsize=16)
    plt.ylabel(r'$c_p$', fontsize=16)
    plt.title('Pressure coefficient distribution', fontsize=16)
    f.tight_layout()
    plt.scatter(x, cp)
    return cd, cl