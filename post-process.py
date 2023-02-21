import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri


def plotMach(U, Mesh, name):
    elem = readgri(Mesh)['E']
    nodes = readgri(Mesh)['V']
    X, Y = nodes.T
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
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel('y (m)', fontsize=16)
    f.tight_layout()
    plt.xlim([-0.3, 1.3])
    plt.ylim([-0.4, 0.4])
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
def outputs(U, Uinf, Mesh, alpha):
    # B2E = getB2E(Mesh)
    E = readgri(Mesh)['E']
    V = readgri(Mesh)['V']
    B = readgri(Mesh)['B']
    IE = readgri()
    Bname = readgri(Mesh)['Bname']

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
    pinf = (gamma - 1) * (Uinf[3] - 0.5 * rinf * qinf ** 2)

    Fx = 0
    Fy = 0
    pBslat = np.array([])
    pBmain = np.array([])
    pBflap = np.array([])
    xslat = np.array([])
    xmain = np.array([])
    xflap = np.array([])
    x = np.array([])
    iElem = len(E)
    for k in range(len(B) - 1):
        for i in range(len(B[k])):
            for j, ne in enumerate(np.isin(E, B[k][i])):
                if np.count_nonzero(ne) == 2:
                    iElem = j
                    break
            node1Coord = V[B[k][i][0] - 1]
            node2Coord = V[B[k][i][1] - 1]
            l = np.sqrt((node2Coord[0] - node1Coord[0]) ** 2 + (node2Coord[1] - node1Coord[1]) ** 2)
            xi = np.append(x, (node1Coord[0] + node2Coord[0]) / 2)
            if k == 1:
                pBslat = np.append(pBslat, p[iElem])
                xslat = np.append(xslat, xi)
            elif k == 2:
                pBmain = np.append(pBmain, p[iElem])
                xmain = np.append(xmain, xi)
            else:
                pBflap = np.append(pBflap, p[iElem])
                xflap = np.append(xflap, xi)
            nx = (node2Coord[1] - node1Coord[1]) / l
            ny = -(node2Coord[0] - node1Coord[0]) / l
            Fx += l * p[iElem] * nx
            Fy += l * p[iElem] * ny

    D = np.cos(alpha) * Fx + np.sin(alpha) * Fy
    L = -np.sin(alpha) * Fx + np.cos(alpha) * Fy

    pB = np.array([pBslat, pBmain, pBflap])
    x = np.array([xslat, xmain, xflap])
    cl = L / (1 / 2 * rinf * qinf ** 2 * c)
    cd = D / (1 / 2 * rinf * qinf ** 2 * c)
    cp = (pB - pinf) / (1 / 2 * rinf * qinf ** 2)
    return cd, cl, x, cp


def computeFreestreamState(Minf, alpha):
    gamma = 1.4
    uInf = np.array([1, Minf * np.cos(alpha), Minf * np.sin(alpha), (1 / (gamma ** 2 - gamma)) + 0.5 * (Minf ** 2)])
    return uInf


def plotCp(x, cp):
    # cp plot
    f = plt.figure(figsize=(8, 8))
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel(r'$c_p$', fontsize=16)
    plt.scatter(x[0], cp[0], label='flap')
    #plt.scatter(x[1], cp[1], label='main')
    #plt.scatter(x[2], cp[2], label='slat')
    f.tight_layout()
    plt.legend()
    plt.show()


def main():
    U = np.loadtxt('outputStates_adapt2.txt')
    Mesh = 'adapt2.gri'
    Minf = 0.5
    alpha = 8 * np.pi / 180
    uInf = computeFreestreamState(Minf, alpha)
    cd, cl, x, cp = outputs(U, uInf, Mesh, alpha)
    print(f"cd = {cd}, cl = {cl}")
    plotCp(x, cp)
    plotMach(U, Mesh, ['8k', '2nd-order'])


if __name__ == "__main__":
    main()
