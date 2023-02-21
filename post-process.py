import numpy as np
import matplotlib.pyplot as plt
from matricesGenerator import getB2E
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

    f = plt.figure(figsize=(12, 6))
    plt.tripcolor(X, Y, triangles=elem, facecolors=p, shading='flat', edgecolor='black')
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
    pB = np.array([])
    x = np.array([])
    iElem = len(E)
    count = 0
    for k in range(1, len(B)):
        for i in range(len(B[k])):
            for j, ne in enumerate(np.isin(E, B[k][i])):
                if np.count_nonzero(ne) == 2:
                    print(ne, )
                    iElem = i
                    count += 1
                    break
            node1Coord = V[B[k][i][0] - 1]
            node2Coord = V[B[k][i][1] - 1]
            l = np.sqrt((node2Coord[0] - node1Coord[0]) ** 2 + (node2Coord[1] - node1Coord[1]) ** 2)
            x = np.append(x, (node1Coord[0] + node2Coord[0]) / 2)
            pB = np.append(pB, p[iElem])

            nx = (node2Coord[1] - node1Coord[1]) / l
            ny = -(node2Coord[0] - node1Coord[0]) / l
            Fx += l * p[iElem] * nx
            Fy += l * p[iElem] * ny
    print(count)
    D = np.cos(alpha) * Fx - np.sin(alpha) * Fy
    L = np.sin(alpha) * Fx + np.cos(alpha) * Fy

    cl = L / (1 / 2 * rinf * qinf ** 2 * c)
    cd = D / (1 / 2 * rinf * qinf ** 2 * c)
    print(pB, len(pB))
    cp = (pB - pinf) / (1 / 2 * rinf * qinf ** 2)
    return cd, cl, [x, cp]


def computeFreestreamState(Minf, alpha):
    gamma = 1.4
    uInf = np.array([1, Minf * np.cos(alpha), Minf * np.sin(alpha), (1 / (gamma ** 2 - gamma)) + 0.5 * (Minf ** 2)])
    return uInf


def plotCp(cp):
    # cp plot
    f = plt.figure(figsize=(8, 8))
    plt.xlabel('x', fontsize=16)
    plt.ylabel(r'$c_p$', fontsize=16)
    plt.title('Pressure coefficient distribution', fontsize=16)
    plt.scatter(cp[0], cp[1])
    f.tight_layout()
    plt.show()


def main():
    U = np.loadtxt('outputStates_2ndOrder_TRANSONIC_32k_jb_rusanov.txt')
    Mesh = 'c2.gri'
    Minf = 0.5
    alpha = 8 * np.pi / 180
    uInf = computeFreestreamState(Minf, alpha)
    cd, cl, cp = outputs(U, uInf, Mesh, alpha)
    print(f"cd = {cd}, cl = {cl}")
    plotCp(cp)
    plotMach(U, Mesh, ['8k', '2nd-order'])


if __name__ == "__main__":
    main()
