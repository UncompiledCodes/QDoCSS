from qutip import *
from numpy import pi

# Identity matrix
si = qeye(2)
# pauli x
sx = sigmax()
# pauli y
sy = sigmay()
# pauli z
sz = sigmaz()
# number of nuclei
n = 3
# number of electrons
m = 1
# external magnetic field of m
extmagfield_m = 1
# external magnetic field of n
extmagfield_n = Qobj([[0, 0, 0], [0, 0, 0], [0.3, 0.3, 0.3]])
# value of exchange interaction constant: J
JJ = Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
# value of exchange interaction constant: A
AA = Qobj([[0, 0, 0], [0, 0, 0], [0.99925, 0.625147, 0.551273]])
# value of exchange interaction constant: Gama
gama = Qobj(
    [
        [0.0124425, 0.0806628, 0.00999575],
        [0.0550028, 0.0758354, 0.07346340],
        [0.0972069, 0.0723954, 0.07405450],
    ]
)
# value of time
t = (10 * 2 * pi) / extmagfield_m

# spin oprators
def spin_operators(N):
    """calculating spin oprators

    Args:
        N (int): number of electrons or atoms

    Returns:
        matrix: matrices of spin oprators
    """
    Sx = []
    Sy = []
    Sz = []

    for n in range(N):
        x_list = []
        y_list = []
        z_list = []

        for m in range(N):
            x_list.append(si)
            y_list.append(si)
            z_list.append(si)

        x_list[n] = sx
        y_list[n] = sy
        z_list[n] = sz

        Sx.append(tensor(x_list))
        Sy.append(tensor(y_list))
        Sz.append(tensor(z_list))

    return [Sx, Sy, Sz]


def hamiltonian(m, n, JJ, AA, gama, extmagfield_m, extmagfield_n):
    """calculating system hamiltonian

    Args:
        m (int): number of electrons
        n (int): number of nuclei
        JJ (matrix): exchange interaction constant
        AA (matrix): exchange interaction constant
        gama (matrix): exchange interaction constant
        extmagfield_m (int): external magnetic field of m
        extmagfield_n (matrix): external magnetic field of n

    Returns:
        matrix: hamiltonian of he entire system
    """
    # calculating spin operators
    Sx, Sy, Sz = spin_operators(m)
    Ix, Iy, Iz = spin_operators(n)
    # HS & HB are bare hamiltonians of the central system and bath
    # V is the system-bath interaction
    HS = HB = V = 0
    for i in range(m):
        HS += extmagfield_m * Sx[i] + extmagfield_m * Sy[i] + extmagfield_m * Sz[i]
        for j in range(m):
            if m != 1:
                if i == j:
                    continue
            HS += JJ * Sx[i] * Sx[j] + JJ * Sy[i] * Sy[j] + JJ * Sz[i] * Sz[j]

    for i in range(n):
        HB += extmagfield_n * Ix[i] + extmagfield_n * Iy[i] + extmagfield_n * Iz[i]
        for j in range(n):
            if n != 1:
                if i == j:
                    continue
            HB += gama * Ix[i] * Ix[j] + gama * Iy[i] * Iy[j] + gama * Iz[i] * Iz[j]

    for i in range(m):
        for j in range(n):
            V += AA * Sx[i] * Ix[j] + AA * Sy[i] * Iy[j] + AA * Sz[i] * Iz[j]

    HH = HS + HB + V
    return HH
