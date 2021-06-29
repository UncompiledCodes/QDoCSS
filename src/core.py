from qutip import *
from numpy import pi
from scipy.special import jv

############## Constants ##################
si = qeye(2)
sx = sigmax()
sy = sigmay()
sz = sigmaz()
n = 2  # number of nuclei
m = 1  # number of electrons
extmagfield_m = 1
extmagfield_n = 1  # Qobj([[0, 0, 0], [0, 0, 0], [0.3, 0.3, 0.3]])
JJ = 1  # Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
AA = 1  # Qobj([[0, 0, 0], [0, 0, 0], [0.99925, 0.625147, 0.551273]])
# gama = Qobj(
#     [
#         [0.0124425, 0.0806628, 0.00999575],
#         [0.0550028, 0.0758354, 0.07346340],
#         [0.0972069, 0.0723954, 0.07405450],
#     ]
# )
gama = 1
t = (10 * 2 * pi) / extmagfield_m
E1 = (
    1 / 2 * JJ
    + 1 / 2 * gama
    + 1 / 2 * AA
    + 1 / 2 * extmagfield_n
    + 1 / 2 * extmagfield_m
)
tau = E1 * t / 2
kappa = 3 / 2 * tau
############## Constants ##################

############## Spin #######################
def spin_op(N):
    """calculates spin oprators

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


############## Spin #######################


############## Hamiltonian ################
def hamiltonian(m, n, JJ, AA, gama, extmagfield_m, extmagfield_n):
    """calculates spin oprators and
       using them, calculates system hamiltonian

    Args:
        m (int): number of electrons
        n (int): number of nuclei
        JJ (matrix): exchange interaction constant
        AA (matrix): exchange interaction constant
        gama (matrix): exchange interaction constant
        extmagfield_m (int): external magnetic field of m
        extmagfield_n (matrix): external magnetic field of n

    Returns:
        matrix: hamiltonian of the entire system
    """
    # calculating spin operators
    Sx, Sy, Sz = spin_op(m + n)
    # HS & HB are bare hamiltonians of the central system and bath
    # V is the system-bath interaction
    HS = HB = V = 0
    for i in range(n, m + 1):
        HS += extmagfield_m * Sx[i] + extmagfield_m * Sy[i] + extmagfield_m * Sz[i]
        for j in range(n, m + 1):
            if m != 1:
                if i == j:
                    continue
            HS += JJ * Sx[i] * Sx[j] + JJ * Sy[i] * Sy[j] + JJ * Sz[i] * Sz[j]

    for i in range(n):
        HB += extmagfield_n * Sx[i] + extmagfield_n * Sy[i] + extmagfield_n * Sz[i]
        for j in range(n):
            if n != 1:
                if i == j:
                    continue
            HB += gama * Sx[i] * Sx[j] + gama * Sy[i] * Sy[j] + gama * Sz[i] * Sz[j]

    for i in range(n, m + n):
        for j in range(n):
            V += AA * Sx[i] * Sx[j] + AA * Sy[i] * Sy[j] + AA * Sz[i] * Sz[j]

    HH = HS + HB + V
    return HH


############## Hamiltonian ################

############ Chebyshev series #############
def TG(k, G):
    """recursive function to calculate T

    Args:
        k (int): number of times to calculate
        G (matrix): oprator G

    Returns:
        matrix: matrix of T in a speccific k
    """
    if k == 1:
        return 1
    elif k == 2:
        return G
    else:
        return 2 * G * TG(k - 1, G) - TG(k - 2, G)


############ Chebyshev series #############


############ Evolution Operator ###########
def Ut(kappa, tau, G):
    """calculating the evolution oprator

    Args:
        kappa ([type]): [description]
        tau ([type]): [description]
        G ([type]): [description]

    Returns:
        [type]: [description]
    """
    UU = 0
    for k in range(1, int(kappa) + 1):
        a = 1
        if k == 0:
            a = 2
        UU += a * ((1j) ** k) * jv(k, tau) * TG(k, G)
    return UU


############ Evolution Operator ###########


def main(m, n, JJ, AA, gama, extmagfield_m, extmagfield_n, kappa, tau):
    HH = hamiltonian(m, n, JJ, AA, gama, extmagfield_m, extmagfield_n)
    G = 2 * HH / E1
    UU = Ut(kappa, tau, G)
    return Qobj(UU)


if __name__ == "__main__":
    main(m, n, JJ, AA, gama, extmagfield_m, extmagfield_n, kappa, tau)
