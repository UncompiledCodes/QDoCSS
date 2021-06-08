from qutip import *
from numpy import pi


si = qeye(2)
sx = sigmax()
sy = sigmay()
sz = sigmaz()
sm = sigmam()
sp = sigmap()
n = 3  # number of nuclei

m = 1  # number of electrons

extmagfield_m = 1  # external magnetic field of m

extmagfield_n = Qobj(
    [[0, 0, 0], [0, 0, 0], [0.3, 0.3, 0.3]]
)  # external magnetic field of n

j = Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 0]])  # value of exchange interaction constant: J

a = Qobj(
    [[0, 0, 0], [0, 0, 0], [0.99925, 0.625147, 0.551273]]
)  # value of exchange interaction constant: A

gama = Qobj(
    [
        [0.0124425, 0.0806628, 0.00999575],
        [0.0550028, 0.0758354, 0.07346340],
        [0.0972069, 0.0723954, 0.07405450],
    ]
)  # value of exchange interaction constant: Gama

t = (10 * 2 * pi) / extmagfield_m  # value of time

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
    Sp = []
    Sm = []

    for n in range(N):
        x_list = []
        y_list = []
        z_list = []
        p_list = []
        m_list = []
        for m in range(N):
            x_list.append(si)
            y_list.append(si)
            z_list.append(si)
            p_list.append(si)
            m_list.append(si)

        x_list[n] = sx
        y_list[n] = sy
        z_list[n] = sz
        p_list[n] = sp
        m_list[n] = sm

        Sx.append(tensor(x_list))
        Sy.append(tensor(y_list))
        Sz.append(tensor(z_list))
        Sp.append(tensor(p_list))
        Sm.append(tensor(m_list))
    # print(Qobj(Sx[1]))
    return [Sx, Sy, Sz, Sp, Sm]
    # return Qobj(Sx)


Sx, Sy, Sz, Sp, Sm = spin_operators(2)

for i in range(2):
    print(Sx[i])
