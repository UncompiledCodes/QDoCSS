from qutip import *
from numpy import pi

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


print(gama)
