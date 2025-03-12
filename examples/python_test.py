import navtools as nt
import numpy as np

# TODO: add more python tests

print(nt.__doc__)

C = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], order="F")
q = np.zeros(4, order="F")
nt.attitude.dcm2quat(q, C)

C2 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]], order="F")
q2 = nt.attitude.dcm2quat(C2)

nt.attitude.RotX(1.0)

print(f"C = \n{C} \nq = \n{q} \n")
print(f"C2 = \n{C2} \nq2 = \n{q2} \n")

print(f"WGS84_OMEGA_SKEW = \n{nt.WGS84_OMEGA_SKEW}\ndtype = {type(nt.WGS84_OMEGA_SKEW)}")
