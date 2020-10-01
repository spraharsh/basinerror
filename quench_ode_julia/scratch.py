from pele.potentials import Harmonic
import numpy as np




a = np.zeros(3)
b = np.ones(3)
pot = Harmonic(a, 2)

print(pot.getEnergyGradientHessian(b))