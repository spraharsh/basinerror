"""
Writes potential to file to compare against julia
"""
from pele.potentials import InversePower
import numpy as np


test_radii = [1.08820262, 1.02000786, 1.0489369,  1.11204466, 1.53072906, 1.33159055,  1.46650619, 1.389405]

test_coords = [6.43109269, 2.55893249,
               5.28365038,3.52962924,
               3.79089799, 6.1770549,
               0.47406571, 0.58146545,
               0.13492935, 5.55656566,
               5.19310115, 5.80610666,
               6.53090015, 5.3332587,
               3.07972527, 5.20893375]
test_radii = np.array(test_radii)
test_coords = np.array(test_coords)
box_length = 6.673592625078725
potential = InversePower(2.5,
                         1.,
                         use_cell_lists=False,
                         ndim=2.0,
                         radii=test_radii,
                         boxvec=[box_length, box_length])

np.savetxt('energy.csv', np.array([potential.getEnergy(test_coords)]), delimiter=',')
np.savetxt('gradient.csv', potential.getEnergyGradient(test_coords)[1], delimiter=',')
np.savetxt('hessian.csv', potential.getEnergyGradientHessian(test_coords)[2], delimiter=',')





print(potential.getEnergy(test_coords))
print(potential.getEnergyGradient(test_coords)[1])
print(potential.getEnergyGradientHessian(test_coords)[2])
